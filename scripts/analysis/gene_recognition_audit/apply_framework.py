"""Mechanical application of DECISION_FRAMEWORK_FROZEN.md. Identity only (quantity 1). No threshold edited.
Two estimators per gene-cancer pair:
 (A) random/all-matched stratum estimate (unbiased sample; targeted oversamples excluded), with Wilson 95% CI;
 (B) confirmed-false count in complete strata / pool credits = exhaustive LOWER BOUND for the mechanisms found (every abstract of the alias stratum was read
     or rule-checked pool-wide); extras = nested-name scan hits that were not in the read set. Uncertain (U) cases reported separately, never counted as wrong."""
import pandas as pd, numpy as np, math, glob
K=[("breast","Breast Cancer"),("lung","Lung Cancer"),("colorectal","Colorectal Cancer"),("prostate","Prostate Cancer"),("melanoma","Melanoma"),("bladder","Bladder Cancer"),("kidney","Kidney Cancer"),("pancreas","Pancreatic Cancer"),("liver","Liver Cancer"),("stomach","Stomach Cancer"),("esophageal","Esophageal Cancer"),("ovarian","Ovarian Cancer"),("endometrial","Endometrial Cancer"),("cervical","Cervical Cancer"),("thyroid","Thyroid Cancer"),("brain","Brain Cancer"),("leukemia","Leukemia"),("lymphoma","Lymphoma"),("myeloma","Myeloma")]
F=pd.concat([pd.read_csv(f"{k}_qc_final.csv",dtype=str).fillna("").assign(key=k) for k,_ in K],ignore_index=True)
P=pd.concat([pd.read_csv(f"{k}_pool_gene_matches.csv",dtype=str).assign(cancer=l) for k,l in K],ignore_index=True)
cred=P.groupby(["cancer","gene"]).size().rename("credits").reset_index()
def wilson(x,n,z=1.96):
    if n==0: return (np.nan,np.nan)
    p=x/n; d=1+z*z/n; c=p+z*z/(2*n); h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n)); return ((c-h)/d,(c+h)/d)
F["stratum"]=F.sampling_method.map(lambda m:"all" if m=="all" else ("targeted-only" if m=="targeted" else ("v1-top15" if "top15" in m else ("v1-alias" if "alias" in m else "random"))))
rnd=F[F.stratum.isin(["random","all"])]
rows=[]
nested_extra={("Breast Cancer","BRCA2"):3,("Prostate Cancer","BRCA1"):1}   # from nested_scan.csv, pool-wide, abstracts not in the read set
for (c,g),n_cred in cred.set_index(["cancer","gene"]).credits.items():
    a=rnd[(rnd.cancer==c)&(rnd.gene==g)]; x=int((a.identity=="X").sum()); n=len(a); lo,hi=wilson(x,n)
    allr=F[(F.cancer==c)&(F.gene==g)]; X_known=int(allr[allr.identity=="X"].pmid.nunique())+nested_extra.get((c,g),0); U_known=int(allr[allr.identity=="U"].pmid.nunique())
    rows.append(dict(cancer=c,gene=g,credits=n_cred,rand_n=n,rand_X=x,rand_pct=round(100*x/n,1) if n else np.nan,ci_lo=round(100*lo,1),ci_hi=round(100*hi,1),
        confirmed_X=X_known,confirmed_pct_lowerbound=round(100*X_known/n_cred,1),uncertain_U=U_known,U_pct=round(100*U_known/n_cred,1)))
D=pd.DataFrame(rows); D.to_csv("ALL19_framework_pair_rates.csv",index=False)
n_pairs=len(D)
pool_rand=100*rnd.identity.eq("X").sum()/len(rnd); pool_rand_U=100*rnd.identity.eq("U").sum()/len(rnd)
tot_cred=D.credits.sum(); pool_conf=100*D.confirmed_X.sum()/tot_cred
print(f"pairs={n_pairs}  total pool credits={tot_cred}  abstracts read={len(F)}  random/all sample n={len(rnd)}")
print(f"POOLED incorrect identity: random/all sample = {pool_rand:.2f}% X (+{pool_rand_U:.2f}% uncertain);  credit-weighted confirmed lower bound = {pool_conf:.2f}% ({D.confirmed_X.sum()} of {tot_cred})")
print(f"pair-weighted mean (random est.) = {D.rand_pct.mean():.2f}% ; pair-weighted mean (confirmed lower bound) = {D.confirmed_pct_lowerbound.mean():.2f}%")
print("pairs with confirmed_X>0:",int((D.confirmed_X>0).sum()),"of",n_pairs)
for t in (0,2,5,10): print(f"  pairs with confirmed lower-bound rate >{t}%: {(D.confirmed_pct_lowerbound>t).sum()}")
print("distribution of confirmed lower-bound rate across 190 pairs (quantiles):",D.confirmed_pct_lowerbound.quantile([.5,.75,.9,.95,1]).round(1).to_dict())
fl=D[(D.confirmed_pct_lowerbound>10)|(D.U_pct>10)].sort_values("confirmed_pct_lowerbound",ascending=False)
print("\nFLAGGED by framework row 4 (>10% incorrect, or >10% uncertain identity):"); print(fl[["cancer","gene","credits","confirmed_X","confirmed_pct_lowerbound","uncertain_U","U_pct","rand_n","rand_X"]].to_string(index=False))
inv=D[(D.confirmed_pct_lowerbound>5)&(D.confirmed_pct_lowerbound<=10)].sort_values("confirmed_pct_lowerbound",ascending=False)
print("\n5-10% tier (framework row 3: investigate):"); print(inv[["cancer","gene","credits","confirmed_X","confirmed_pct_lowerbound"]].to_string(index=False))
mid=D[(D.confirmed_pct_lowerbound>2)&(D.confirmed_pct_lowerbound<=5)]
print("\n2-5% tier (framework row 2):"); print(mid[["cancer","gene","credits","confirmed_X","confirmed_pct_lowerbound"]].to_string(index=False))
