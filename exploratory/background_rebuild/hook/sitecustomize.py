import os
import sys

if os.environ.get("BGFIX_ENABLE") == "1":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import bgfix_hook
    bgfix_hook.install()
