import sys
from pathlib import Path

# test the source tree, not whatever copy of arcsv happens to be installed
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
