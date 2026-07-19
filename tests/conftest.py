import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

CATALOG_PATH = ROOT / "catalogs" / "somatiq_grch38.json"
