from pathlib import Path

# --- CRITICAL FIX FOR SRC LAYOUT ---
# Old: pomaar/config.py          -> parents[1] = root
# New: src/pomaar/config.py      -> parents[2] = root
PROJ_ROOT = Path(__file__).resolve().parents[2]

# Standard Data Science structure
DATA_DIR = PROJ_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
INTERIM_DATA_DIR = DATA_DIR / "interim"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
EXTERNAL_DATA_DIR = DATA_DIR / "external"

# Reports
REPORTS_DIR = PROJ_ROOT / "reports"
FIGURES_DIR = REPORTS_DIR / "figures"

# If you ever need to debug paths:
if __name__ == "__main__":
    print(f"Project root is: {PROJ_ROOT}")
