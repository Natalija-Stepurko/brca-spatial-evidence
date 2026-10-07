"""Fixed study parameters. Every value here is stated in docs/DESIGN.md; tests check that they agree."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RESULTS = ROOT / "results"

# §3 data and candidates
STUDY2_REPO = "https://github.com/Natalija-Stepurko/brca-target-evidence"
STUDY2_COMMIT = None                 # pinned by stage `data` when the tables are copied in
SUBTYPES = ("basal", "her2", "luminal")
SUBTYPE_TO_SECTION = {"basal": "TNBC", "her2": "HER2+", "luminal": "ER+"}

# §2 control genes
CONTROLS = {"tumour": ("EPCAM", "KRT8", "KRT18"), "immune": ("PTPRC",), "stroma": ("COL1A1", "DCN"),
            "housekeeping": ("ACTB", "GAPDH")}

# §4 compartment score
COMPARTMENT_LFC = 0.5                # tumour log-fold-change for a "tumour" call
N_PERMUTATIONS = 1000
N_MATCHED_RANDOM = 1000
MATCH_STRATA = {"expression": 5, "detection": 3}

# §6 histology
TILE_PX = 224
TILE_MPP = 0.5
ENCODERS = {"phikon": "hf-hub:owkin/phikon", "imagenet": "vit_base_patch16_224.augreg_in21k_ft_in1k"}
PCA_DIM = 256
N_TOP_VARIABLE = 50

# §7 statistics
SEOI_LFC = 0.3
SEOI_PEARSON = 0.05
SEOI_ODDS_RATIO = 1.5
N_BOOTSTRAP = 1000
SEED_LADDER = 0
SEED_SPLITS = 1
