"""The numbers in src/bse/config.py must be the ones docs/DESIGN.md pre-registers."""
import re
from pathlib import Path

from bse import config as C
from bse.cli import STAGES

DESIGN = (Path(__file__).resolve().parents[1] / "docs" / "DESIGN.md").read_text()


def test_eight_predictions_are_pre_registered():
    assert re.findall(r"^\| \*\*P(\d)\*\*", DESIGN, flags=re.M) == [str(i) for i in range(1, 9)]


def test_thresholds_match():
    assert f"TLFC ≥ {C.COMPARTMENT_LFC}" in DESIGN
    assert f"{C.SEOI_LFC} in TLFC" in DESIGN
    assert f"{C.SEOI_PEARSON} in Pearson" in DESIGN
    assert f"odds ratio of {C.SEOI_ODDS_RATIO}" in DESIGN
    assert f"PCA to {C.PCA_DIM} dimensions" in DESIGN
    assert f"{C.TILE_PX} × {C.TILE_PX} px tile at {C.TILE_MPP} µm/px" in DESIGN


def test_counts_and_seeds():
    assert f"{C.N_PERMUTATIONS:,} times" in DESIGN
    assert f"{C.N_MATCHED_RANDOM:,} random gene sets" in DESIGN
    assert f"{C.N_BOOTSTRAP:,}-resample bootstrap" in DESIGN
    assert f"Seeds fixed ({C.SEED_LADDER} for the ladder, {C.SEED_SPLITS} for the histology splits)" in DESIGN


def test_controls_named_in_design():
    for genes in C.CONTROLS.values():
        for g in genes:
            assert g in DESIGN, g


def test_stages_named_in_design():
    for stage, _ in STAGES:
        assert f"`{stage}`" in DESIGN, stage


def test_no_banned_phrases():
    for path in ("docs/DESIGN.md", "README.md", "research/literature.md"):
        p = Path(__file__).resolve().parents[1] / path
        if p.exists():
            assert not re.search(r"rather than|instead of", p.read_text(), flags=re.I), path
