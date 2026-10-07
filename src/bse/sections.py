"""Load every Visium section into one shape: an AnnData with raw counts, spot positions, a compartment
label per spot (tumour / stroma / immune / other) and section metadata (dataset, patient, subtype,
full-resolution image path or None).

Wu 2021: compartment from the shipped pathologist classification. Li 2025: the deposit carries no
per-spot tissue-annotation table, so the compartment comes from a fixed marker panel that excludes the
control genes and the candidates (DESIGN §11). Janesick: the serial Visium section, marker panel too.
"""
import gzip
import json
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp

from bse import config as C

UNPACKED = C.DATA / "unpacked"
RAW = C.DATA / "raw"

# Wu's pathologist classes -> compartment
WU_SUBTYPE = {"TNBC": "TNBC", "ER": "ER+", "HER2": "HER2+"}
# markers for the expression-based labels (Li, Janesick); none is a control gene or a study-2 candidate
MARKERS = {"tumour": ["KRT19", "KRT7", "CDH1", "CLDN4", "ERBB3"],
           "immune": ["CD3E", "CD2", "CD79A", "MS4A1", "CD68", "CD14"],
           "stroma": ["COL1A2", "COL3A1", "LUM", "PDGFRB", "SPARC"]}


@dataclass
class Section:
    id: str
    dataset: str
    patient: str
    subtype: str              # TNBC / HER2+ / ER+
    adata: ad.AnnData         # .obs: compartment, array_row, array_col, pxl_row, pxl_col ; .X raw counts
    image: Path | None        # full-resolution H&E, or None
    label_source: str         # "pathologist" or "markers"


def wu_compartment(cls: str) -> str:
    c = str(cls).lower() if isinstance(cls, str) else ""
    if "normal" in c:                      # normal glands express epithelial genes; not a comparator
        return "other"
    if "cancer" in c or "dcis" in c:
        return "tumour"
    if c in ("lymphocytes", "tls"):
        return "immune"
    if "stroma" in c or "adipose" in c:
        return "stroma"
    return "other"


def _open(path: Path):
    """Some deposits name plain-text files .gz; open by content, not by suffix."""
    with path.open("rb") as f:
        magic = f.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else path.open("rt")


def read_mtx_dir(d: Path) -> ad.AnnData:
    m = scipy.io.mmread(_open(d / "matrix.mtx.gz")).tocsr().T.tocsr()
    feats = pd.read_csv(_open(d / "features.tsv.gz"), sep="\t", header=None)
    bcs = pd.read_csv(_open(d / "barcodes.tsv.gz"), sep="\t", header=None)[0].values
    a = ad.AnnData(X=sp.csr_matrix(m, dtype=np.float32), obs=pd.DataFrame(index=bcs),
                   var=pd.DataFrame(index=feats[feats.shape[1] - 1 if feats.shape[1] == 1 else 1].values))
    a.var_names_make_unique()
    return a


def read_positions(spatial_dir: Path) -> pd.DataFrame:
    p = spatial_dir / "tissue_positions_list.csv"
    if not p.exists():
        p = spatial_dir / "tissue_positions.csv"
    df = pd.read_csv(p, header=None if "list" in p.name else 0)
    df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    return df.set_index("barcode")


def marker_compartment(a: ad.AnnData) -> pd.Series:
    """Score each compartment's marker panel on log-normalised counts; label = highest score if it
    exceeds the others by 0.5, otherwise 'other'. Spots with a tumour score close to a second score are
    mixtures and go to 'other' too."""
    b = a.copy()
    sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b)
    scores = {}
    for k, genes in MARKERS.items():
        g = [x for x in genes if x in b.var_names]
        scores[k] = np.asarray(b[:, g].X.mean(axis=1)).ravel()
    S = pd.DataFrame(scores, index=a.obs_names)
    S = (S - S.mean()) / S.std()
    top = S.idxmax(axis=1)
    sorted_vals = np.sort(S.values, axis=1)
    margin = sorted_vals[:, -1] - sorted_vals[:, -2]
    return pd.Series(np.where(margin >= 0.5, top, "other"), index=a.obs_names)


# ----------------------------------------------------------------------------- loaders
def load_wu() -> list[Section]:
    base = UNPACKED / "wu"
    out = []
    for meta in sorted((base / "metadata").glob("*_metadata.csv")):
        sid = meta.name.replace("_metadata.csv", "")
        md = pd.read_csv(meta, index_col=0)
        a = read_mtx_dir(base / "filtered_count_matrices" / f"{sid}_filtered_count_matrix")
        pos = read_positions(base / "spatial" / f"{sid}_spatial")
        a = a[a.obs_names.isin(md.index)].copy()
        a.obs["compartment"] = md.loc[a.obs_names, "Classification"].map(wu_compartment).values
        a.obs["pathologist_class"] = md.loc[a.obs_names, "Classification"].values
        for c in ("array_row", "array_col", "pxl_row", "pxl_col"):
            a.obs[c] = pos.reindex(a.obs_names)[c].values
        out.append(Section(sid, "wu2021", sid, WU_SUBTYPE[md.subtype.iloc[0]], a, None, "pathologist"))
    return out


def _find_image(folder: Path, name, visium_id: str) -> Path | None:
    """The metadata's file name, or (one deposited name has a typo in the slide id) the file that carries
    the slide number and capture area, e.g. '113.A1' for V10F24-113_A1."""
    if not folder.exists():
        return None
    if isinstance(name, str) and (folder / name).exists():
        return folder / name
    num, area = visium_id.split("-")[1].split("_")          # V10F24-113_A1 -> 113, A1
    for f in folder.iterdir():
        if f"{num}.{area}" in f.name or f"{num}_{area}" in f.name:
            return f
    if isinstance(name, str):                                 # BCSA1: ..._V19T26-012_KT_V1-Spot000001.jpg
        tail = name.split("_KT_")[-1] if "_KT_" in name else None
        for f in folder.iterdir():
            if tail and f.name.endswith(tail):
                return f
    return None


def load_li() -> list[Section]:
    base = UNPACKED / "li" / "spaceranger_output"
    meta = pd.read_excel(base / "Visium_metadata.xlsx")
    img_dir = UNPACKED / "li"
    out = []
    for r in meta.itertuples(index=False):
        sid = r[2]                                            # the "Visium ID" column
        d = base / sid / "outs"
        a = sc.read_10x_h5(d / "filtered_feature_bc_matrix.h5")
        a.var_names_make_unique()
        pos = read_positions(d / "spatial")
        a = a[a.obs_names.isin(pos.index[pos.in_tissue == 1])].copy()
        for c in ("array_row", "array_col", "pxl_row", "pxl_col"):
            a.obs[c] = pos.reindex(a.obs_names)[c].values
        a.obs["compartment"] = marker_compartment(a).values
        img = _find_image(img_dir / "Images" / "raw_images", r.Images, sid)
        out.append(Section(sid, "li2025", r.Patientid, r.type, a, img, "markers"))
    return out


def load_janesick_visium() -> list[Section]:
    h5 = RAW / "janesick_vis_filtered_feature_bc_matrix.h5"
    if not h5.exists() or not (UNPACKED / "janesick").exists():
        return []
    a = sc.read_10x_h5(h5)
    a.var_names_make_unique()
    spatial = next((UNPACKED / "janesick").rglob("tissue_positions*.csv")).parent
    pos = read_positions(spatial)
    a = a[a.obs_names.isin(pos.index[pos.in_tissue == 1])].copy()
    for c in ("array_row", "array_col", "pxl_row", "pxl_col"):
        a.obs[c] = pos.reindex(a.obs_names)[c].values
    a.obs["compartment"] = marker_compartment(a).values
    img = RAW / "janesick_vis_tissue_image.tif"
    return [Section("janesick_vis", "janesick2023", "block1", "HER2+", a, img if img.exists() else None,
                    "markers")]


def load_all() -> list[Section]:
    return load_wu() + load_li() + load_janesick_visium()


def scalefactors(section: Section) -> dict:
    if section.dataset == "wu2021":
        p = UNPACKED / "wu" / "spatial" / f"{section.id}_spatial" / "scalefactors_json.json"
    elif section.dataset == "li2025":
        p = (UNPACKED / "li" / "spaceranger_output" / section.id / "outs" / "spatial"
             / "scalefactors_json.json")
    else:
        p = next((UNPACKED / "janesick").rglob("scalefactors_json.json"))
    return json.loads(p.read_text())


def normalise(a: ad.AnnData, log: bool = True) -> ad.AnnData:
    b = a.copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    if log:
        sc.pp.log1p(b)
    return b

