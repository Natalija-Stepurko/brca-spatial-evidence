"""The section inventory: one row per section with dataset, patient, subtype, spots, compartment
counts, label source and whether a full-resolution H&E image is available for the histology arm."""
import pandas as pd

from bse import config as C
from bse import sections as S

OUT = C.RESULTS / "data"


def write_inventory() -> pd.DataFrame:
    rows = []
    for s in S.load_all():
        comp = s.adata.obs.compartment.value_counts()
        rows.append({"section": s.id, "dataset": s.dataset, "patient": s.patient, "subtype": s.subtype,
                     "spots": int(s.adata.n_obs), "genes": int(s.adata.n_vars),
                     "label_source": s.label_source,
                     **{f"n_{k}": int(comp.get(k, 0)) for k in ("tumour", "stroma", "immune", "other")},
                     "full_res_image": str(s.image) if s.image else "", "histology_arm": bool(s.image)})
    inv = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    inv.to_csv(OUT / "sections.csv", index=False)
    print(inv[["section", "dataset", "subtype", "spots", "n_tumour", "n_stroma", "n_immune", "n_other",
               "label_source", "histology_arm"]].to_string(index=False))
    return inv
