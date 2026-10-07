"""Every input file: where it comes from, its licence, and the local name it is saved under.

Nothing is downloaded until `bse run data` is called, and that stage records the sha256 of every file
in results/MANIFEST.sha256. Sizes are what Zenodo reported on 2026-10-07, for progress only.
"""
from dataclasses import dataclass

Z = "https://zenodo.org/api/records"


@dataclass(frozen=True)
class Source:
    key: str
    url: str
    local: str
    release: str
    licence: str
    size_mb: float
    note: str = ""


SOURCES = [
    # --- Wu et al. 2021 (Nat Genet): 6 Visium sections, pathologist annotation per spot -----------------
    Source("wu_filtered", f"{Z}/4739739/files/filtered_count_matrices.tar.gz/content", "wu2021_filtered_count_matrices.tar.gz",
           "Zenodo 4739739", "CC BY 4.0", 160.0, "Space Ranger filtered matrices, 6 sections"),
    Source("wu_spatial", f"{Z}/4739739/files/spatial.tar.gz/content", "wu2021_spatial.tar.gz",
           "Zenodo 4739739", "CC BY 4.0", 50.0, "spot positions, scale factors, low-resolution images"),
    Source("wu_metadata", f"{Z}/4739739/files/metadata.tar.gz/content", "wu2021_metadata.tar.gz",
           "Zenodo 4739739", "CC BY 4.0", 1.0, "pathologist annotation per spot"),
    # --- Li et al. 2025 (npj Precis Oncol): 23 Visium sections, tissue annotation ---------------------
    Source("li_spaceranger", f"{Z}/15211538/files/spaceranger_output.zip/content", "li2025_spaceranger_output.zip",
           "Zenodo 15211538", "CC BY 4.0", 1550.0, "Space Ranger outputs for 23 sections"),
    Source("li_images", f"{Z}/15211538/files/Images.zip/content", "li2025_images.zip",
           "Zenodo 15211538", "CC BY 4.0", 10140.0, "full-resolution H&E images"),
    # --- Janesick et al. 2023 (Nat Commun): Xenium + serial Visium CytAssist --------------------------
    Source("jan_xe1_matrix", f"{Z}/10076046/files/xe1_cell_feature_matrix.h5/content", "janesick_xe1_cell_feature_matrix.h5",
           "Zenodo 10076046", "CC BY 4.0", 10.0),
    Source("jan_xe1_cells", f"{Z}/10076046/files/xe1_cells.csv.gz/content", "janesick_xe1_cells.csv.gz",
           "Zenodo 10076046", "CC BY 4.0", 10.0),
    Source("jan_xe2_matrix", f"{Z}/10076046/files/xe2_cell_feature_matrix.h5/content", "janesick_xe2_cell_feature_matrix.h5",
           "Zenodo 10076046", "CC BY 4.0", 10.0),
    Source("jan_xe2_cells", f"{Z}/10076046/files/xe2_cells.csv.gz/content", "janesick_xe2_cells.csv.gz",
           "Zenodo 10076046", "CC BY 4.0", 10.0),
    Source("jan_celltypes", f"{Z}/10076046/files/Cell_Barcode_Type_Matrices.xlsx/content", "janesick_cell_barcode_types.xlsx",
           "Zenodo 10076046", "CC BY 4.0", 10.0, "the authors' cell-type annotation per Xenium cell"),
    Source("jan_vis_h5", f"{Z}/10076046/files/vis_filtered_feature_bc_matrix.h5/content", "janesick_vis_filtered_feature_bc_matrix.h5",
           "Zenodo 10076046", "CC BY 4.0", 30.0, "serial Visium CytAssist section"),
    Source("jan_vis_spatial", f"{Z}/10076046/files/vis_spatial.tar.gz/content", "janesick_vis_spatial.tar.gz",
           "Zenodo 10076046", "CC BY 4.0", 30.0),
    Source("jan_vis_image", f"{Z}/10076046/files/vis_tissue_image.tif/content", "janesick_vis_tissue_image.tif",
           "Zenodo 10076046", "CC BY 4.0", 1270.0, "full-resolution H&E of the Visium section"),
]
BY_KEY = {s.key: s for s in SOURCES}

# study 2's tables, copied at a pinned commit
STUDY2 = "https://raw.githubusercontent.com/Natalija-Stepurko/brca-target-evidence"
STUDY2_FILES = ["results/dossier/candidates.csv", "results/replicate/overlap.csv", "results/truth/t1_gene_effect.csv",
                "results/nominate/lists_discovery.csv", "results/replicate/lists_krug.csv"]
