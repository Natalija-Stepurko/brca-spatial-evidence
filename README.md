# Where are the candidates expressed in tissue, and can histology see them?

Study 3 of a breast-cancer series. Study 1, [single-cell cell states and survival](https://github.com/Natalija-Stepurko/single-cell-fm-probing),
asked whether foundation-model cell states predict patient survival. Study 2,
[does the protein layer pick better drug targets?](https://github.com/Natalija-Stepurko/brca-target-evidence),
nominated drug-target candidates per subtype from tumour omics and found that subtype over-expression is a
weak guide to dependency. This study takes those candidates to intact tissue.

**Status (2026-10-07):** design pre-registered in [`docs/DESIGN.md`](docs/DESIGN.md) before any data were
downloaded. No results yet. Project page: <https://natalija-stepurko.github.io/brca-spatial-evidence/>.

## The questions

1. **Compartment.** For each candidate, is its expression in tumour cells, stroma or immune infiltrate, and
   does nomination enrich for tumour-cell genes beyond random genes of the same abundance?
2. **Compartment and the study-2 outcomes.** Does a candidate's tumour-cell compartment predict whether its
   nomination replicated between cohorts, or how strongly subtype-matched cell lines depend on it?
3. **Histology.** Can a frozen pathology foundation model predict the candidates' local expression from the
   H&E image alone, and is what it sees simply the compartment structure?

The literature has no study of question 2, and none that asks question 3 of nominated drug targets.

## What will be done

- **Data.** Wu et al. 2021 Visium (6 sections, three subtypes, pathologist labels per spot; discovery);
  Li et al. 2025 Visium (23 sections, tumour / immune / stroma cell counts per spot; replication);
  Janesick et al. 2023 Xenium for per-cell validation of the candidates on its panel. All CC BY 4.0, no
  access agreement. Study 2's candidate tables are copied in at a pinned commit.
- **Compartment score.** Per section, a gene's tumour log-fold-change (tumour spots vs the rest) and
  tumour fraction, against a within-section label-permutation null and 1,000 abundance-matched random gene
  sets; control genes for each compartment.
- **Compartment and outcomes.** Logistic regression of replication status on the score; Spearman
  correlation with subtype-matched dependency; the score by nominating arm.
- **Histology.** One 224-px tile at 0.5 µm/px per spot; CLS features from Phikon and from an ImageNet
  ViT-B; PCA-256 and ridge per gene, every section held out once; correlations raw and section-centred,
  against a within-section feature permutation and a compartment-composition baseline.
- **Predictions P1–P8** with smallest effects of interest, fixed in `docs/DESIGN.md` §8.

## Repository

| Path | What |
|---|---|
| [`docs/DESIGN.md`](docs/DESIGN.md) | the pre-registered design; §11 logs any later change |
| [`research/literature.md`](research/literature.md) | the survey behind the design, with a strength tag on every entry |
| [`docs/site/build.py`](docs/site/build.py) → `docs/index.html` | the project page; result slots fill from `results/` once stages have run |
| `src/`, `tests/` | the pipeline (to come, one pull request per stage) |

## Data and licences

| Dataset | Route | Licence |
|---|---|---|
| Wu et al. 2021 Visium, 6 sections, pathologist annotation | Zenodo 4739739 | CC BY 4.0 |
| Li et al. 2025 Visium, 23 sections, tissue annotation | Zenodo 15211538 | CC BY 4.0 |
| Janesick et al. 2023 Xenium + Visium CytAssist | Zenodo 10076046 | CC BY 4.0 |
| Study 2 candidate tables | brca-target-evidence, pinned commit | MIT |
| Phikon encoder weights | Hugging Face hub, fetched with the user's token, never redistributed | non-commercial |
| ImageNet ViT-B/16 (timm) | Hugging Face hub | Apache-2.0 |

Every downloaded file is recorded with URL and sha256 in `results/MANIFEST.sha256`.

## Licence

MIT — see [`LICENSE`](LICENSE).
