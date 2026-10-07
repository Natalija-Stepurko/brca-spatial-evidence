# Related work and data access

Survey run 2026-10-07 as two parallel searches (public breast spatial-transcriptomics data and
compartment methods; histology-to-expression prediction and open pathology encoders), then
cross-checked. The design that follows from it is [`../docs/DESIGN.md`](../docs/DESIGN.md).

Strength tags: **[F]** read in the full text or primary page · **[Ab]** abstract or record only ·
**[S]** secondary (snippet, blog) · **[U]** unverified. Entries dated 2026 rest entirely on the search;
re-check before citing anything in public.

---

## 1. Public breast spatial transcriptomics with H&E (no access agreement)

| Dataset | n | Subtypes | Technology | Genes | Labels shipped | Licence, route |
|---|---|---|---|---|---|---|
| Wu et al. 2021, Nat Genet | 6 tumours, 6 sections | TNBC, ER+, HER2+ | Visium, 55 µm spots | whole transcriptome | pathologist annotation per spot | CC BY 4.0; Zenodo 4739739, ~0.9 GB [F] |
| Li et al. 2025, npj Precis Oncol | 4 patients, 23 sections | TNBC, HER2+ | Visium, fresh-frozen | whole transcriptome | QuPath tissue annotation (tumour / immune / stroma per cell, pathologist-reviewed) per spot | CC BY 4.0; Zenodo 15211538, 11.7 GB [F] |
| Janesick et al. 2023, Nat Commun | 2 FFPE blocks; 2 Xenium replicates + serial Visium CytAssist + scFFPE | HER2+ (one also ER+) | Xenium, subcellular; Visium 55 µm | Xenium 313-gene panel; Visium whole transcriptome | authors' cell types; post-Xenium H&E | CC BY 4.0; Zenodo 10076046, 7.6 GB [F] |
| Andersson et al. 2021, Nat Commun | 8 patients, 36 sections | HER2+ | legacy ST, 100 µm | whole transcriptome | one section per patient annotated | CC BY 4.0; Zenodo 4751624 [F] |
| He et al. 2020, Nat Biomed Eng | 23 patients, 68 sections | LumA, LumB, HER2+, TNBC | legacy ST, 100 µm | whole transcriptome | tumour annotation files | CC BY 4.0; Mendeley 29ntw7sh4r [F] |
| Wang et al. 2024, Nat Commun | 92 patients, 281 arrays | TNBC | legacy ST, 100 µm | whole transcriptome | 17-class pathologist annotation | open; Zenodo 8135721, 58 GB [F]; licence text [U] |
| Bassiouni et al. 2023 / 2025 | 14 / 10 patients | TNBC | Visium | whole transcriptome | images with coordinates | GEO GSE210616 [U]; Zenodo 15252874 CC BY 4.0 [F] |
| 10x demo sets (Visium FFPE breast; Visium HD DCIS; Xenium Prime 5K breast) | 1 block each | mixed | Visium / Visium HD / Xenium 5K | whole transcriptome / 5,001-gene panel | some pathologist regions (Zenodo 15411357) | CC BY 4.0 per 10x policy [S] |
| HEST-1k (Jaume et al. 2024) | 1,229 profiles, 26 organs, aligned WSIs | mixed | mixed | mixed | — | CC BY-NC-SA 4.0; Hugging Face, gated form [F] |

Not open: the MOSAIC window release holds no breast data. [F]

## 2. Compartment assignment in spatial data

- **Shipped annotations.** Wu 2021 and Andersson 2021 (pathologist, per spot or region); Wang 2024
  (17 classes); Li 2025 (QuPath random-trees classifier on H&E nuclei, pathologist-reviewed, giving
  per-cell tumour / immune / stroma counts in each spot). [F]
- **Reference-based deconvolution.** stereoscope (Wu, Andersson), cell2location, RCTD, Tangram. Li 2025
  benchmarked seven methods against its tissue annotation: cell2location, RCTD and stereoscope most
  concordant (median Spearman > 0.65 for tumour cells). RCTD runs on a 4-core CPU in minutes per
  section; cell2location is ~17× slower on CPU. [F]
- **CNV-based tumour calling.** SpatialInferCNV (Erickson et al. 2022, Nature); CalicoST (2024). [Ab]
- **Xenium cell typing.** Janesick: marker-based clusters, pathologist-verified ROIs, kNN label transfer
  from scFFPE. Cheng et al. 2025 (BMC Bioinformatics) benchmarked SingleR, Azimuth, RCTD and others on
  the 10x breast Xenium data; SingleR best and fastest. [F]
- **Controls used in these papers:** enrichment of deconvolved types inside pathologist regions;
  concordance of deconvolution with H&E cell counts; Visium-vs-Xenium marker agreement. None uses a
  permutation null for a per-gene compartment score. [F]

## 3. Drug targets in spatial data

Nothing breast-specific. ADC-target co-expression (ERBB2, TACSTD2, NECTIN4) mapped with Visium HD in
bladder cancer (2026) [Ab]; 72 ADC targets in 909 breast metastases by bulk RNA-seq and IHC, no spatial
compartment split (Nat Commun 2026) [Ab]. **No study scores a candidate's tumour-cell compartment and
tests whether it predicts replication between bulk cohorts or CRISPR dependency.** [F, by absence]

## 4. Predicting spatial expression from H&E

- **ST-Net** (He et al. 2020): 23 patients, fine-tuned DenseNet-121, 250 genes, leave-one-patient-out;
  top genes median Pearson 0.29–0.34 (GNAS, ACTG1, FASN, DDX5, XBP1); 102 of 250 genes positively
  correlated in ≥ 20 of 23 patients. [F]
- **BLEEP** (Xie et al. 2023, NeurIPS): frozen ResNet50 + contrastive retrieval, liver Visium, 4 sections;
  mean r on the held-out section ~0.17–0.22; best zonation genes 0.5–0.74. [F]
- **HEST-1k** (Jaume et al. 2024, NeurIPS D&B): 224-px tiles at 0.5 µm/px per spot, log1p top-50 HVGs,
  frozen encoder → PCA(256) → ridge, patient-stratified folds, mean Pearson over 50 genes. Breast IDC
  (4 Xenium samples): ResNet50-ImageNet 0.474, CTransPath 0.511, Phikon 0.533, CONCH 0.536, UNI 0.570,
  Virchow2 0.592, H-Optimus-0 0.598. The spread between encoders is smaller than between tasks; the
  ImageNet baseline is close behind. [F]
- **Nguyen et al. 2026 (arXiv:2609.32857):** aggregate Pearson is dominated by between-slide mean
  differences; a mean-only predictor scores 0.69 against 0.72 for ridge on frozen features; they
  recommend a slide-centred "local R²". [F]
- **Jang et al. 2026 (BMC Bioinformatics):** genes with coherent tissue-level patterns are inferable;
  genes without spatial organisation are not. [F]
- **Path2Space** (Shulman et al. 2026, Cell): CTransPath features + MLP ensemble predicting ~14,000 genes
  in breast cancer; links predicted spatial expression to treatment response. [Ab]
- Summary: predictable = compartment and lineage markers, zonation, proliferation structure, ECM;
  poorly predicted = genes without spatial structure, lowly expressed genes, within-compartment
  variation. Typical held-out per-gene r 0.1–0.3; best genes 0.5–0.7. [F]

## 5. Histology and nominated targets

No study asks whether H&E predicts nominated drug targets specifically, with target validation as the
endpoint. Closest: Path2Space (treatment response), ST-Net (predictable genes enriched in
pharmacogenomic pathways). [F, by absence]

## 6. Open pathology encoders on CPU

| Encoder | Size | Licence | Gated | CPU, 224-px tiles/s (4 threads) |
|---|---|---|---|---|
| timm ViT-B/16 ImageNet (baseline) | 86M | Apache-2.0 | no | ~13 measured on this machine under load |
| Phikon (owkin/phikon) | ViT-B/16 iBOT, 86M | non-commercial | no | ~13 |
| Phikon-v2 | ViT-L | CC BY-NC-ND 4.0 | [U] | ~3 |
| H0-mini | ViT-B/14 | CC BY-NC-ND 4.0 | yes | ~10 |
| CONCH, UNI, UNI2-h | ViT-B / L / H | CC BY-NC-ND 4.0 | yes | 13 / 3 / 1 |
| Virchow, Virchow2 | ViT-H | Apache-2.0 / CC BY-NC-ND | yes | ~1 |
| H-optimus-0, Prov-GigaPath | ViT-g | Apache-2.0 | yes | ~0.5–1 [U] |
| CTransPath | Swin-T | GPL-3, academic | n/a | [U] |

A public repository never redistributes weights; it depends on a model through the Hugging Face hub
with the user's own token. Ungated and benchmarked on breast: Phikon (non-commercial, stated in the
README). Fully clean: the ImageNet baseline only. [F]

## 7. Honest expectations and the controls that make a null informative

With 3–30 Visium sections, a frozen ViT-B and ridge: mean r over the 50 most variable genes 0.2–0.45 on a
held-out section; best compartment genes 0.5–0.7; many candidate genes near zero. Controls: hold out
whole sections (random spots are optimistic under spatial autocorrelation); report section-centred
correlations beside raw ones; within-section feature permutation per gene; an ImageNet encoder and a
compartment-composition baseline; matched random genes; positive (KRT/EPCAM, PTPRC, COL1A1, MKI67) and
negative (housekeeping) panels; a pre-registered gene list and thresholds. [F]

## 8. CPU feasibility

Visium: 2–4k spots under tissue per section. ~30 sections × ~3k tiles ≈ 90k tiles at ~13 tiles/s per
ViT-B ≈ 2 h per encoder; tile reading from the image is comparable; ridge on 90k × 256 is seconds.
Compartment scores and permutation nulls: minutes. [F]

## References

1. Wu SZ et al. 2021. Nat Genet 53:1334. doi:10.1038/s41588-021-00911-1; data doi:10.5281/zenodo.4739739
2. Li et al. 2025. npj Precis Oncol 9:310. doi:10.1038/s41698-025-01104-3; data doi:10.5281/zenodo.15211538
3. Janesick A et al. 2023. Nat Commun 14:8353. doi:10.1038/s41467-023-43458-x; data doi:10.5281/zenodo.10076046
4. Andersson A et al. 2021. Nat Commun 12:6012. doi:10.1038/s41467-021-26271-2
5. He B et al. 2020. Nat Biomed Eng 4:827. doi:10.1038/s41551-020-0578-x
6. Wang X et al. 2024. Nat Commun 15:10232. doi:10.1038/s41467-024-54145-w
7. Bassiouni R et al. 2023. Cancer Res 83:34. doi:10.1158/0008-5472.CAN-22-2682; 2025 Nat Commun doi:10.1038/s41467-025-61034-3
8. Jaume G et al. 2024. HEST-1k. NeurIPS Datasets and Benchmarks. arXiv:2406.16192
9. Xie R et al. 2023. BLEEP. NeurIPS. arXiv:2306.01859
10. Nguyen et al. 2026. arXiv:2609.32857
11. Jang et al. 2026. BMC Bioinformatics 27:168. doi:10.1186/s12859-026-06447-7
12. Shulman E et al. 2026. Cell. doi:10.1016/j.cell.2026.04.023
13. Filiot A et al. 2023. Phikon. medRxiv 10.1101/2023.07.21.23292757
14. Cable DM et al. 2022. RCTD. Nat Biotechnol. doi:10.1038/s41587-021-00830-w
15. Kleshchevnikov V et al. 2022. cell2location. Nat Biotechnol. doi:10.1038/s41587-021-01139-4
16. Erickson A et al. 2022. Nature 608:360. doi:10.1038/s41586-022-05023-2
17. Cheng J et al. 2025. BMC Bioinformatics. doi:10.1186/s12859-025-06044-0
