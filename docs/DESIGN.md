# Design — where are the candidates expressed in tissue, and can histology see them?

*This document is the experimental design: what is measured, against what, and in what order. Related
work is in [`../research/literature.md`](../research/literature.md). Sections 1–9 are committed before
any dataset is downloaded for this study (the merge of this file into `main` is what "pre-specified"
refers to). Any later change to §1–§9 is listed in §11 with the date and the reason. §10 records what was
known before this commit.*

## 1. The questions

[Study 2](https://github.com/Natalija-Stepurko/brca-target-evidence) nominated breast-cancer drug-target
candidates per PAM50 subtype from tumour omics and found that subtype over-expression, in any layer, is
a weak guide to CRISPR dependency; candidates picked with DNA evidence were the most dependency-enriched
and the least reproducible between cohorts. Bulk tumour measurements mix malignant cells with stroma and
immune infiltrate, so a gene can look "up in the subtype" because the tumour's composition differs, not
its cancer cells. This study takes the candidates to intact tissue and asks:

> **Q1.** For each candidate, is its expression in tumour cells, stroma or immune infiltrate — and does
> nomination enrich for tumour-cell genes beyond what random genes of the same abundance show?
>
> **Q2.** Does a candidate's tumour-cell compartment predict whether its nomination replicated between
> cohorts, or how strongly subtype-matched cell lines depend on it?
>
> **Q3.** Can a pathology foundation model predict the candidates' local expression from the H&E image
> alone — and is what it can see simply the compartment structure?

The literature has no study of Q2, and none that asks Q3 of nominated drug targets specifically
(`literature.md` §3, §5).

## 2. Why the ladder again

Every compartment score and every histology correlation is read against a null, because both have
well-known inflations: spots are spatially autocorrelated, lineage genes dominate any "predictable"
list, and aggregate correlations are driven by between-section means (`literature.md` §4).

| Rung | What it is | What it establishes |
|---|---|---|
| **Matched-random floor** | 1,000 random gene sets matched to the candidates on mean expression and detection rate in the same sections | what any genes of that abundance score |
| **Label-permutation null** | compartment labels shuffled across spots within each section | the compartment score of a gene with no spatial structure |
| **Section hold-out** | histology models trained on all sections but one, scored on the held-out one, per section | what transfers to tissue the model never saw |
| **Within-section permutation** | image features shuffled across spots within a section | the histology correlation with no image information |
| **Plain baselines** | an ImageNet ViT-B encoder; a compartment-composition model | whether a pathology encoder adds anything to generic image features and to knowing the compartment |
| **Positive and negative controls** | EPCAM, KRT8, KRT18 (tumour); PTPRC (immune); COL1A1, DCN (stroma); housekeeping ACTB, GAPDH | the pipeline recovers what is known |

## 3. Data (all open, no access agreement)

| Role | Dataset | Sections | Subtypes | Compartment labels | Licence |
|---|---|---|---|---|---|
| **Discovery** | Wu et al. 2021, Visium, whole transcriptome | 6 tumours, 6 sections | TNBC, ER+, HER2+ | pathologist annotation per spot | CC BY 4.0, Zenodo 4739739 |
| **Replication** | Li et al. 2025, Visium, whole transcriptome | 4 patients, 23 sections | TNBC, HER2+ | computational tissue annotation: tumour / immune / stroma cell counts per spot, pathologist-reviewed | CC BY 4.0, Zenodo 15211538 |
| **Per-cell validation** | Janesick et al. 2023, Xenium (313-gene panel) with serial Visium CytAssist | 2 blocks | HER2+ (one ER+) | cell types from the authors' annotation | CC BY 4.0, Zenodo 10076046 |
| **Candidates** | study 2: `results/dossier/candidates.csv` (150 genes: 50 per subtype from the best arm), `results/replicate/overlap.csv`, `results/truth/t1_gene_effect.csv` | — | — | — | copied in with the study-2 commit pinned |

Subtype matching: a candidate list is scored in sections of its own subtype where the dataset has them
(basal-like ↔ TNBC; HER2-enriched ↔ HER2+; luminal ↔ ER+), and in all sections as a secondary analysis.
Xenium covers only the candidates on its panel; the number is recorded when stage `data` runs.

## 4. Compartment score (Q1)

Per section, spots are labelled tumour, stroma, immune or other from the shipped annotation (Li: the
majority cell class per spot; spots with no class are "other"). Counts are library-size normalised and
log1p-transformed. For gene *g* in section *k*:

- **tumour log-fold-change** TLFC(*g*,*k*) = mean log expression in tumour spots − mean in all
  non-tumour spots;
- **tumour fraction** TF(*g*,*k*) = share of *g*'s total counts that fall in tumour spots, divided by the
  share of all counts that fall in tumour spots (1 = no enrichment).

A gene's score is the median TLFC across the sections it is scored in. Its **compartment call** is tumour
if TLFC ≥ 0.5, stroma or immune if the analogous contrast for that compartment is ≥ 0.5, otherwise
"mixed". Two nulls: (i) spot labels permuted within section, 1,000 times, giving a p-value per gene;
(ii) 1,000 random gene sets matched on mean log expression (quintile) and detection rate (tercile),
giving the distribution of the candidates' mean TLFC under no selection.

## 5. Compartment and the study-2 outcomes (Q2)

Across the 150 candidates, with subtype as a covariate:

- **Replication.** Logistic regression of "in the replication cohort's list for the same arm" on the
  compartment score; reported as the odds ratio per SD of TLFC with a profile-likelihood interval, and as
  the difference in median TLFC between replicated and non-replicated candidates.
- **Dependency.** Spearman correlation between TLFC and the mean Chronos effect in subtype-matched lines
  (study 2's primary endpoint), with a 1,000-resample bootstrap interval over candidates.
- **By arm.** Median TLFC of candidates by the arm that nominated them (R, R+P, R+D, R-all), against the
  matched-random floor.

## 6. Histology (Q3)

- **Sections.** Every Visium section with a full-resolution H&E image at or near 0.5 µm/px; those with
  only the Space Ranger low-resolution image are excluded and listed. The count is recorded at stage `data`.
- **Tiles.** One 224 × 224 px tile at 0.5 µm/px centred on each spot (the HEST-1k protocol).
- **Encoders, frozen.** Phikon (owkin/phikon, ViT-B/16, ungated, non-commercial licence; weights are
  not redistributed) and an ImageNet ViT-B/16 from timm as the plain baseline. CLS features, cached once.
- **Model.** PCA to 256 dimensions, then ridge regression per gene (HEST-1k settings), trained on all
  sections but one and scored on the held-out section; every section is held out once.
- **Targets.** The 150 candidates, the control genes of §2, and the 50 most variable genes per section as
  context. Expression is log1p-normalised counts.
- **Metrics.** Pearson correlation per gene on the held-out section, raw and after centring per section
  (the "local" correlation); the median over candidates; the fraction of candidates whose correlation
  exceeds the 95th percentile of their within-section permutation null.
- **Composition baseline.** Ridge on the spot's compartment one-hot (tumour / stroma / immune / other)
  alone: how much of the histology signal is just the compartment.

## 7. Statistics

- Intervals are percentile intervals over resamples, never standard errors of a mean.
- Smallest effects of interest: 0.3 in TLFC; 0.05 in Pearson correlation; an odds ratio of 1.5 per SD.
- Multiplicity: everything is reported in full per subtype and per section; no claim rests on one cell.
- Seeds fixed (0 for the ladder, 1 for the histology splits).

## 8. Predictions, fixed before any data are downloaded

| | Prediction | If it fails |
|---|---|---|
| **P1** | The control genes land in their compartments in every section: EPCAM, KRT8, KRT18 tumour; PTPRC immune; COL1A1, DCN stroma; ACTB, GAPDH mixed | the compartment labels or the normalisation are wrong; debug before anything else |
| **P2** | The candidates' mean TLFC exceeds the matched-random floor's 95th percentile in each subtype; at least half the candidates are called "tumour" | nomination by subtype over-expression does not enrich for tumour-cell genes |
| **P3** | Candidates from the DNA-informed arm (R+D) have the highest median TLFC of the arms, by at least 0.3 over R | amplification-driven candidates are not more tumour-cell-specific |
| **P4** | Replicated candidates have a higher median TLFC than non-replicated ones, by at least 0.3; odds ratio per SD ≥ 1.5 | compartment does not explain which nominations transfer between cohorts |
| **P5** | TLFC and subtype-matched dependency correlate weakly: \|ρ\| < 0.2 | tumour-cell specificity tracks dependency, which would make spatial data a usable filter for targets |
| **P6** | Histology predicts the candidates poorly: median held-out Pearson < 0.3 raw and < 0.15 section-centred; the control compartment genes exceed 0.5 | the candidates are more visible in H&E than the field's benchmarks suggest |
| **P7** | Phikon exceeds the ImageNet ViT-B by at least 0.05 in median correlation over candidates, and by less than 0.05 over the composition baseline | a pathology encoder adds nothing over generic features, or adds more than the compartment structure explains |
| **P8** | A candidate's histology correlation is explained by its compartment score: Spearman(r, \|TLFC\|) ≥ 0.5 | histology sees something in the candidates that the compartment structure does not carry |

## 9. Stages

`data` (download with sha256 provenance, candidates copied from study 2 at a pinned commit, section
inventory incl. image resolution) → `compartment` (labels, normalisation, TLFC/TF, nulls; P1–P5) →
`tiles` (tile extraction, encoder features cached) → `histology` (ridge, hold-out, nulls; P6–P8) →
`report` (verdicts and figures) → the page. Compute: tiles for ~30 sections at ~13 tiles/s per ViT-B
is about 2 h per encoder on this machine; everything else is minutes.

## 10. What was known before this commit

A literature and data survey on 2026-10-07 established the datasets, licences and the HEST-1k protocol;
no spatial or image data were downloaded. Study 2's candidate tables already exist on this machine and
were produced without reference to spatial data.

## 11. Deviations from §1–§9

Each entry carries the date, what changed and why. No prediction or threshold has changed.

1. **2026-10-07 — the deposited data differ from the papers' counts (§3).** Wu et al.'s Zenodo deposit holds
   6 Visium sections: 4 TNBC and 2 ER+, no HER2+. Li et al.'s deposit holds 14 Visium sections from 3
   patients (4 TNBC, 10 HER2+), not 23 from 4. HER2-enriched candidates are therefore scored in Li's
   HER2+ sections (and Janesick's), luminal candidates in Wu's two ER+ sections only.
2. **2026-10-07 — Li et al. compartment labels (§3, §4).** The deposit carries no per-spot tissue-annotation
   table (the companion repository has one example section). Li's spots are labelled from a fixed marker
   panel scored on log-normalised counts (tumour: KRT19, KRT7, CDH1, CLDN4, ERBB3; immune: CD3E, CD2,
   CD79A, MS4A1, CD68, CD14; stroma: COL1A2, COL3A1, LUM, PDGFRB, SPARC), label = the highest standardised
   panel score when it leads the runner-up by 0.5 SD, otherwise "other". No control gene and no candidate
   is in the panel. Wu's pathologist labels remain the expression-independent source; Li's are flagged
   `markers` in `results/data/sections.csv`.
3. **2026-10-07 — the comparator and the scale of the compartment score (§4), triggered by P1.** On the
   first run the control genes pointed the right way but with magnitudes of 0.1–0.6, so the pre-registered
   threshold of 0.5 was mis-calibrated for the statistic as written (a difference of mean log1p values
   across mixed spots). Two changes: (i) each compartment is contrasted with the other two *named*
   compartments only — Wu's necrosis, artefact, uncertain and normal-gland spots are "other" and no longer
   a comparator, since necrotic spots still carry tumour RNA and normal glands express epithelial genes;
   (ii) the score is the log2 ratio of mean normalised expression (counts per 10,000, pseudocount 0.1)
   in the compartment to that in the comparator, the conventional pseudobulk log-fold-change. The 0.5
   threshold is unchanged and now applies to a log2 ratio. The permutation null is computed the same way.
4. **2026-10-07 — Wu's "Normal + stroma + lymphocytes" class** had been mapped to stroma by the word
   "stroma"; any class containing "normal" is "other".
5. **2026-10-07 — the histology arm's sections (§6).** Wu et al. ship their H&E only as a PDF and the
   Space Ranger low-resolution images, so by the rule in §6 their six sections are in the compartment
   analysis and not in the histology arm. Li et al. (full-resolution tif/jpg at 0.172 µm/px) and
   Janesick et al. (full-resolution tif) carry the histology arm.
6. **2026-10-07 — P1 as found.** 136 of 158 control checks pass. The failures are KRT8 and KRT18 in TNBC
   sections (luminal keratins; TNBC tumours express KRT5/14), the two sections whose "tumour" spots are
   almost all mixed "cancer + stroma + lymphocytes" (1142243F, CID44971), and a few housekeeping or COL1A1
   calls near the threshold. The labels and normalisation were checked and kept; the failures describe
   the data (spot-level mixing; the control choice) and the full table is `results/compartment/controls.csv`.

## 12. Limitations, stated in advance

- Six discovery sections across three subtypes, and a replication cohort with two subtypes only; the
  luminal arm rests on the ER+ sections of Wu et al. alone.
- Spot-level compartments are mixtures; a 55 µm spot can hold both tumour and stroma. The Xenium
  per-cell check covers only the candidates on a 313-gene panel.
- Histology correlations from a frozen encoder with a few thousand spots per section are modest by
  construction; the comparison that matters is against the baselines and the nulls, not the absolute value.
- Phikon's licence is non-commercial; the repository depends on it through the Hugging Face hub and
  never redistributes weights.
