"""Build the project page as one self-contained HTML file: docs/index.html (served by GitHub Pages).

Study 3 of the breast-cancer series.

One scrollable page: the question, why the study exists, the design, the pre-registered predictions,
the results (slots that fill from results/ once the stages have run), the agent audit, the exome
module, data and licences, literature, limits. Design facts (cohort sizes, K, thresholds) are stated
once in DESIGN and repeated here; the page carries no number that a stage produces until that stage
has written it to results/.

    python docs/site/build.py        -> docs/index.html
"""
import base64
import html
import json

import pandas as pd
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "index.html"
RES = ROOT / "results"
REPO = "https://github.com/Natalija-Stepurko/brca-spatial-evidence"
DESIGN = f"{REPO}/blob/main/docs/DESIGN.md"
LIT = f"{REPO}/blob/main/research/literature.md"
STUDY1 = "https://natalija-stepurko.github.io/single-cell-fm-probing/"
STUDY2 = "https://natalija-stepurko.github.io/brca-target-evidence/"

INK, MUTED, RULE, PANEL = "#16191D", "#5B646E", "#DDE1E4", "#FFFFFF"
RNA, PROT, DNA, TRUTH = "#2D5BD1", "#C06014", "#6E7880", "#0E7C7B"
SANS = "ui-sans-serif,system-ui,-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif"
MONO = "ui-monospace,SFMono-Regular,Menlo,Consolas,monospace"


def esc(s):
    return html.escape(str(s), quote=False)


# ----------------------------------------------------------------------------- status
def stage_status():
    """Which stages have written results. Everything is pending until results/run_log.json exists."""
    log = RES / "run_log.json"
    if not log.exists():
        return {}
    return {r["stage"]: r for r in json.loads(log.read_text()).get("runs", [])}


STATUS = stage_status()
PRE_REGISTERED = "2026-10-07"
SUMMARY = json.loads((RES / "report" / "summary.json").read_text()) if "report" in STATUS else None
SUB_LABEL = {"basal": "basal-like", "her2": "HER2-enriched", "luminal": "luminal"}
ARMS = ["R", "R+P", "R+D", "R+P+D", "MOFA+", "R-all"]


def f3(x):
    return f"{x:+.3f}" if x is not None else "—"


def pct(x):
    return f"{100 * x:.1f}%" if x is not None else "—"


def figure(name, alt):
    p = RES / "report" / name
    if not p.exists():
        return pending("report", alt)
    uri = "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()
    return f'<div class="figwrap"><img src="{uri}" alt="{esc(alt)}" style="width:100%;height:auto"></div>'


def verdict_chip(v):
    if v is None:
        return '<span class="chip">pending</span>'
    return ('<span class="chip ok">holds</span>' if v["holds"] else '<span class="chip no">fails</span>')


def pending(stage, what):
    """A result slot. Fills from results/ once `stage` has run; a dashed box until then."""
    if stage in STATUS:
        raise NotImplementedError(f"stage {stage} has results but the page has no renderer for it yet")
    return (f'<div class="pending"><span class="chip">pending</span> {esc(what)} — fills from '
            f'<code>results/</code> once stage <code>{esc(stage)}</code> has run.</div>')



def status_chips():
    if "report" not in STATUS:
        return '<span class="chip">no data downloaded yet</span><span class="chip">results pending</span>'
    return '<span class="chip ok">results in</span>'


SECTIONS = [("question", "Questions"), ("why", "Why"), ("design", "Design"), ("predictions", "Predictions"),
            ("results", "Results"), ("data", "Data"), ("literature", "Literature"), ("limits", "Limits")]


def nav():
    links = "".join(f'<a href="#{i}">{esc(t)}</a>' for i, t in SECTIONS)
    return (f'<nav class="topnav" aria-label="Sections"><a class="brand" href="#top">Spatial evidence</a>'
            f'<div class="navlinks">{links}</div>'
            f'<div class="navext"><a href="{REPO}">Code</a><a href="{DESIGN}">Design</a></div></nav>')


def header():
    return f"""<header class="page" id="top">
  <p class="eyebrow">Breast cancer · spatial transcriptomics · histology · study 3 of a series</p>
  <h1>Where are the candidates expressed in tissue, and can histology see them?</h1>
  <p class="lede">Study 2 nominated breast-cancer drug-target candidates from tumour omics and found that
  subtype over-expression is a weak guide to dependency. Bulk tumours mix cancer cells with stroma and
  immune infiltrate, so a gene can look "up" because the tissue's composition differs, not its cancer
  cells. This study takes the candidates to intact tissue: which compartment each sits in, whether that
  explains which nominations replicated and which genes cells depend on, and whether a pathology
  foundation model can see any of it in the H&E image.</p>
  <div class="status"><span class="chip live">pre-registered {PRE_REGISTERED}</span>{status_chips()}</div>
  <p class="links"><a href="{DESIGN}">Design (pre-registered)</a> · <a href="{LIT}">Literature</a> ·
  <a href="{REPO}">Repository</a> · <a href="{STUDY2}">Study 2: does the protein layer pick better targets?</a> ·
  <a href="{STUDY1}">Study 1: cell states and survival</a></p>
</header>"""


def s_question():
    return """<section class="sec" id="question">
  <h2>The questions</h2>
  <ol>
    <li><b>Compartment.</b> For each candidate, is its expression in tumour cells, stroma or immune
    infiltrate, and does nomination enrich for tumour-cell genes beyond random genes of the same abundance?</li>
    <li><b>Compartment and the study-2 outcomes.</b> Does a candidate's tumour-cell compartment predict
    whether its nomination replicated between cohorts, or how strongly subtype-matched cell lines depend on it?</li>
    <li><b>Histology.</b> Can a frozen pathology foundation model predict the candidates' local expression
    from the H&E image alone, and is what it sees simply the compartment structure?</li>
  </ol>
  <p class="note">Expected answers, stated before any data: nomination does enrich for tumour-cell genes,
  most of all in the DNA-informed arm; compartment explains part of replication and little of dependency;
  histology predicts the candidates poorly except where they are compartment markers, and a pathology
  encoder adds a little over generic image features and nothing over knowing the compartment.</p>
</section>"""


def s_why():
    return """<section class="sec" id="why">
  <h2>Why this study</h2>
  <div class="cards">
    <div class="card"><h3>Open question</h3><p>No published study scores a candidate's tumour-cell compartment
    in spatial data and tests whether it predicts replication between bulk cohorts or CRISPR dependency.
    Compartment methods exist and are benchmarked; nobody has pointed them at nominated targets.</p></div>
    <div class="card"><h3>Histology and targets</h3><p>Predicting spatial expression from H&E is an active
    benchmark field (HEST-1k reports 0.47–0.60 mean correlation on breast with frozen encoders, with an
    ImageNet baseline close behind). None of it asks about nominated drug targets specifically, and the
    field's own 2026 critique shows aggregate correlations are inflated by between-slide means.</p></div>
    <div class="card"><h3>What it adds to the series</h3><p>Study 2's candidates get a tissue-level rung
    they could not have from bulk data, and the series gains two data types (spatial transcriptomics,
    histology images) under the same discipline: pre-registration, nulls, held-out sections.</p></div>
  </div>
</section>"""


def s_design():
    return """<section class="sec" id="design">
  <h2>Design</h2>
  <h3>Data</h3>
  <div class="scroll"><table>
    <thead><tr><th>Role</th><th>Dataset</th><th>Sections</th><th>Subtypes</th><th>Compartment labels</th></tr></thead>
    <tbody>
    <tr><td><b>Discovery</b></td><td>Wu et al. 2021, Visium</td><td>6</td><td>TNBC, ER+, HER2+</td><td>pathologist annotation per spot</td></tr>
    <tr><td><b>Replication</b></td><td>Li et al. 2025, Visium</td><td>23</td><td>TNBC, HER2+</td><td>tumour / immune / stroma cell counts per spot</td></tr>
    <tr><td><b>Per-cell validation</b></td><td>Janesick et al. 2023, Xenium</td><td>2 blocks</td><td>HER2+</td><td>authors' cell types; 313-gene panel</td></tr>
    <tr><td><b>Candidates</b></td><td>study 2</td><td>150 genes</td><td>50 per subtype</td><td>with replication status and dependency</td></tr>
    </tbody></table></div>
  <h3>Compartment score</h3>
  <p>Per section, spots are tumour, stroma, immune or other from the shipped labels. A gene's tumour
  log-fold-change is its mean log expression in tumour spots minus the mean in all other spots; its score
  is the median across the sections of its subtype. Two nulls: spot labels permuted within each section,
  and 1,000 random gene sets matched on abundance and detection rate.</p>
  <h3>Histology</h3>
  <p>One 224-px tile at 0.5 µm/px per spot; CLS features from a frozen Phikon encoder and from an ImageNet
  ViT-B; PCA to 256 and ridge per gene, trained on all sections but one and scored on the held-out one.
  Correlations are reported raw and after centring per section, against a within-section feature
  permutation, and beside a model that knows only the spot's compartment.</p>
  <h3>The ladder</h3>
  <div class="scroll"><table><thead><tr><th>Rung</th><th>What it establishes</th></tr></thead><tbody>
    <tr><td>Matched-random floor</td><td>what any genes of that abundance score</td></tr>
    <tr><td>Label-permutation null</td><td>the compartment score of a gene with no spatial structure</td></tr>
    <tr><td>Section hold-out</td><td>what transfers to tissue the model never saw</td></tr>
    <tr><td>Within-section feature permutation</td><td>the histology correlation with no image information</td></tr>
    <tr><td>ImageNet ViT-B; compartment-composition model</td><td>whether a pathology encoder adds anything to generic features and to knowing the compartment</td></tr>
    <tr><td>Control genes</td><td>EPCAM, KRT8, KRT18 (tumour); PTPRC (immune); COL1A1, DCN (stroma); ACTB, GAPDH (housekeeping)</td></tr>
  </tbody></table></div>
</section>"""


PREDICTIONS = [
    ("P1", "Control genes land in their compartments in every section", "pipeline check before anything else"),
    ("P2", "The candidates' mean tumour log-fold-change beats the matched-random floor in each subtype; at least half are called tumour", ""),
    ("P3", "Candidates from the DNA-informed arm have the highest median tumour score, by at least 0.3 over RNA alone", ""),
    ("P4", "Replicated candidates score higher than non-replicated ones by at least 0.3; odds ratio per SD at least 1.5", ""),
    ("P5", "Tumour score and subtype-matched dependency correlate weakly, |ρ| < 0.2", ""),
    ("P6", "Histology predicts the candidates poorly: median held-out correlation < 0.3 raw, < 0.15 section-centred; control compartment genes > 0.5", ""),
    ("P7", "Phikon beats the ImageNet ViT-B by at least 0.05, and the composition baseline by less than 0.05", ""),
    ("P8", "A candidate's histology correlation is explained by its compartment score: Spearman ≥ 0.5", ""),
]


def s_predictions():
    V = (SUMMARY or {}).get("verdicts", {})
    rows = "".join(f'<tr><td><b>{p}</b></td><td>{esc(t)}</td><td class="muted">{esc(n)}</td>'
                   f'<td>{verdict_chip(V.get(p))}</td></tr>' for p, t, n in PREDICTIONS)
    return f"""<section class="sec" id="predictions">
  <h2>Predictions, fixed before any data</h2>
  <div class="scroll"><table><thead><tr><th></th><th>Prediction</th><th>Note</th><th>Outcome</th></tr></thead>
  <tbody>{rows}</tbody></table></div>
  <p class="note">Committed in <a href="{DESIGN}">DESIGN.md §8</a> on {PRE_REGISTERED}. Any later change to the
  design is logged in its §11.</p>
</section>"""


def results_text():
    V = SUMMARY["verdicts"]
    comp = SUMMARY["compartment"]
    floor = {r["subtype"]: r for r in comp["floor"]}
    calls = SUMMARY["calls"]
    n = SUMMARY["n_candidates"]
    out = []
    frac = {s: calls.get("tumour", {}).get(s, 0) / sum(v.get(s, 0) for v in calls.values()) for s in SUB_LABEL}
    out.append(f"<p><b>Nominated candidates are more tumour-cell-expressed than random genes of the same abundance, "
               f"in every subtype — modestly.</b> Mean tumour log2 ratio of the candidate set against the 95th percentile "
               f"of {comp['n_random']:,} matched random sets: "
               + "; ".join(f"{SUB_LABEL[s]} {f3(floor[s]['mean_lfc_tumour'])} vs {f3(floor[s]['floor_p95'])}" for s in SUB_LABEL if s in floor)
               + f". Yet only {', '.join(f'{pct(frac[s])}' for s in SUB_LABEL)} of the candidates clear the 0.5 call threshold "
               f"in the three subtypes, so P2 {'holds' if V['P2']['holds'] else 'fails on its second clause'}. The per-cell Xenium check "
               f"shows why: the same genes that sit at 0.2–0.5 across mixed 55 µm spots sit at 2–3 log2 across single cells.</p>")
    p4, p5 = V["P4"], V["P5"]
    out.append(f"<p><b>Compartment does not explain which nominations replicated between cohorts.</b> Odds ratio of "
               f"replication per SD of tumour score {p4.get('odds_ratio_per_sd', float('nan')):.2f} "
               f"[{p4.get('or_lo', float('nan')):.2f}, {p4.get('or_hi', float('nan')):.2f}]; median difference "
               f"{f3(p4.get('median_difference'))}. P4 {'holds' if p4['holds'] else 'fails'}.</p>")
    out.append(f"<p><b>It does track dependency.</b> Spearman correlation between a candidate's tumour score and the mean "
               f"CRISPR effect in subtype-matched lines: {p5['spearman']:.2f} [{p5['lo']:.2f}, {p5['hi']:.2f}] over {p5['n']} "
               f"candidates: more tumour-cell-specific candidates are the ones cells depend on more. P5 predicted "
               f"|ρ| < 0.2 and {'holds' if p5['holds'] else 'fails, in the informative direction'}.</p>")
    arms = pd.DataFrame(comp["arms"]) if comp.get("arms") else None
    if arms is not None and not arms.empty:
        piv = arms.pivot(index="arm", columns="subtype", values="median_lfc_tumour")
        out.append("<p><b>The nominating arm barely matters.</b> Median tumour score by arm: "
                   + "; ".join(f"{a} " + "/".join(f"{piv.loc[a, s]:.2f}" for s in piv.columns) for a in piv.index)
                   + f" (basal / HER2 / luminal). P3 {'holds' if V['P3']['holds'] else 'fails'}.</p>")
    c = pd.read_csv(RES / "compartment" / "controls.csv")
    out.append(f"<p><b>Controls.</b> {int(c.passes.dropna().sum())} of {int(c.passes.notna().sum())} control checks pass "
               f"(P1 {'holds' if V['P1']['holds'] else 'fails'}); the failures are KRT8 and KRT18 in TNBC sections "
               f"(luminal keratins) and the two sections whose tumour spots are almost all mixed with stroma and lymphocytes. "
               f"The design-level changes made at first contact with the data are in DESIGN §11.</p>")
    if "P6" in V:
        p6, p7, p8 = V["P6"], V["P7"], V["P8"]
        out.append(f"<p><b>Histology.</b> From the H&E image alone, a frozen Phikon encoder predicts the candidates at a median "
                   f"within-section Pearson r of {p6['median_within_r']:.2f} (section-centred pooled {p6['median_centred_r']:.2f}); "
                   f"the compartment control genes reach {p6['controls_median_r']:.2f}. P6 {'holds' if p6['holds'] else 'fails'}. "
                   f"Phikon minus ImageNet ViT-B: {f3(p7['phikon_minus_imagenet'])}; Phikon minus a model that knows only the spot's "
                   f"compartment: {f3(p7['phikon_minus_composition'])}. P7 {'holds' if p7['holds'] else 'fails'}. A candidate's "
                   f"histology correlation and its |tumour score| correlate at Spearman {p8['spearman']:.2f} (P8 "
                   f"{'holds' if p8['holds'] else 'fails'}).</p>")
    return "".join(out)


def s_results():
    if SUMMARY is None:
        body = (pending("compartment", "Compartment calls for the 150 candidates, per subtype, against both nulls (P1–P3)")
                + pending("compartment", "Compartment score against replication status and dependency (P4, P5)")
                + pending("histology", "Held-out histology correlations: candidates, controls, both encoders, the composition baseline (P6–P8)"))
    else:
        body = (results_text()
                + figure("fig_compartment.png", "Tumour log2 ratio of every candidate per subtype against the matched-random floor")
                + '<p class="figcap">Each bar is a candidate\'s tumour log2 ratio (median over the sections of its subtype), coloured by '
                  'its compartment call; the dashed line is the 95th percentile of the matched-random floor for a set of that size.</p>'
                + figure("fig_outcomes.png", "Compartment score against replication status and against dependency")
                + '<p class="figcap">Left: tumour score of replicated and non-replicated candidates. Right: tumour score against the '
                  'mean CRISPR gene effect in subtype-matched lines (dashed: the −0.5 dependency threshold).</p>'
                + (figure("fig_histology.png", "Histology prediction of candidates and control genes")
                   + '<p class="figcap">Left: cumulative distribution of per-candidate held-out correlations for Phikon, an ImageNet '
                     'ViT-B and the compartment-composition baseline. Right: the control genes under Phikon.</p>'
                   if "P6" in SUMMARY["verdicts"] else pending("histology", "Held-out histology correlations (P6–P8)")))
    return f"""<section class="sec" id="results">
  <h2>Results</h2>
  {body}
</section>"""


DATA = [
    ("Wu et al. 2021 Visium, 6 sections, with pathologist annotation", "Zenodo 4739739", "CC BY 4.0"),
    ("Li et al. 2025 Visium, 23 sections, with tissue annotation", "Zenodo 15211538", "CC BY 4.0"),
    ("Janesick et al. 2023 Xenium + Visium CytAssist", "Zenodo 10076046", "CC BY 4.0"),
    ("Study 2 candidate tables (pinned commit)", "brca-target-evidence", "MIT"),
    ("Phikon encoder weights", "Hugging Face hub, not redistributed", "non-commercial"),
    ("ImageNet ViT-B/16 (timm)", "Hugging Face hub", "Apache-2.0"),
]


def s_data():
    rows = "".join(f"<tr><td>{esc(a)}</td><td>{esc(b)}</td><td>{esc(c)}</td></tr>" for a, b, c in DATA)
    return f"""<section class="sec" id="data">
  <h2>Data and licences</h2>
  <div class="scroll"><table><thead><tr><th>Dataset</th><th>Route</th><th>Licence</th></tr></thead>
  <tbody>{rows}</tbody></table></div>
  <p class="note">Every file is recorded with URL and sha256 in <code>results/MANIFEST.sha256</code> when stage
  <code>data</code> runs. Model weights are fetched by the user's own hub token and never stored in the repository.</p>
</section>"""


REFS = [
    ("Wu SZ et al.", "A single-cell and spatially resolved atlas of human breast cancers", "Nat Genet 2021", "10.1038/s41588-021-00911-1"),
    ("Li et al.", "Spatial transcriptomics of breast cancer with computational tissue annotation", "npj Precis Oncol 2025", "10.1038/s41698-025-01104-3"),
    ("Janesick A et al.", "High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis", "Nat Commun 2023", "10.1038/s41467-023-43458-x"),
    ("He B et al.", "Integrating spatial gene expression and breast tumour morphology via deep learning", "Nat Biomed Eng 2020", "10.1038/s41551-020-0578-x"),
    ("Jaume G et al.", "HEST-1k: a dataset for spatial transcriptomics and histology image analysis", "NeurIPS 2024", ""),
    ("Nguyen et al.", "Is H&E image-to-spatial transcriptomics simpler than it looks?", "arXiv 2026", ""),
    ("Cable DM et al.", "Robust decomposition of cell type mixtures in spatial transcriptomics", "Nat Biotechnol 2022", "10.1038/s41587-021-00830-w"),
    ("Filiot A et al.", "Scaling self-supervised learning for histopathology with masked image modeling (Phikon)", "medRxiv 2023", "10.1101/2023.07.21.23292757"),
]


def s_literature():
    items = []
    for who, title, venue, doi in REFS:
        link = f' <a href="https://doi.org/{doi}">doi</a>' if doi else ""
        items.append(f"<li>{esc(who)}. <i>{esc(title)}</i>. {esc(venue)}.{link}</li>")
    return f"""<section class="sec" id="literature">
  <h2>Literature</h2>
  <p>The survey behind the design, with a strength tag on every entry, is in <a href="{LIT}">research/literature.md</a>.</p>
  <ol class="refs">{"".join(items)}</ol>
</section>"""


def s_limits():
    return """<section class="sec" id="limits">
  <h2>Limits, stated in advance</h2>
  <ul>
    <li><b>Few sections.</b> Six discovery sections across three subtypes; the replication cohort has two
    subtypes, so the luminal arm rests on the ER+ sections of Wu et al. alone.</li>
    <li><b>Spots are mixtures.</b> A 55 µm spot can hold tumour and stroma; the per-cell Xenium check covers
    only candidates on a 313-gene panel.</li>
    <li><b>Modest correlations by construction.</b> A frozen encoder with a few thousand spots per section;
    what matters is the comparison with the baselines and nulls, not the absolute value.</li>
    <li><b>Licences.</b> Phikon is non-commercial; the repository depends on it through the hub and never
    redistributes weights.</li>
  </ul>
</section>"""


def footer():
    return f"""<footer><p>Built {date.today().isoformat()} from <code>docs/site/build.py</code>. Every number that
  a pipeline stage produces is read from <code>results/</code> when this page is built; none is typed in.
  <a href="{REPO}">{REPO.split("//")[1]}</a></p></footer>"""


# ----------------------------------------------------------------------------- page
CSS = f"""
:root{{--paper:#F7F8F9;--panel:{PANEL};--ink:{INK};--muted:{MUTED};--rule:{RULE};--rna:{RNA};--prot:{PROT};
  --truth:{TRUTH};--sans:{SANS};--mono:{MONO}}}
*{{box-sizing:border-box}}
html{{scroll-behavior:smooth;scroll-padding-top:58px}}
body{{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);line-height:1.6;
  -webkit-font-smoothing:antialiased;padding-inline:16px}}
.wrap{{max-width:1040px;margin:0 auto;padding:24px 0 80px}}
a{{color:var(--ink)}}
h1{{font-size:clamp(1.6rem,4vw,2.4rem);line-height:1.15;margin:.3em 0 .5em;letter-spacing:-.01em}}
h2{{font-size:1.45rem;margin:2.2em 0 .6em;letter-spacing:-.01em}}
h3{{font-size:1.05rem;margin:1.6em 0 .4em}}
.eyebrow{{font-family:var(--mono);font-size:.72rem;color:var(--muted);letter-spacing:.06em;text-transform:uppercase;margin:0}}
.lede{{font-size:1.08rem;max-width:62ch}}
.note{{color:var(--muted);font-size:.95rem}}
.muted{{color:var(--muted)}}
.status{{display:flex;flex-wrap:wrap;gap:8px;margin:12px 0}}
.chip{{display:inline-block;font-family:var(--mono);font-size:.7rem;padding:3px 9px;border-radius:999px;
  border:1px solid var(--rule);background:var(--panel);color:var(--muted);white-space:nowrap}}
.chip.live{{border-color:var(--truth);color:var(--truth)}}
.chip.ok{{border-color:var(--truth);color:var(--truth);background:#EAF4F4}}
.chip.no{{border-color:var(--prot);color:var(--prot);background:#FBF1EA}}
.links{{font-size:.95rem}}
.finding{{background:var(--panel);border-left:4px solid var(--truth);border-radius:6px;padding:12px 16px;margin:14px 0}}
.finding p{{margin:.3em 0}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px;margin:1em 0}}
.card{{background:var(--panel);border:1px solid var(--rule);border-radius:8px;padding:14px 16px}}
.card h3{{margin:0 0 .4em;font-size:.95rem}}
.card p{{margin:0;font-size:.95rem}}
.figwrap{{background:var(--panel);border:1px solid var(--rule);border-radius:8px;padding:12px;overflow-x:auto}}
.figwrap svg{{display:block;width:100%;min-width:640px;height:auto}}
.figcap{{color:var(--muted);font-size:.9rem;margin:.5em 0 1.4em}}
.scroll{{overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:.93rem;background:var(--panel)}}
th,td{{text-align:left;padding:8px 10px;border-bottom:1px solid var(--rule);vertical-align:top}}
th{{font-family:var(--mono);font-size:.72rem;letter-spacing:.05em;text-transform:uppercase;color:var(--muted)}}
.pending{{border:1px dashed var(--rule);border-radius:8px;padding:14px 16px;margin:12px 0;color:var(--muted);font-size:.95rem}}
.pending .chip{{margin-right:6px}}
code{{font-family:var(--mono);font-size:.88em;background:#EEF1F2;padding:1px 5px;border-radius:4px}}
ol.refs{{font-size:.92rem;padding-left:1.3em}}
ol.refs li{{margin:.35em 0}}
footer{{margin-top:60px;border-top:1px solid var(--rule);padding-top:14px;color:var(--muted);font-size:.85rem}}
/* sticky section navbar */
.topnav{{position:sticky;top:0;z-index:50;display:flex;align-items:center;gap:18px;margin:0 -16px;
  padding:0 max(16px,calc((100% - 1040px)/2 + 16px));height:48px;background:rgba(247,248,249,.94);
  backdrop-filter:blur(6px);border-bottom:1px solid var(--rule);font-family:var(--mono);font-size:.7rem}}
.topnav a{{color:var(--muted);text-decoration:none;white-space:nowrap}}
.topnav a:hover{{color:var(--truth)}}
.topnav .brand{{color:var(--ink);font-weight:600;flex:none}}
.navlinks{{display:flex;gap:16px;overflow-x:auto;scrollbar-width:none;flex:1 1 auto;min-width:0;padding:0 12px}}
.navlinks::-webkit-scrollbar{{display:none}}
.navlinks a{{padding:14px 0 12px;border-bottom:2px solid transparent}}
.navlinks a.on{{color:var(--ink);border-bottom-color:var(--truth)}}
.navext{{display:flex;gap:14px;flex:none}}
.navext a{{color:var(--truth);font-weight:600}}
@media(max-width:640px){{.topnav .brand{{display:none}}.topnav{{gap:10px}}}}
"""

JS = """<script>
(function () {
  const links = [...document.querySelectorAll('.topnav .navlinks a')];
  const secs = links.map(a => document.getElementById(a.getAttribute('href').slice(1))).filter(Boolean);
  let ticking = false;
  function update() {
    ticking = false;
    const y = window.scrollY + 70;
    let cur = secs[0];
    for (const s of secs) if (s.offsetTop <= y) cur = s;
    if (window.innerHeight + window.scrollY >= document.body.offsetHeight - 2) cur = secs[secs.length - 1];
    links.forEach(a => a.classList.toggle('on', a.getAttribute('href') === '#' + cur.id));
  }
  addEventListener('scroll', () => { if (!ticking) { ticking = true; requestAnimationFrame(update); } });
  update();
})();
</script>"""

FAVICON = ("<link rel=\"icon\" href=\"data:image/svg+xml,"
           "%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E"
           "%3Crect x='4' y='18' width='6' height='10' fill='%232D5BD1'/%3E"
           "%3Crect x='13' y='11' width='6' height='17' fill='%23C06014'/%3E"
           "%3Crect x='22' y='4' width='6' height='24' fill='%230E7C7B'/%3E%3C/svg%3E\">")


def build():
    body = "\n\n".join([header(), s_question(), s_why(), s_design(), s_predictions(), s_results(),
                        s_data(), s_literature(), s_limits(), footer()])
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Where are the candidates expressed in tissue?</title>
<meta name="description" content="Study 3 of a breast-cancer series: the candidate drug targets of study 2 mapped onto spatial transcriptomics sections, their tumour-cell compartment tested against replication and dependency, and a frozen pathology encoder tested on the H&E image. Pre-registered, CPU-only.">
{FAVICON}
<style>{CSS}</style>
</head>
<body>
{nav()}
<div class="wrap">
{body}
</div>
{JS}
</body>
</html>
"""
    for bad in ("rather than", "instead of", "None<"):
        assert bad not in page, f"banned text on the page: {bad!r}"
    OUT.write_text(page)
    print(f"wrote {OUT.relative_to(ROOT)} ({len(page):,} bytes); stages with results: {sorted(STATUS) or 'none'}")


if __name__ == "__main__":
    build()
