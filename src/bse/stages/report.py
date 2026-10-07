"""Stage `report`: the verdicts on P1-P8 gathered from the compartment and histology summaries, and the
figures for the page. The page reads results/report/summary.json and never types a verdict in."""
import json
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from bse import config as C  # noqa: E402
from bse import provenance as P  # noqa: E402

OUT = C.RESULTS / "report"
COMP, HIST = C.RESULTS / "compartment", C.RESULTS / "histology"
INK, MUTED = "#16191D", "#5B646E"
COL = {"tumour": "#C06014", "stroma": "#6E7880", "immune": "#2D5BD1", "mixed": "#B9C0C6", "other": "#B9C0C6"}
ENC = {"phikon": "#0E7C7B", "imagenet": "#6E7880", "composition": "#C06014"}
SUB_LABEL = {"basal": "basal-like", "her2": "HER2-enriched", "luminal": "luminal"}
plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": "#B9C0C6", "xtick.color": MUTED,
                     "ytick.color": MUTED, "figure.dpi": 150})


def fig_compartment(cd: pd.DataFrame, floor: pd.DataFrame):
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    for ax, s in zip(axes, C.SUBTYPES, strict=True):
        d = cd[cd.subtype == s].sort_values("lfc_tumour")
        if d.empty:
            ax.set_title(SUB_LABEL[s]); continue
        colors = [COL.get(c, "#B9C0C6") for c in d.call]
        ax.barh(range(len(d)), d.lfc_tumour, color=colors, height=0.8)
        f = floor[floor.subtype == s]
        if not f.empty:
            ax.axvline(float(f.floor_p95.iloc[0]), color=INK, lw=.9, ls="--")
            ax.axvline(float(f.floor_median.iloc[0]), color="#B9C0C6", lw=.9)
        ax.axvline(C.COMPARTMENT_LFC, color=COL["tumour"], lw=.8, ls=":")
        ax.set_yticks([])
        ax.set_title(f"{SUB_LABEL[s]} · {len(d)} candidates", fontsize=10, color=INK)
        ax.set_xlabel("tumour log2 ratio (median over sections)")
    fig.text(0.5, -0.03, "bars = candidates, coloured by call (orange tumour, grey stroma/mixed, blue immune); dashed = matched-random "
             "floor 95th percentile for a set of this size, solid grey = its median; dotted = the 0.5 call threshold",
             ha="center", color=MUTED, fontsize=8.3)
    fig.tight_layout(); fig.savefig(OUT / "fig_compartment.png", bbox_inches="tight"); plt.close(fig)


def fig_outcomes(cd: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    d = cd.dropna(subset=["lfc_tumour"])
    ax = axes[0]
    for rep, lab, col in ((True, "replicated in Krug", "#0E7C7B"), (False, "not replicated", "#B9C0C6")):
        v = d[d.replicated == rep].lfc_tumour
        ax.hist(v, bins=25, alpha=.75, color=col, label=f"{lab} (n={len(v)})")
    ax.set_xlabel("tumour log2 ratio"); ax.set_ylabel("candidates"); ax.legend(frameon=False, fontsize=8)
    ax.set_title("compartment and replication", fontsize=10, color=INK)
    ax = axes[1]
    ax.scatter(d.lfc_tumour, d.mean_effect_subtype, s=14, color="#2D5BD1", alpha=.7)
    ax.axhline(-0.5, color="#B9C0C6", lw=.8, ls="--")
    ax.set_xlabel("tumour log2 ratio"); ax.set_ylabel("mean gene effect, subtype-matched lines")
    ax.set_title("compartment and dependency", fontsize=10, color=INK)
    fig.tight_layout(); fig.savefig(OUT / "fig_outcomes.png", bbox_inches="tight"); plt.close(fig)


def fig_histology(res: dict):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
    ax = axes[0]
    for kind, df in res.items():
        c = df[df.is_candidate].median_within_section_r.dropna().sort_values().values
        ax.plot(c, np.linspace(0, 1, len(c)), color=ENC[kind], lw=1.6, label=kind)
    ax.axvline(0.3, color="#B9C0C6", lw=.8, ls="--")
    ax.set_xlabel("median within-section Pearson r per candidate"); ax.set_ylabel("cumulative share of candidates")
    ax.legend(frameon=False, fontsize=8); ax.set_title("how well histology predicts the candidates", fontsize=10, color=INK)
    ax = axes[1]
    ph = res["phikon"]
    ctrl = ph[ph.is_control].set_index("gene").median_within_section_r
    ax.barh(range(len(ctrl)), ctrl.values, color=["#0E7C7B"] * len(ctrl))
    ax.set_yticks(range(len(ctrl)), ctrl.index)
    ax.axvline(0.5, color="#B9C0C6", lw=.8, ls="--")
    ax.set_xlabel("median within-section Pearson r (Phikon)"); ax.set_title("control genes", fontsize=10, color=INK)
    fig.tight_layout(); fig.savefig(OUT / "fig_histology.png", bbox_inches="tight"); plt.close(fig)


def run():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    comp = json.loads((COMP / "summary.json").read_text())
    cd = pd.read_csv(COMP / "candidates.csv")
    floor = pd.read_csv(COMP / "floor.csv")
    V = dict(comp["verdicts"])
    fig_compartment(cd, floor)
    fig_outcomes(cd)
    hist = None
    if (HIST / "summary.json").exists():
        hist = json.loads((HIST / "summary.json").read_text())
        V.update(hist["verdicts"])
        res = {k: pd.read_csv(HIST / f"genes_{k}.csv") for k in ("phikon", "imagenet", "composition")}
        fig_histology(res)
    for k in sorted(V):
        print(f"  {k}: {'holds' if V[k]['holds'] else 'fails'}")
    summary = {"verdicts": V, "compartment": {k: v for k, v in comp.items() if k != "verdicts"},
               "histology": {k: v for k, v in (hist or {}).items() if k != "verdicts"},
               "calls": cd.groupby("subtype").call.value_counts().unstack(fill_value=0).to_dict(),
               "n_candidates": int(len(cd))}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    P.log_run("report", {}, [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*"))], time.time() - t0)
    print(f"done in {time.time() - t0:,.0f}s")


if __name__ == "__main__":
    sys.exit(run())
