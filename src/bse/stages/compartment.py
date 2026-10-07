"""Stage `compartment`: where each candidate is expressed (DESIGN §4) and how that relates to study 2's
outcomes (§5). Predictions P1-P5.

Per section: log-normalised counts; for every gene the tumour, stroma and immune log-fold-changes
(compartment spots vs all other spots) and the tumour fraction. A candidate's score is the median
tumour LFC over the sections of its subtype. Nulls: spot labels permuted within section; random gene
sets matched on mean expression and detection rate.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from bse import config as C
from bse import provenance as P
from bse import sections as S

OUT = C.RESULTS / "compartment"
STUDY2 = C.DATA / "study2"
SMOKE = os.environ.get("BSE_SMOKE") == "1"
N_PERM = 20 if SMOKE else C.N_PERMUTATIONS
N_RAND = 20 if SMOKE else C.N_MATCHED_RANDOM
N_BOOT = 20 if SMOKE else C.N_BOOTSTRAP
COMPARTMENTS = ("tumour", "stroma", "immune")
CONTROL_GENES = sorted({g for gs in C.CONTROLS.values() for g in gs})


# ----------------------------------------------------------------------------- per-section statistics
PSEUDO = 0.1          # counts per 10,000, added to both means before the log2 ratio


def log2_ratio(X, m, rest):
    a = np.asarray(X[m].mean(axis=0)).ravel()
    b = np.asarray(X[rest].mean(axis=0)).ravel()
    return np.log2((a + PSEUDO) / (b + PSEUDO))


def lfc_all(X, labels: np.ndarray) -> pd.DataFrame:
    """Per gene and compartment: log2 ratio of the mean normalised expression (counts per 10,000) in the
    compartment's spots to the mean in the other named compartments' spots (DESIGN §11); plus mean,
    detection and tumour fraction. X is normalised (not log-transformed), spots x genes (sparse)."""
    out = {}
    tot = np.asarray(X.sum(axis=0)).ravel()
    named = np.isin(labels, COMPARTMENTS)              # 'other' (necrosis, artefact, uncertain) is no comparator
    for comp in COMPARTMENTS:
        m = labels == comp
        rest = named & ~m
        if m.sum() < 10 or rest.sum() < 10:
            out[f"lfc_{comp}"] = np.full(X.shape[1], np.nan); continue
        out[f"lfc_{comp}"] = log2_ratio(X, m, rest)
    t = labels == "tumour"
    share_tumour_spots = t.mean()
    out["tumour_fraction"] = (np.asarray(X[t].sum(axis=0)).ravel() / np.maximum(tot, 1e-9)) / max(share_tumour_spots, 1e-9)
    out["mean_log"] = np.log1p(np.asarray(X.mean(axis=0)).ravel())
    out["detection"] = np.asarray((X > 0).mean(axis=0)).ravel()
    return pd.DataFrame(out)


def permuted_tlfc(X_sub, labels: np.ndarray, rng, n: int) -> np.ndarray:
    """Tumour LFC of the given genes (dense spots x genes) under n within-section label permutations."""
    Xd = X_sub.toarray() if hasattr(X_sub, "toarray") else np.asarray(X_sub)
    named = np.isin(labels, COMPARTMENTS)
    Xd = Xd[named]                                     # permute labels among the named compartments only
    k = int((labels[named] == "tumour").sum())
    out = np.empty((n, Xd.shape[1]))
    for i in range(n):
        idx = rng.permutation(Xd.shape[0])
        a, b = Xd[idx[:k]], Xd[idx[k:]]
        out[i] = np.log2((a.mean(axis=0) + PSEUDO) / (b.mean(axis=0) + PSEUDO))
    return out


# ----------------------------------------------------------------------------- study-2 inputs
def study2_tables():
    cand = pd.read_csv(STUDY2 / "results__dossier__candidates.csv")
    lists = pd.read_csv(STUDY2 / "results__nominate__lists_discovery.csv")
    lists = lists[lists.k == 50]
    t1 = pd.read_csv(STUDY2 / "results__truth__t1_gene_effect.csv", index_col=0)
    return cand, lists, t1


def call(row) -> str:
    hits = [c for c in COMPARTMENTS if row.get(f"lfc_{c}", np.nan) >= C.COMPARTMENT_LFC]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        return max(hits, key=lambda c: row[f"lfc_{c}"])
    return "mixed"


def profile_ci(y, x, cov, level=3.84):
    """Profile-likelihood 95% interval for the coefficient of x in a logistic model with covariates."""
    import statsmodels.api as sm
    Xc = sm.add_constant(cov)
    full = sm.Logit(y, np.column_stack([Xc, x])).fit(disp=0)
    b, ll = full.params[-1], full.llf
    def ll_at(beta):
        return sm.Logit(y, Xc, offset=beta * x).fit(disp=0).llf
    grid = np.linspace(b - 4, b + 4, 161)
    ok = [g for g in grid if 2 * (ll - ll_at(g)) <= level]
    return float(b), float(min(ok)), float(max(ok))


def run():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(C.SEED_LADDER)
    cand, lists, t1 = study2_tables()
    secs = S.load_all() if not SMOKE else S.load_wu()[:2]
    sub_of = {v: k for k, v in C.SUBTYPE_TO_SECTION.items()}         # TNBC -> basal, ...
    genes_of_interest = sorted(set(cand.gene) | set(lists.gene) | set(CONTROL_GENES))

    # 1. per-section statistics for every gene; permutation distributions for the genes of interest
    per_section, perm = {}, {}
    for s in secs:
        a = S.normalise(s.adata, log=False)
        labels = a.obs.compartment.values.astype(str)
        df = lfc_all(a.X.tocsr(), labels)
        df.index = a.var_names
        per_section[s.id] = df
        df.to_csv(OUT / f"genes_{s.id}.csv")
        goi = [g for g in genes_of_interest if g in a.var_names]
        perm[s.id] = pd.DataFrame(permuted_tlfc(a[:, goi].X, labels, rng, N_PERM), columns=goi)
        print(f"  {s.id:<16} {s.subtype:<5} spots={a.n_obs:>5} tumour={int((labels == 'tumour').sum()):>5}", flush=True)
    inv = pd.DataFrame([{"section": s.id, "dataset": s.dataset, "subtype": s.subtype, "study2_subtype": sub_of[s.subtype]}
                        for s in secs]).set_index("section")

    # 2. controls (P1): each control gene in its compartment in every section where that compartment exists
    ctrl_rows = []
    for sid, df in per_section.items():
        for comp, genes in C.CONTROLS.items():
            for g in genes:
                if g not in df.index:
                    continue
                r = df.loc[g]
                if comp == "housekeeping":
                    ok = all(not (r[f"lfc_{c}"] >= C.COMPARTMENT_LFC) for c in COMPARTMENTS if pd.notna(r[f"lfc_{c}"]))
                else:
                    ok = bool(r[f"lfc_{comp}"] >= C.COMPARTMENT_LFC) if pd.notna(r[f"lfc_{comp}"]) else None
                ctrl_rows.append({"section": sid, "gene": g, "expected": comp, "lfc_tumour": r.lfc_tumour,
                                  "lfc_stroma": r.lfc_stroma, "lfc_immune": r.lfc_immune, "passes": ok})
    ctrl = pd.DataFrame(ctrl_rows)
    ctrl.to_csv(OUT / "controls.csv", index=False)
    p1 = bool(ctrl.passes.dropna().all())
    print(f"  P1 controls: {int(ctrl.passes.dropna().sum())} of {int(ctrl.passes.notna().sum())} checks pass")

    # 3. candidate scores: median tumour LFC over the sections of the candidate's subtype
    def score(gene, subtype, col="lfc_tumour"):
        ids = inv.index[inv.study2_subtype == subtype]
        vals = [per_section[i].loc[gene, col] for i in ids if gene in per_section[i].index]
        return float(np.nanmedian(vals)) if vals and not all(np.isnan(vals)) else np.nan, len(vals)

    rows = []
    for r in cand.itertuples(index=False):
        ids = inv.index[inv.study2_subtype == r.subtype]
        tl, n = score(r.gene, r.subtype)
        st, _ = score(r.gene, r.subtype, "lfc_stroma")
        im, _ = score(r.gene, r.subtype, "lfc_immune")
        tf, _ = score(r.gene, r.subtype, "tumour_fraction")
        # permutation p: median over the subtype's sections of each permutation's TLFC
        pm = [perm[i][r.gene].values for i in ids if r.gene in perm[i].columns]
        p = np.nan
        if pm:
            med = np.nanmedian(np.vstack(pm), axis=0)
            p = float((np.sum(np.abs(med) >= abs(tl)) + 1) / (len(med) + 1)) if not np.isnan(tl) else np.nan
        row = {"subtype": r.subtype, "arm": r.arm, "rank": r.rank, "gene": r.gene, "n_sections": n,
               "lfc_tumour": tl, "lfc_stroma": st, "lfc_immune": im, "tumour_fraction": tf, "perm_p": p,
               "replicated": bool(r.in_krug_list_same_arm) if hasattr(r, "in_krug_list_same_arm") else None,
               "mean_effect_subtype": r.mean_effect_subtype, "hit_subtype": r.hit_subtype}
        row["call"] = call(row)
        rows.append(row)
    cd = pd.DataFrame(rows)
    cd.to_csv(OUT / "candidates.csv", index=False)
    print("  compartment calls per subtype:\n" + cd.groupby("subtype").call.value_counts().unstack(fill_value=0).to_string())

    # 4. matched-random floor per subtype (P2) and arms (P3)
    floor_rows, arm_rows = [], []
    for subtype in C.SUBTYPES:
        ids = inv.index[inv.study2_subtype == subtype]
        if len(ids) == 0:
            continue
        pooled = pd.concat([per_section[i][["lfc_tumour", "mean_log", "detection"]].assign(section=i) for i in ids])
        g = pooled.groupby(level=0).median(numeric_only=True)
        g = g.dropna(subset=["lfc_tumour"])
        g["q_expr"] = pd.qcut(g.mean_log.rank(method="first"), C.MATCH_STRATA["expression"], labels=False)
        g["q_det"] = pd.qcut(g.detection.rank(method="first"), C.MATCH_STRATA["detection"], labels=False)
        g["stratum"] = g.q_expr.astype(str) + "-" + g.q_det.astype(str)
        pool = {st: idx.values for st, idx in g.groupby("stratum").groups.items()}
        for arm, d in lists[lists.subtype == subtype].groupby("arm"):
            genes = [x for x in d.gene if x in g.index]
            obs = float(g.loc[genes, "lfc_tumour"].mean())
            dist = np.array([g.loc[[rng.choice(pool[g.loc[x, "stratum"]]) for x in genes], "lfc_tumour"].mean()
                             for _ in range(N_RAND)])
            arm_rows.append({"subtype": subtype, "arm": arm, "n_genes": len(genes), "mean_lfc_tumour": obs,
                             "median_lfc_tumour": float(g.loc[genes, "lfc_tumour"].median()),
                             "floor_p95": float(np.percentile(dist, 95)), "floor_median": float(np.median(dist)),
                             "floor_pct": float(np.mean(dist >= obs)), "beats_floor": bool(obs > np.percentile(dist, 95))})
        best_arm = cand[cand.subtype == subtype].arm.iloc[0]
        cg = [x for x in cand[cand.subtype == subtype].gene if x in g.index]
        obs = float(g.loc[cg, "lfc_tumour"].mean())
        dist = np.array([g.loc[[rng.choice(pool[g.loc[x, "stratum"]]) for x in cg], "lfc_tumour"].mean()
                         for _ in range(N_RAND)])
        frac_tumour = float((cd[cd.subtype == subtype].call == "tumour").mean())
        floor_rows.append({"subtype": subtype, "arm": best_arm, "n_sections": len(ids), "mean_lfc_tumour": obs,
                           "floor_p95": float(np.percentile(dist, 95)), "floor_median": float(np.median(dist)),
                           "floor_pct": float(np.mean(dist >= obs)), "beats_floor": bool(obs > np.percentile(dist, 95)),
                           "fraction_called_tumour": frac_tumour})
    floor = pd.DataFrame(floor_rows); floor.to_csv(OUT / "floor.csv", index=False)
    arms = pd.DataFrame(arm_rows); arms.to_csv(OUT / "arms.csv", index=False)
    print("  candidates vs matched floor:\n" + floor.round(3).to_string(index=False))
    print("  arms (median tumour LFC):\n" + arms.pivot(index="arm", columns="subtype", values="median_lfc_tumour").round(3).to_string())

    # 5. replication (P4) and dependency (P5)
    d = cd.dropna(subset=["lfc_tumour"]).copy()
    rep = {}
    if d.replicated.notna().any():
        y = d.replicated.astype(int).values
        x = ((d.lfc_tumour - d.lfc_tumour.mean()) / d.lfc_tumour.std()).values
        cov = pd.get_dummies(d.subtype, drop_first=True).astype(float).values
        try:
            b, lo, hi = profile_ci(y, x, cov)
            rep = {"odds_ratio_per_sd": float(np.exp(b)), "or_lo": float(np.exp(lo)), "or_hi": float(np.exp(hi))}
        except Exception as e:
            rep = {"error": str(e)[:200]}
        rep["median_lfc_replicated"] = float(d[d.replicated].lfc_tumour.median())
        rep["median_lfc_not_replicated"] = float(d[~d.replicated].lfc_tumour.median())
        rep["n_replicated"], rep["n_not"] = int(d.replicated.sum()), int((~d.replicated).sum())
        rep["median_difference"] = rep["median_lfc_replicated"] - rep["median_lfc_not_replicated"]
    dep = {}
    dd = d.dropna(subset=["mean_effect_subtype"])
    rho = float(spearmanr(dd.lfc_tumour, dd.mean_effect_subtype).statistic)
    boots = [spearmanr(dd.lfc_tumour.values[i], dd.mean_effect_subtype.values[i]).statistic
             for i in (rng.integers(0, len(dd), len(dd)) for _ in range(N_BOOT))]
    dep = {"spearman": rho, "lo": float(np.nanpercentile(boots, 2.5)), "hi": float(np.nanpercentile(boots, 97.5)), "n": int(len(dd))}

    # 6. verdicts
    seoi = C.SEOI_LFC
    arms_med = arms.pivot(index="arm", columns="subtype", values="median_lfc_tumour")
    p3 = all(("R+D" in arms_med.index and arms_med.loc["R+D", s] == arms_med[s].max()
              and arms_med.loc["R+D", s] - arms_med.loc["R", s] >= seoi) for s in arms_med.columns)
    V = {"P1": {"holds": p1, "rule": "every control gene in its compartment in every section where it exists"},
         "P2": {"holds": bool(floor.beats_floor.all() and (floor.fraction_called_tumour >= 0.5).all()),
                "rule": "candidates' mean tumour LFC above the matched floor's p95 in each subtype; >= half called tumour"},
         "P3": {"holds": bool(p3), "rule": f"R+D has the highest median tumour LFC in every subtype, >= {seoi} over R"},
         "P4": {"holds": bool(rep.get("median_difference", -9) >= seoi and rep.get("odds_ratio_per_sd", 0) >= C.SEOI_ODDS_RATIO),
                "rule": f"replicated minus non-replicated median LFC >= {seoi}; OR per SD >= {C.SEOI_ODDS_RATIO}", **rep},
         "P5": {"holds": bool(abs(rho) < 0.2), "rule": "|Spearman(tumour LFC, dependency)| < 0.2", **dep}}
    for k, v in V.items():
        print(f"  {k}: {'holds' if v['holds'] else 'fails'}")
    summary = {"verdicts": V, "sections": inv.reset_index().to_dict("records"), "floor": floor.to_dict("records"),
               "arms": arms.to_dict("records"), "calls": cd.groupby("subtype").call.value_counts().unstack(fill_value=0).to_dict(),
               "n_perm": N_PERM, "n_random": N_RAND}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    if SMOKE:
        print("smoke run; nothing logged"); return
    P.log_run("compartment", {"n_perm": N_PERM, "n_random": N_RAND, "seed": C.SEED_LADDER},
              [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*.csv")) + [OUT / "summary.json"]], time.time() - t0)
    print(f"done in {time.time() - t0:,.0f}s")


if __name__ == "__main__":
    sys.exit(run())
