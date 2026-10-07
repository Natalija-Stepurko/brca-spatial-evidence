"""Stage `histology`: can a frozen pathology encoder predict the candidates from H&E (DESIGN §6)? P6-P8.

For each encoder (Phikon, ImageNet ViT-B) and for the compartment-composition baseline: PCA to 256
dimensions fitted on the training sections, ridge per gene (alpha chosen by inner leave-one-section-out
on the training sections), every section held out once. Targets: the 150 candidates, the control
genes, and the 50 most variable genes of each section. Metrics per gene: Pearson within the held-out
section (median over sections), pooled Pearson over all held-out predictions raw and after centring
per section, and a within-section permutation null of the predictions.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge

from bse import config as C
from bse import provenance as P
from bse import sections as S

OUT = C.RESULTS / "histology"
FEAT = C.DATA / "features"
STUDY2 = C.DATA / "study2"
SMOKE = os.environ.get("BSE_SMOKE") == "1"
N_PERM = 50 if SMOKE else C.N_PERMUTATIONS
ALPHAS = (1.0, 10.0, 100.0, 1000.0)
CONTROL_GENES = sorted({g for gs in C.CONTROLS.values() for g in gs})


def targets(secs, cand_genes):
    """Union of candidates, controls and each section's top variable genes, restricted to genes present
    in every histology section."""
    common = set.intersection(*[set(s.adata.var_names) for s in secs])
    hv = set()
    for s in secs:
        a = S.normalise(s.adata)
        sc.pp.highly_variable_genes(a, n_top_genes=C.N_TOP_VARIABLE)
        hv |= set(a.var_names[a.var.highly_variable])
    genes = sorted((set(cand_genes) | set(CONTROL_GENES) | hv) & common)
    return genes, sorted(set(cand_genes) & common), sorted(set(CONTROL_GENES) & common)


def expression(secs, genes):
    """Log-normalised expression of the target genes per section, spots in the feature file's order."""
    out = {}
    for s in secs:
        z = np.load(FEAT / f"{s.id}.npz", allow_pickle=True)
        a = S.normalise(s.adata)
        a = a[list(z["barcodes"]), genes]
        out[s.id] = np.asarray(a.X.todense()) if hasattr(a.X, "todense") else np.asarray(a.X)
    return out


def features(secs, kind):
    out = {}
    for s in secs:
        z = np.load(FEAT / f"{s.id}.npz", allow_pickle=True)
        if kind == "composition":
            lab = s.adata.obs.compartment.reindex(z["barcodes"]).values.astype(str)
            out[s.id] = np.column_stack([(lab == c).astype(float) for c in ("tumour", "stroma", "immune", "other")])
        else:
            out[s.id] = z[kind]
    return out


def fit_predict(F, Y, train_ids, test_id, use_pca):
    Xtr = np.concatenate([F[i] for i in train_ids]); Ytr = np.concatenate([Y[i] for i in train_ids])
    Xte = F[test_id]
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-8
    Xtr, Xte = (Xtr - mu) / sd, (Xte - mu) / sd
    if use_pca and Xtr.shape[1] > C.PCA_DIM:
        pca = PCA(n_components=C.PCA_DIM, random_state=C.SEED_SPLITS).fit(Xtr)
        Xtr, Xte = pca.transform(Xtr), pca.transform(Xte)
    # inner leave-one-section-out on the training sections to choose alpha
    best, best_score = ALPHAS[0], -np.inf
    if len(train_ids) >= 3:
        sizes = [len(F[i]) for i in train_ids]
        bounds = np.cumsum([0] + sizes)
        for alpha in ALPHAS:
            scores = []
            for k in range(len(train_ids)):
                m = np.ones(len(Xtr), bool); m[bounds[k]:bounds[k + 1]] = False
                r = Ridge(alpha=alpha).fit(Xtr[m], Ytr[m])
                pred = r.predict(Xtr[~m])
                scores.append(np.nanmean([_pearson(Ytr[~m][:, j], pred[:, j]) for j in range(Ytr.shape[1])]))
            if np.nanmean(scores) > best_score:
                best, best_score = alpha, np.nanmean(scores)
    model = Ridge(alpha=best).fit(Xtr, Ytr)
    return model.predict(Xte), best


def _pearson(y, p):
    if np.std(y) == 0 or np.std(p) == 0:
        return np.nan
    return float(pearsonr(y, p).statistic)


def evaluate(kind, secs, genes, Y, rng, use_pca=True):
    F = features(secs, kind)
    ids = [s.id for s in secs]
    per_section = {}                       # section -> array of Pearson per gene
    preds, truth, sect = [], [], []
    alphas = {}
    for test in ids:
        train = [i for i in ids if i != test]
        pred, alpha = fit_predict(F, Y, train, test, use_pca)
        alphas[test] = alpha
        per_section[test] = np.array([_pearson(Y[test][:, j], pred[:, j]) for j in range(len(genes))])
        preds.append(pred); truth.append(Y[test]); sect.append(np.full(len(pred), test))
    preds, truth, sect = np.concatenate(preds), np.concatenate(truth), np.concatenate(sect)
    pooled_raw = np.array([_pearson(truth[:, j], preds[:, j]) for j in range(len(genes))])
    # centre per section: remove each section's mean from truth and prediction
    tc, pc = truth.copy(), preds.copy()
    for sid in ids:
        m = sect == sid
        tc[m] -= tc[m].mean(0); pc[m] -= pc[m].mean(0)
    pooled_centred = np.array([_pearson(tc[:, j], pc[:, j]) for j in range(len(genes))])
    # within-section permutation null of the predictions (per gene): the 95th percentile
    null95 = np.zeros(len(genes))
    for j in range(len(genes)):
        vals = []
        for _ in range(N_PERM):
            perm = pc[:, j].copy()
            for sid in ids:
                m = np.where(sect == sid)[0]
                perm[m] = perm[rng.permutation(m)]
            vals.append(_pearson(tc[:, j], perm))
        null95[j] = np.nanpercentile(vals, 95)
    med_within = np.nanmedian(np.vstack(list(per_section.values())), axis=0)
    df = pd.DataFrame({"gene": genes, "median_within_section_r": med_within, "pooled_r": pooled_raw,
                       "pooled_centred_r": pooled_centred, "null95_centred": null95,
                       "beats_null": pooled_centred > null95})
    return df, per_section, alphas


def run():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(C.SEED_SPLITS)
    secs = [s for s in S.load_all() if s.image is not None and (FEAT / f"{s.id}.npz").exists()]
    if SMOKE:
        secs = secs[:3]
    assert len(secs) >= 3, "the histology arm needs at least three sections with features"
    cand = pd.read_csv(STUDY2 / "results__dossier__candidates.csv")
    genes, cand_genes, ctrl_genes = targets(secs, cand.gene.unique())
    Y = expression(secs, genes)
    print(f"  {len(secs)} sections, {len(genes)} target genes ({len(cand_genes)} candidates, {len(ctrl_genes)} controls)")

    results = {}
    for kind, use_pca in (("phikon", True), ("imagenet", True), ("composition", False)):
        df, per_section, alphas = evaluate(kind, secs, genes, Y, rng, use_pca)
        df["is_candidate"] = df.gene.isin(cand_genes); df["is_control"] = df.gene.isin(ctrl_genes)
        df.to_csv(OUT / f"genes_{kind}.csv", index=False)
        pd.DataFrame(per_section, index=genes).to_csv(OUT / f"per_section_{kind}.csv")
        results[kind] = df
        c = df[df.is_candidate]
        print(f"  {kind:<12} candidates: median within-section r {c.median_within_section_r.median():.3f}, "
              f"pooled {c.pooled_r.median():.3f}, centred {c.pooled_centred_r.median():.3f}; "
              f"{int(c.beats_null.sum())}/{len(c)} beat the null; alphas {sorted(set(alphas.values()))}", flush=True)

    # compartment scores for P8
    comp = pd.read_csv(C.RESULTS / "compartment" / "candidates.csv").groupby("gene").lfc_tumour.median()
    ph = results["phikon"]; im = results["imagenet"]; co = results["composition"]
    cph = ph[ph.is_candidate].set_index("gene"); cim = im[im.is_candidate].set_index("gene"); cco = co[co.is_candidate].set_index("gene")
    ctrl_comp = ph[ph.gene.isin(C.CONTROLS["tumour"] + C.CONTROLS["immune"] + C.CONTROLS["stroma"])]
    med_raw = float(cph.median_within_section_r.median()); med_centred = float(cph.pooled_centred_r.median())
    ctrl_med = float(ctrl_comp.median_within_section_r.median())
    d_imagenet = float((cph.median_within_section_r - cim.median_within_section_r.reindex(cph.index)).median())
    d_comp = float((cph.median_within_section_r - cco.median_within_section_r.reindex(cph.index)).median())
    both = cph.join(comp.rename("lfc_tumour"), how="inner").dropna()
    rho = float(spearmanr(both.median_within_section_r, both.lfc_tumour.abs()).statistic)
    V = {"P6": {"holds": bool(med_raw < 0.3 and med_centred < 0.15 and ctrl_med > 0.5),
                "rule": "candidates: median within-section r < 0.3 and median centred r < 0.15; compartment controls > 0.5",
                "median_within_r": med_raw, "median_centred_r": med_centred, "controls_median_r": ctrl_med},
         "P7": {"holds": bool(d_imagenet >= C.SEOI_PEARSON and d_comp < C.SEOI_PEARSON),
                "rule": f"Phikon minus ImageNet >= {C.SEOI_PEARSON} and Phikon minus composition < {C.SEOI_PEARSON} (median over candidates)",
                "phikon_minus_imagenet": d_imagenet, "phikon_minus_composition": d_comp},
         "P8": {"holds": bool(rho >= 0.5), "rule": "Spearman(candidate r, |tumour LFC|) >= 0.5", "spearman": rho,
                "n": int(len(both))}}
    for k, v in V.items():
        print(f"  {k}: {'holds' if v['holds'] else 'fails'}  " + json.dumps({a: round(b, 3) for a, b in v.items() if isinstance(b, float)}))
    summary = {"verdicts": V, "sections": [s.id for s in secs], "n_genes": len(genes), "n_candidates": len(cand_genes),
               "n_perm": N_PERM, "controls": ph[ph.is_control][["gene", "median_within_section_r", "pooled_centred_r"]].to_dict("records")}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    if SMOKE:
        print("smoke run; nothing logged"); return
    P.log_run("histology", {"pca_dim": C.PCA_DIM, "alphas": list(ALPHAS), "n_perm": N_PERM, "seed": C.SEED_SPLITS},
              [str(p.relative_to(C.ROOT)) for p in sorted(OUT.glob("*"))], time.time() - t0)
    print(f"done in {(time.time() - t0) / 60:.0f} min")


if __name__ == "__main__":
    sys.exit(run())
