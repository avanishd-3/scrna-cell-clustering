"""Step 5: re-cluster each top-level compartment at finer resolution.

usage: 05_subcluster.py global.h5ad counts.h5ad compartments.json outdir

For each compartment (a set of global clusters):
  * HVGs recomputed within the compartment (seurat_v3 on counts, batch=library) with the same
    exclusions as globally -> the subspace reflects within-compartment variation.
  * PCA (30) -> Harmony(library) -> kNN(k=20) -> Leiden.
  * Resolution: sweep, keep the highest resolution whose mean pairwise ARI across 5 Leiden seeds
    is >= 0.75 (stability criterion; boundaries along continua shift between seeds).
  * Over-splitting guard: iteratively merge a cluster into its nearest cluster (Harmony centroid
    distance) if the pair has < MIN_DE genes with log2FC > 1, padj < 0.01 and detection
    difference >= 20 percentage points in either direction.
"""
import sys, json, re
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, harmonypy
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
sys.path.insert(0, str(Path(__file__).parent)); from knn import NNDTransformer

G, C, CJ, OUT = sys.argv[1], sys.argv[2], sys.argv[3], Path(sys.argv[4]); OUT.mkdir(parents=True, exist_ok=True)
COMP = json.load(open(CJ))      # {"B_lineage": ["8","9",...], ...}
GKEY = COMP.pop("_key", "leiden_0.2")
RES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2]
ARI_MIN, MIN_DE = 0.75, 10
EXCL = re.compile(r"^(IG[HKL][VDJC]|IGHM$|IGHD$|IGHG|IGHA|IGHE|IGLL5$|TR[ABDG][VDJC]|MT-|RP[SL]\d|"
                  r"XIST$|TSIX$|RPS4Y|DDX3Y$|UTY$|KDM5D$|EIF1AY$|USP9Y$|ZFY$|MALAT1$|NEAT1$)")

g = ad.read_h5ad(G, backed="r")
gobs = g.obs[[GKEY]].copy()
cnt = ad.read_h5ad(C)
EXCLUDE_LIBS = {"Bcell_2_2"}   # compromised library: labelled later by kNN transfer (06)
gobs = gobs[~g.obs.library.isin(EXCLUDE_LIBS).values]
cnt = cnt[gobs.index].copy()
cnt.obs[GKEY] = gobs[GKEY].astype(str).values

def n_de(a, c1, c2, key):
    """# genes distinguishing c1 vs c2 (both directions): log2FC>1, padj<0.01, detection diff>=20pp."""
    s = a[a.obs[key].isin([c1, c2])]
    sc.tl.rank_genes_groups(s, key, groups=[c1], reference=c2, method="wilcoxon")
    d = sc.get.rank_genes_groups_df(s, c1).set_index("names")
    X = s.X; in1 = (s.obs[key] == c1).values
    p1 = pd.Series(np.asarray((X[in1] > 0).mean(0)).ravel(), index=s.var_names)
    p2 = pd.Series(np.asarray((X[~in1] > 0).mean(0)).ravel(), index=s.var_names)
    dp = (p1 - p2).reindex(d.index)
    up = (d.logfoldchanges > 1) & (d.pvals_adj < 0.01) & (dp >= 0.2)
    dn = (d.logfoldchanges < -1) & (d.pvals_adj < 0.01) & (dp <= -0.2)
    return int(up.sum() + dn.sum())

all_labels = []; report = []
for name, clusters in COMP.items():
    if (OUT / f"{name}.h5ad").exists():   # cache: compartment already done
        all_labels.append(ad.read_h5ad(OUT / f"{name}.h5ad", backed="r").obs[["subcluster"]].copy())
        print(f"== {name}: cached", flush=True); continue
    a = cnt[cnt.obs[GKEY].isin(clusters)].copy()
    print(f"== {name}: {a.n_obs} cells from global clusters {clusters}", flush=True)
    try:
        sc.pp.highly_variable_genes(a, n_top_genes=3000, flavor="seurat_v3", batch_key="library", span=0.5)
    except ValueError:  # loess singularity in small compartments / tiny per-library batches
        sc.pp.highly_variable_genes(a, n_top_genes=3000, flavor="seurat_v3", batch_key="library", span=1.0)
    hv = a.var.highly_variable & ~a.var_names.str.match(EXCL)
    hv_genes = a.var.loc[hv].sort_values("highly_variable_rank").index[:2000]
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    X = a[:, hv_genes].X.toarray().astype(np.float32)
    X -= X.mean(0); sd = X.std(0); sd[sd == 0] = 1; X /= sd; np.clip(X, -10, 10, out=X)
    pcs = PCA(30, svd_solver="randomized", random_state=0).fit_transform(X); del X
    ho = harmonypy.run_harmony(pcs, a.obs, ["library"], random_state=0, max_iter_harmony=20)
    Z = np.asarray(ho.Z_corr); Z = Z if Z.shape[0] == a.n_obs else Z.T
    a.obsm["X_harmony"] = np.ascontiguousarray(Z, dtype=np.float32)
    sc.pp.neighbors(a, n_neighbors=20, use_rep="X_harmony", transformer=NNDTransformer(20))
    chosen = None
    for res in RES:
        labs = []
        for seed in range(5):
            sc.tl.leiden(a, resolution=res, key_added="tmp", random_state=seed, flavor="igraph", n_iterations=2)
            labs.append(a.obs.tmp.astype(str).values.copy())
        ari = np.mean([adjusted_rand_score(labs[i], labs[j]) for i in range(5) for j in range(i + 1, 5)])
        k = len(set(labs[0])); print(f"  res={res} k={k} ARI={ari:.3f}", flush=True)
        report.append(dict(compartment=name, resolution=res, k=k, mean_ARI=round(ari, 3)))
        if ari >= ARI_MIN and k > 1: chosen = (res, labs[0])
    res, lab = chosen if chosen else (RES[0], labs[0])
    a.obs["sub"] = lab
    # merge clusters without distinctive markers vs nearest neighbour cluster
    while True:
        cl = sorted(a.obs["sub"].unique(), key=lambda x: int(x) if x.isdigit() else x)
        if len(cl) < 2: break
        cen = np.stack([a.obsm["X_harmony"][a.obs["sub"].values == c].mean(0) for c in cl])
        D = np.linalg.norm(cen[:, None] - cen[None], axis=2); np.fill_diagonal(D, np.inf)
        merged = False
        for i in np.argsort(D.min(1)):
            j = int(np.argmin(D[i])); nd = n_de(a, cl[i], cl[j], "sub")
            if nd < MIN_DE:
                print(f"  merge {cl[i]} -> {cl[j]} (DE genes={nd})", flush=True)
                a.obs.loc[a.obs["sub"] == cl[i], "sub"] = cl[j]; merged = True; break
        if not merged: break
    # relabel by size
    order = a.obs["sub"].value_counts().index
    m = {c: str(i) for i, c in enumerate(order)}
    a.obs["sub"] = a.obs["sub"].map(m)
    sc.tl.umap(a, random_state=0)
    a.obs["subcluster"] = name + "_" + a.obs["sub"]
    print(f"  chosen res={res}; final k={a.obs['sub'].nunique()}", flush=True)
    a.raw = a
    keep = a[:, a.var.highly_variable | a.var_names.isin(hv_genes)].copy()
    keep.write_h5ad(OUT / f"{name}.h5ad", compression="gzip")
    all_labels.append(a.obs[["subcluster"]])

rp = OUT / "subcluster_resolution_sweep.csv"
pd.concat(([pd.read_csv(rp)] if rp.exists() else []) + [pd.DataFrame(report)]).drop_duplicates().to_csv(rp, index=False)
pd.concat(all_labels).to_csv(OUT / "subcluster_labels.csv")
