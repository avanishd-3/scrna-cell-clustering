"""Step 3: global embedding + clustering of all QC-passing cells.

Decisions:
- log-normalisation (CP10k + log1p); HVGs = seurat_v3 on raw counts, batch-aware (library).
- HVG exclusions: immunoglobulin / TCR V/D/J/C segment genes (clonotype-specific; would
  cluster cells by receptor gene usage instead of state), mito, ribosomal, sex-chromosome
  genes (donor identity), and hashtag-unrelated technical genes (MALAT1, NEAT1).
  Cell-cycle genes are KEPT: proliferation is a real state in BM progenitors / pre-B cells.
- PCA (50 PCs) -> Harmony on `library` (22 independent 10x reactions in 2 runs; run and donor
  are nested in library). Composition differs by enrichment (BM/CD34/B), which Harmony tolerates
  because it corrects within soft clusters rather than globally aligning distributions.
- kNN (k=20) on Harmony PCs -> Leiden at a range of resolutions; resolution chosen by
  seed-stability (mean pairwise ARI across 5 seeds) and by marker interpretability.
"""
import sys, re
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import harmonypy
sys.path.insert(0, str(Path(__file__).parent)); from knn import NNDTransformer
from sklearn.metrics import adjusted_rand_score

IN, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
NPC = int(sys.argv[3]) if len(sys.argv) > 3 else 30
sc.settings.figdir = OUT; sc.settings.set_figure_params(dpi=90, frameon=False)

adata = sc.read_h5ad(IN)
print(adata)

EXCL = re.compile(r"^(IG[HKL][VDJC]|IGHM$|IGHD$|IGHG|IGHA|IGHE|IGLL5$|TR[ABDG][VDJC]|MT-|RP[SL]\d|"
                  r"XIST$|TSIX$|RPS4Y|DDX3Y$|UTY$|KDM5D$|EIF1AY$|USP9Y$|ZFY$|MALAT1$|NEAT1$)")
# NB: Ig constant genes are excluded too (light-chain isotype / ambient plasma-cell Ig would
# dominate); IGLL1 (surrogate light chain, non-clonal) is kept. All still used for markers.
sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor="seurat_v3", batch_key="library", span=0.5)
hv = adata.var.highly_variable & ~adata.var_names.str.match(EXCL)
hv_genes = adata.var.loc[hv].sort_values("highly_variable_rank").index[:3000]
adata.var["hvg_used"] = adata.var_names.isin(hv_genes)
print("HVGs used:", adata.var.hvg_used.sum())

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # log-normalised, all genes, for markers

X = adata[:, adata.var.hvg_used].X.toarray().astype(np.float32)
mu, sd = X.mean(0), X.std(0); sd[sd == 0] = 1
X -= mu; X /= sd; np.clip(X, -10, 10, out=X)  # in place (RAM)
from sklearn.decomposition import PCA
pca = PCA(n_components=50, svd_solver="randomized", random_state=0)
adata.obsm["X_pca"] = pca.fit_transform(X).astype(np.float32)
adata.uns["pca_var_ratio"] = pca.explained_variance_ratio_
del X
np.savetxt(OUT / "pca_variance_ratio.txt", pca.explained_variance_ratio_)

ho = harmonypy.run_harmony(adata.obsm["X_pca"][:, :NPC], adata.obs, ["library"], random_state=0, max_iter_harmony=20)
Z = np.asarray(ho.Z_corr)
if Z.shape[0] != adata.n_obs: Z = Z.T
adata.obsm["X_harmony"] = np.ascontiguousarray(Z, dtype=np.float32)

sc.pp.neighbors(adata, n_neighbors=20, use_rep="X_harmony", random_state=0, transformer=NNDTransformer(20))
sc.tl.umap(adata, random_state=0)

stab = []
for res in [0.2, 0.4, 0.6, 0.8, 1.0, 1.5]:
    labs = []
    for seed in range(5):
        sc.tl.leiden(adata, resolution=res, key_added="tmp", random_state=seed, flavor="igraph", n_iterations=2)
        labs.append(adata.obs.tmp.values.copy())
    aris = [adjusted_rand_score(labs[i], labs[j]) for i in range(5) for j in range(i + 1, 5)]
    adata.obs[f"leiden_{res}"] = pd.Categorical(labs[0])
    stab.append(dict(resolution=res, n_clusters=len(set(labs[0])), mean_ARI=np.mean(aris), min_ARI=np.min(aris)))
    print(stab[-1], flush=True)
pd.DataFrame(stab).to_csv(OUT / "leiden_stability.csv", index=False)
del adata.obs["tmp"]

# lightweight output: drop X (keep raw) to save disk; keep obsm/obs
adata.write_h5ad(OUT / "global.h5ad", compression="gzip")
