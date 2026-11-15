"""Evaluate a clustering: UMAPs, canonical-marker dotplot, per-cluster QC/composition, DE markers.
usage: 04_evaluate.py <h5ad> <cluster_key> <outdir> [embedding]
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

H5, KEY, OUT = sys.argv[1], sys.argv[2], Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = OUT; sc.settings.set_figure_params(dpi=80, frameon=False)
a = sc.read_h5ad(H5)
a.obs[KEY] = a.obs[KEY].astype(str).astype("category")

MARKERS = {
    "HSPC": ["CD34", "CRHBP", "AVP", "HLF", "PROM1", "SPINK2"],
    "LMPP/CLP": ["FLT3", "DNTT", "CD7", "ACY3"],
    "proB": ["IGLL1", "VPREB1", "EBF1", "PAX5", "CD19", "MME", "RAG1", "RAG2", "CD79A", "CD79B"],
    "cycling": ["MKI67", "TOP2A", "CDK1"],
    "imm/mat B": ["MS4A1", "IGHM", "IGHD", "TCL1A", "CD24", "FCER2", "CR2", "CD27", "TNFRSF13B"],
    "plasma": ["JCHAIN", "MZB1", "XBP1", "PRDM1", "SDC1"],
    "T/NK": ["CD3E", "CD4", "CD8A", "IL7R", "NKG7", "GNLY"],
    "myeloid": ["MPO", "ELANE", "AZU1", "LTF", "CAMP", "CD14", "LYZ", "FCGR3A", "CSF1R"],
    "DC": ["FCER1A", "CLEC10A", "LILRA4", "IRF8"],
    "ery/mk/baso": ["GATA1", "KLF1", "HBB", "CA1", "PF4", "ITGA2B", "CPA3", "HDC"],
    "stroma": ["CXCL12", "LEPR", "COL1A1"],
}
MARKERS = {k: [g for g in v if g in a.raw.var_names] for k, v in MARKERS.items()}

for c in [KEY, "library", "enrichment", "run", "hash_call"]:
    if c in a.obs:
        sc.pl.umap(a, color=c, legend_loc="on data" if c == KEY else "right margin",
                   legend_fontsize=7, s=1, show=False, save=f"_{KEY}_{c}.png")
sc.pl.umap(a, color=["doublet_score", "log1p_total_counts", "pct_counts_mt", "n_genes_by_counts"],
           s=1, show=False, save=f"_{KEY}_qc.png", vmax="p99")
sc.pl.dotplot(a, MARKERS, groupby=KEY, standard_scale="var", show=False, save=f"{KEY}_markers.png")

o = a.obs
summ = o.groupby(KEY, observed=True).agg(n=("library", "size"), med_umi=("total_counts", "median"),
        med_genes=("n_genes_by_counts", "median"), mt=("pct_counts_mt", "median"),
        dbl=("doublet_score", "mean"))
comp = pd.crosstab(o[KEY], o.enrichment, normalize="index").round(2)
nlib = pd.crosstab(o[KEY], o.library)
summ["max_lib_frac"] = (nlib.max(1) / nlib.sum(1)).round(2)
summ["n_libs_ge1pct"] = (nlib.div(nlib.sum(1), axis=0) > 0.01).sum(1)
# DE on a subsample (≤1500 cells / cluster) for speed
rng = np.random.default_rng(0)
idx = np.concatenate([rng.choice(np.where(o[KEY].values == g)[0], min((o[KEY] == g).sum(), 1500), replace=False)
                      for g in o[KEY].cat.categories])
s = a[idx].copy()
sc.tl.rank_genes_groups(s, KEY, method="wilcoxon", use_raw=True, pts=True)
df = sc.get.rank_genes_groups_df(s, None)
df = df[(df.logfoldchanges > 0.5) & (df.pvals_adj < 0.01) & ~df.names.str.match(r"^(MT-|RP[SL])")]
top = df.groupby("group").head(12).groupby("group").names.apply(lambda x: ",".join(x))
summ = summ.join(comp).join(top.rename("top_markers"))
summ.to_csv(OUT / f"{KEY}_summary.csv")
df.groupby("group").head(50).to_csv(OUT / f"{KEY}_DE_top50.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print(summ.to_string())
