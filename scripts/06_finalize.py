"""Step 6: assemble final labels.

* Cells from the 21 reference libraries get their compartment sub-cluster label directly.
* Cells from the excluded library (Bcell_2_2) are labelled by kNN majority vote (k=30) against
  reference cells in the global Harmony space; vote fraction is reported as confidence.
* Writes results/cell_clusters.csv (one row per retained cell) and results/removed_cells.csv
  (cells removed by QC / doublet filters, with reason).
usage: 06_finalize.py global.h5ad sub_labels.csv annotations.csv qc_dir lineage_calls.csv outdir
"""
import sys, glob
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
from sklearn.neighbors import NearestNeighbors

G, SUB, ANN, QCD, LD, OUT = sys.argv[1:7]; OUT = Path(OUT)
g = ad.read_h5ad(G)
sub = pd.read_csv(SUB, index_col=0)["subcluster"]
ann = pd.read_csv(ANN, index_col=0)["annotation"] if Path(ANN).exists() else pd.Series(dtype=str)

o = g.obs[["library", "barcode", "enrichment", "run", "hash_call", "leiden_0.2"]].copy()
o["cluster"] = sub.reindex(o.index)
ref = o.cluster.notna().values
Z = g.obsm["X_harmony"]
nn = NearestNeighbors(n_neighbors=30).fit(Z[ref])
_, idx = nn.kneighbors(Z[~ref])
rl = np.asarray(o.cluster.values[ref], dtype=object)
votes = pd.DataFrame(rl[idx])
lab = np.asarray(votes.mode(axis=1)[0].values, dtype=object)
conf = (votes.values == lab[:, None]).mean(1)
o.loc[~ref, "cluster"] = lab
o["label_source"] = np.where(ref, "clustered", "knn_transfer")
o["transfer_confidence"] = np.nan
o.loc[~ref, "transfer_confidence"] = conf.round(3)
o["qc_flag"] = np.where(o.library == "Bcell_2_2", "compromised_library_ambient_or_doublet_contamination", "")
o["compartment"] = o.cluster.str.rsplit("_", n=1).str[0]
o["annotation"] = o.cluster.map(ann) if len(ann) else ""
o["donor"] = o.hash_call.str.replace("Cellhash", "donor")
o.index.name = "cell_id"
o = o.rename(columns={"leiden_0.2": "global_cluster"})
cols = ["cluster", "annotation", "compartment", "global_cluster", "library", "barcode", "enrichment",
        "donor", "label_source", "transfer_confidence", "qc_flag"]
o[cols].to_csv(OUT / "cell_clusters.csv")
print(o.cluster.value_counts().sort_index().to_string())
print("transfer confidence (Bcell_2_2): median", np.median(conf), "frac>=0.6", (conf >= 0.6).mean())

# removed cells with reasons
rem = []
for f in glob.glob(f"{QCD}/obs_*.csv.gz"):
    d = pd.read_csv(f, index_col=0)
    reason = np.select([~d.qc_pass, d.hash_call == "Doublet", d.scrublet_doublet],
                       ["qc_fail", "hashtag_doublet", "scrublet_doublet"], "")
    d["removal_reason"] = reason
    rem.append(d.loc[reason != "", ["library", "barcode", "removal_reason"]])
rem = pd.concat(rem)
ldc = pd.read_csv(LD, index_col=0)
ldc = ldc[ldc.lineage_doublet]
rem = pd.concat([rem, pd.DataFrame({"library": ldc.index.str.rsplit("_", n=1).str[0],
                                    "barcode": ldc.index.str.rsplit("_", n=1).str[1],
                                    "removal_reason": "lineage_coexpression_doublet"}, index=ldc.index)])
rem.index.name = "cell_id"
rem.to_csv(OUT / "removed_cells.csv")
print(rem.removal_reason.value_counts())

# final UMAP figure
g.obs["final_cluster"] = o.cluster.astype("category")
sc.settings.figdir = OUT / "figures"; sc.settings.set_figure_params(dpi=100, frameon=False)
sc.pl.umap(g, color="final_cluster", s=1, legend_fontsize=6, show=False, save="_final_clusters.png")
g.obs["compartment"] = o.compartment.astype("category")
sc.pl.umap(g, color="compartment", s=1, show=False, save="_compartments.png")
