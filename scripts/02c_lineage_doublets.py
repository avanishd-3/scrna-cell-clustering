"""Step 2c: heterotypic doublet filter based on co-expression of mutually exclusive lineage programs.

Motivation (see README): after the first global clustering, library Bcell_2_2 (largest library,
33k cells; hashtag demux unreliable because of very high Cellhash4 ambient) contributed a discrete
population co-expressing B-cell genes (CD79A/B, VPREB1, IGHM, MS4A1) with granulocytic genes
(MPO, ELANE, AZU1, LYZ). Such cells were ~0-1% in all other libraries. Scrublet missed them.

Rule: score = mean log-normalised expression of each program; a cell is flagged if any two
mutually exclusive programs both exceed 1.0. IGLL1 and CD24 are deliberately excluded from the
B score (expressed in myeloid progenitors / neutrophils too).
usage: 02c_lineage_doublets.py merged_qc.h5ad out.h5ad report_dir
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp

IN, OUTH5, REP = sys.argv[1], sys.argv[2], Path(sys.argv[3]); REP.mkdir(parents=True, exist_ok=True)
PROGRAMS = {
    "B": ["CD79A", "CD79B", "VPREB1", "VPREB3", "CD19", "IGHM", "MS4A1", "CD22"],
    "My": ["MPO", "ELANE", "AZU1", "PRTN3", "CTSG", "LYZ", "S100A8", "S100A9"],
    "Ery": ["HBB", "HBA1", "HBA2", "HBD", "CA1", "AHSP"],
    "TNK": ["CD3E", "CD3D", "NKG7", "GNLY", "CCL5", "IL7R"],
}
THR = 1.0
a = sc.read_h5ad(IN)
lib = a.X.sum(1).A1
for k, genes in PROGRAMS.items():
    X = a[:, genes].X
    X = X.multiply(1e4 / lib[:, None]).tocsr(); X.data = np.log1p(X.data)
    a.obs[f"score_{k}"] = np.asarray(X.mean(1)).ravel()
hi = {k: a.obs[f"score_{k}"].values > THR for k in PROGRAMS}
pairs = [("B", "My"), ("B", "Ery"), ("B", "TNK"), ("My", "Ery"), ("My", "TNK"), ("Ery", "TNK")]
flag = np.zeros(a.n_obs, bool)
for p, q in pairs:
    f = hi[p] & hi[q]; a.obs[f"conf_{p}_{q}"] = f; flag |= f
a.obs["lineage_doublet"] = flag
rep = a.obs.groupby("library", observed=True)[[f"conf_{p}_{q}" for p, q in pairs] + ["lineage_doublet"]].mean().round(4)
rep["n"] = a.obs.library.value_counts()
rep.to_csv(REP / "lineage_doublet_rates.csv"); print(rep.to_string())
a.obs[["lineage_doublet"] + [f"score_{k}" for k in PROGRAMS]].to_csv(REP / "lineage_doublet_calls.csv.gz")
a = a[~flag].copy()
print("kept", a.n_obs, "removed", int(flag.sum()))
a.write_h5ad(OUTH5, compression="gzip")
