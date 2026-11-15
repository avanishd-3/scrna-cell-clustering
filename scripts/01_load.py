"""Step 1: load each 10x library, call cells for unfiltered matrices, save per-library h5ad.

Some GEO matrices are Cell Ranger *filtered* (~10-25k barcodes), others are *raw*
(~730k barcodes, incl. empty droplets). For raw ones we call cells ourselves with a
conservative knee/inflection rule on the total-UMI barcode-rank curve (GEX counts only),
with a hard floor of 500 UMIs. Downstream QC + clustering removes residual junk.
"""
import sys, gzip, json
from pathlib import Path
import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp
import anndata as ad

RAW = Path(sys.argv[1]); OUT = Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
libs = sorted({p.name.replace("GSE337533_", "").rsplit("_", 1)[0] for p in RAW.glob("*_matrix.mtx.gz")})
summary = []

def knee_threshold(umis, floor=500):
    """Inflection point of log(rank) vs log(UMI), searched among barcodes above `floor`.
    Same idea as DropletUtils::barcodeRanks inflection."""
    u = np.sort(umis)[::-1]
    u = u[u >= 50]
    r = np.arange(1, len(u) + 1)
    lx, ly = np.log10(r), np.log10(u)
    # smooth on a uniform log-rank grid, take most negative derivative
    grid = np.linspace(lx[0], lx[-1], 2000)
    yi = np.interp(grid, lx, ly)
    k = 15
    ys = np.convolve(yi, np.ones(k) / k, mode="same")
    d = np.gradient(ys, grid)
    # restrict to the region after the first 1000 barcodes (avoid top plateau artifacts)
    valid = (grid > np.log10(1000)) & (grid < lx[-1] - 0.2)
    i = np.argmin(np.where(valid, d, np.inf))
    thr = 10 ** yi[i]
    return max(thr, floor)

for lib in libs:
    f = lambda s: RAW / f"GSE337533_{lib}_{s}"
    feats = pd.read_csv(f("features.tsv.gz"), sep="\t", header=None, names=["id", "name", "type"])
    bcs = pd.read_csv(f("barcodes.tsv.gz"), header=None)[0].values
    X = scipy.io.mmread(f("matrix.mtx.gz")).T.tocsr().astype(np.float32)  # cells x features
    gex = (feats.type == "Gene Expression").values
    umis = np.asarray(X[:, gex].sum(1)).ravel()
    is_raw = len(bcs) > 100_000
    if is_raw:
        thr = knee_threshold(umis)
        keep = umis >= thr
    else:
        thr = float("nan"); keep = np.ones(len(bcs), bool)
    X = X[keep]
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=[f"{lib}_{b}" for b in bcs[keep]]),
                   var=pd.DataFrame({"gene_id": feats.id.values, "feature_type": feats.type.values},
                                    index=pd.Index(feats.name.values)))
    a.var_names_make_unique()
    a.obs["library"] = lib
    a.obs["barcode"] = bcs[keep]
    a.obs["enrichment"] = lib.split("_")[0]
    a.write_h5ad(OUT / f"{lib}.h5ad", compression="gzip")
    s = dict(library=lib, input_barcodes=len(bcs), raw_matrix=is_raw, umi_threshold=thr,
             cells_kept=int(keep.sum()), median_umi=float(np.median(umis[keep])))
    print(s, flush=True); summary.append(s)

pd.DataFrame(summary).to_csv(OUT / "load_summary.csv", index=False)
