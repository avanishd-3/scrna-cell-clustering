"""Diagnostics: barcode-rank curves for raw libraries + barcode overlap between libraries."""
import sys, itertools
from pathlib import Path
import numpy as np, pandas as pd, scipy.io, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, anndata as ad
RAW, LIBS, OUT = map(Path, sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
summ = pd.read_csv(LIBS/"load_summary.csv")
raw = summ[summ.raw_matrix]
fig, axs = plt.subplots(2, 4, figsize=(18, 8))
for ax, (_, r) in zip(axs.ravel(), raw.iterrows()):
    lib = r.library
    feats = pd.read_csv(RAW/f"GSE337533_{lib}_features.tsv.gz", sep="\t", header=None)
    X = scipy.io.mmread(RAW/f"GSE337533_{lib}_matrix.mtx.gz").tocsr()
    g = (feats[2]=="Gene Expression").values
    u = np.sort(np.asarray(X[g].sum(0)).ravel())[::-1]; u = u[u>0]
    ax.loglog(np.arange(1,len(u)+1), u); ax.axhline(r.umi_threshold, color="r")
    ax.set_title(f"{lib}: thr={r.umi_threshold:.0f}, n={r.cells_kept}")
    np.save(OUT/f"rank_{lib}.npy", u[:200000])
plt.tight_layout(); plt.savefig(OUT/"barcode_rank_raw_libs.png", dpi=80)
# barcode overlap (strip -1 suffix)
bc = {l: set(ad.read_h5ad(LIBS/f"{l}.h5ad", backed="r").obs.barcode.str[:16]) for l in summ.library}
rows=[]
for a,b in itertools.combinations(bc,2):
    o=len(bc[a]&bc[b]); rows.append((a,b,o, o/min(len(bc[a]),len(bc[b]))))
df=pd.DataFrame(rows,columns=["a","b","overlap","frac_of_smaller"]).sort_values("frac_of_smaller",ascending=False)
df.to_csv(OUT/"barcode_overlap.csv",index=False); print(df.head(15).to_string())
