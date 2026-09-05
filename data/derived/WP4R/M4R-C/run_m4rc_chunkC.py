#!/usr/bin/env python3
"""Chunk runner for C-step DHA groups. Usage: chunkC.py START END (indices into sorted groups)."""
import sys
src = open("data/derived/WP4R/M4R-C/run_m4rc_salvage.py").read().split("STEP = ")[0]
exec(src)
import pandas as pd, numpy as np, time
t0 = time.time()
OUT = "data/derived/WP4R/M4R-C"
s, e = int(sys.argv[1]), int(sys.argv[2])
ALL = load_state(); have = set(r["contrast"] for r in ALL)
sv = pd.read_csv(f"{OUT}/M4RC_stagevec.tsv", sep="\t", index_col=0).iloc[:, 0].dropna()
RESP = pd.read_csv("data/derived/WP4/L1_codon_x_stress/L1_response_by_gene.tsv", sep="\t")
grps = sorted(RESP.groupby(["bg", "time"]).size().index.tolist())[s:e]
dha_rows = []
n = 0
for (bg, tm) in grps:
    j = RESP[(RESP.bg == bg) & (RESP.time == tm)]
    y = j.set_index("gene")["resp"].dropna(); y = y[y.index.isin(FEAT.index)]
    if len(y) < 100: continue
    lab = f"GSE151189_DHA-CTRL_{bg}_{tm}h"
    if lab not in have:
        ALL += ftests(y, lab); ALL += ftests(no_schiz(y), lab + "_noSchiz"); n += 1
    dha_rows.append(dict(contrast=lab, n=len(y), median_resp=float(y.median()),
                         median_noSchiz=float(no_schiz(y).median())))
    common = y.index.intersection(sv.index)
    if len(common) > 100:
        yy = y.loc[common].values; xx = sv.loc[common].values
        X = np.column_stack([np.ones(len(common)), xx])
        beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
        resid = pd.Series(yy - X @ beta, index=common)
        r2 = float(1 - np.var(yy - X @ beta) / np.var(yy))
        if lab + "_STAGEADJ" not in have:
            ALL += ftests(resid, lab + "_STAGEADJ")
            ALL += ftests(no_schiz(resid), lab + "_STAGEADJ_noSchiz"); n += 1
        dha_rows.append(dict(contrast=lab + "_STAGEADJ", n=len(common),
                             median_resp=float(resid.median()), stage_R2=r2,
                             median_noSchiz=float(no_schiz(resid).median())))
d = pd.DataFrame(dha_rows)
old = f"{OUT}/M4RC_stage_adjusted_DHA.tsv"
import os
if os.path.exists(old):
    prev = pd.read_csv(old, sep="\t")
    d = pd.concat([prev, d]).drop_duplicates("contrast", keep="last")
d.to_csv(old, sep="\t", index=False)
save_state(ALL); checkpoint("C_chunk", ALL)
print(f"CHUNK {s}:{e} new={n} T={round(time.time()-t0,1)}", flush=True)
