#!/usr/bin/env python
"""M4R-A ROUND2 step3: lifecycle tripole (gametocyte + zygote + liver) with BH."""
import gzip, os, re
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4R/M4R-A"

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)

def mwu(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3 or np.unique(np.concatenate([a, b])).size < 10:
        return float("nan"), float("nan"), float("nan"), float("nan")
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic")
    return float(u), float(p), float(np.median(a)), float(np.median(b))

# ---- universe + features (same recipe) ----
FAA = f"{ROOT}/data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7_AnnotatedProteins.fasta"
def gkey(s):
    m = re.search(r"PF3D7_\d+", str(s)); return m.group(0) if m else None
prots = {}
with open(FAA) as fh:
    hid, seq = None, []
    def flush():
        if hid is None: return
        g = gkey(hid)
        if not g: return
        s = "".join(seq).upper().replace("*", "")
        if len(s) > len(prots.get(g, "")): prots[g] = s
    for line in fh:
        if line.startswith(">"): flush(); hid, seq = line.rstrip(), []
        else: seq.append(line.strip())
    flush()
def feats(s):
    n = len(s); asn = s.count("N") / n if n else 0.0
    mt, cur = 0, 0
    for ch in s:
        cur = cur + 1 if ch == "N" else 0
        mt = max(mt, cur)
    w, th, low = 20, 1.8, 0
    if n >= w:
        arr = np.frombuffer(s.encode(), dtype=np.uint8)
        cnt = np.zeros(256, int)
        for c in arr[:w]: cnt[c] += 1
        def ent():
            p = cnt[cnt > 0] / w
            return float(-(p * np.log2(p)).sum())
        low = int(ent() < th)
        for i in range(w, n):
            cnt[arr[i]] += 1; cnt[arr[i - w]] -= 1
            p = cnt[cnt > 0] / w
            if float(-(p * np.log2(p)).sum()) < th: low += 1
        lcr = low / (n - w + 1)
    else:
        lcr = 0.0
    return n, asn, mt, lcr
F = pd.DataFrame({g: feats(s) for g, s in prots.items()},
                 index=["length", "asn_frac", "polyn_max", "lcr_frac"]).T
F.index.name = "gene"; F = F.reset_index()
q90 = F.asn_frac.quantile(0.9)
F["asn_top10"] = F.asn_frac >= q90
print("universe:", len(F), "q90:", round(float(q90), 4), flush=True)

M = pd.read_csv(f"{OUT}/M4RA_interpro_map.tsv", sep="\t")
g2ipr = dict(zip(M.gene, M.interpro.fillna("").apply(lambda s: set(s.split(";")) if s else set())))
F["iprset"] = F.gene.map(lambda g: g2ipr.get(g, set()))
F["AP2_ip"] = F.iprset.map(lambda s: "IPR001471" in s)
F["PUF_strict_ip"] = F.iprset.map(lambda s: "IPR001313" in s)

def batt(lfc, pole, contrast, min_n=3):
    """Test battery on a per-gene logFC series. Returns rows."""
    d = pd.DataFrame({"lfc": lfc}).join(F.set_index("gene"), how="inner")
    d = d[np.isfinite(d.lfc)]
    rows = []
    v = d.asn_frac.values; y = d.lfc.values
    r, p = stats.spearmanr(v, y)
    rows.append({"pole": pole, "contrast": contrast, "test": "spearman(asn_frac,logFC)",
                 "n": len(d), "stat": float(r), "p": float(p), "med_hit": float("nan"),
                 "med_rest": float("nan")})
    u, p2, mh, mr = mwu(d[d.asn_top10].lfc, d[~d.asn_top10].lfc)
    rows.append({"pole": pole, "contrast": contrast, "test": "MW Asn-top10 vs rest",
                 "n": len(d), "stat": u, "p": p2, "med_hit": mh, "med_rest": mr})
    for col, nm in [("AP2_ip", "MW AP2(IPR001471) vs rest"),
                    ("PUF_strict_ip", "MW PUF-strict vs rest")]:
        h = d[d[col]].lfc; o = d[~d[col]].lfc
        if len(h) < min_n:
            rows.append({"pole": pole, "contrast": contrast, "test": nm, "n": len(d),
                         "stat": float("nan"), "p": float("nan"),
                         "med_hit": float(h.median()) if len(h) else float("nan"),
                         "med_rest": float(o.median())})
            continue
        u3, p3, mh3, mr3 = mwu(h, o)
        rows.append({"pole": pole, "contrast": contrast, "test": nm, "n": len(d),
                     "stat": u3, "p": p3, "med_hit": mh3, "med_rest": mr3})
    return rows

ALL = []
# ---- pole 1: gametocyte (reuse round-1 counts) ----
G = pd.read_csv(f"{OUT}/M4RA_gametocyte_counts.tsv", sep="\t")
G["expressed"] = (G.male_cpm + G.female_cpm) > 2
Ge = G[G.expressed]
print("gametocyte expressed:", len(Ge), flush=True)
ALL += batt(Ge.set_index("gene").log2FC_FvsM, "gametocyte", "female-vs-male")

# ---- pole 2: zygote T0/T2/T4 ----
Z = pd.read_csv(f"{ROOT}/data/raw/geo/GSE222586/GSE222586_zygote_raw_expressionmatrix.txt.gz",
                sep=";", index_col=0)
print("zygote matrix:", Z.shape, flush=True)
meta = pd.read_csv(f"{ROOT}/data/raw/geo/GSE222586/GSE222586_metadata.txt.gz", sep=";",
                   index_col=0)
print("zygote meta Timepoints:", meta.Timepoint.value_counts().to_dict(), flush=True)
meta = meta[meta.index.isin(Z.columns)]
Z = Z.loc[:, meta.index]
tp = meta.Timepoint
pb = {}
for t, cols in tp.groupby(tp).groups.items():
    pb[t] = Z.loc[:, list(cols)].sum(axis=1)
PB = pd.DataFrame(pb)
PB = PB.loc[(PB > 0).any(axis=1)]
lib = PB.sum(axis=0)
CP = np.log2(PB.div(lib / 1e6, axis=1) + 0.5)
ts = sorted(PB.columns)
print("zygote pseudobulk:", PB.shape, "groups:", ts, flush=True)
for a, b, nm in [(ts[1], ts[0], f"{ts[1]}-vs-{ts[0]}"), (ts[2], ts[0], f"{ts[2]}-vs-{ts[0]}"),
                 (ts[2], ts[1], f"{ts[2]}-vs-{ts[1]}")] if len(ts) >= 3 else []:
    ALL += batt(CP[a] - CP[b], "zygote", nm)

# ---- pole 3: liver Pf Day0/2/4/6 GFP/NoGFP ----
LP = f"{ROOT}/data/raw/geo/GSE220039/GSE220039_PfHsMmu_GRCh38_GRCm39_Pf3D7v58_raw_counts.csv.gz"
cols_needed = None
pf_parts = []
for ch in pd.read_csv(LP, chunksize=20000):
    ch = ch[ch.gene_name.str.startswith("PF3D7_", na=False)]
    if cols_needed is None:
        cols_needed = [c for c in ch.columns if c != "gene_name"]
    pf_parts.append(ch)
L = pd.concat(pf_parts, ignore_index=True)
print("liver Pf rows:", len(L), flush=True)
L = L.groupby("gene_name", as_index=False).sum(numeric_only=True)
L = L.set_index("gene_name")
lib2 = L.sum(axis=0)
CL = np.log2(L.div(lib2 / 1e6, axis=1) + 0.5)
gfp = [c for c in CL.columns if ".GFP_" in c and "NoGFP" not in c]
nog = [c for c in CL.columns if ".NoGFP_" in c]
d0 = [c for c in CL.columns if c.startswith("Day0.")]
print("liver GFP:", len(gfp), "NoGFP:", len(nog), "Day0:", len(d0), flush=True)
if gfp and nog:
    ALL += batt(CL[gfp].mean(axis=1) - CL[nog].mean(axis=1), "liver", "GFP-vs-NoGFP(all-days)")
for day in ["Day2.", "Day4.", "Day6."]:
    gg = [c for c in gfp if c.startswith(day)]; nn = [c for c in nog if c.startswith(day)]
    if gg and nn:
        ALL += batt(CL[gg].mean(axis=1) - CL[nn].mean(axis=1), "liver", f"GFP-vs-NoGFP({day.rstrip('.')})")
if d0 and gfp:
    ALL += batt(CL[gfp].mean(axis=1) - CL[d0].mean(axis=1), "liver", "GFP-vs-Day0naive")

T = pd.DataFrame(ALL)
T["q_bh"] = bh(T.p.fillna(1.0).values)
T.to_csv(f"{OUT}/M4RA_lifecycle_tripole.tsv", sep="\t", index=False)
print(T.to_string(), flush=True)
print("STEP3 DONE", flush=True)
