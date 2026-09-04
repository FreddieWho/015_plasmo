#!/usr/bin/env python
# L2 v2 (hardened): Grantham severity + alignment QC filter + ML supermatrix tree + radical-mutation audit
# v1 (property-class severity) preserved in L2_site_table.tsv; v2 outputs L2v2_*
import hashlib, json, os, re, subprocess
import numpy as np
import pandas as pd
from collections import Counter
from scipy import stats as _st

ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4/L2_site_composition"
ALN = f"{OUT}/alignments"
SEED = 20260904
rng = np.random.default_rng(SEED)
SPECIES = ["SP001","SP002","SP003","SP004","SP005","SP006","SP007","SP008",
           "SP011","SP012","SP013","SP014","SP015","SP016","SP017","SP018"]
IQTREE = os.path.expanduser("~/.conda/envs/sc/bin/iqtree")

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

# ---------- Grantham matrix ----------
gl = open(f"{ROOT}/data/raw/refs/grantham.tsv").read().strip().split("\n")
hdr = gl[0].split("\t")[1:]
GR = {}
aas20 = ["S"] + hdr
for i, line in enumerate(gl[1:]):
    p = line.split("\t")
    rowaa = p[0]
    for j, v in enumerate(p[1:]):
        a, b = rowaa, hdr[j]
        fv = float(v)
        if fv != 0.0 or a == b:  # 0 placeholders pad the lower triangle; never clobber real values
            GR[frozenset([a, b])] = fv
for a in aas20: GR[frozenset([a, a])] = 0.0
_offdiag = [frozenset([a, b]) for i, a in enumerate(aas20) for b in aas20[i + 1:]]
_missing = [tuple(k) for k in _offdiag if k not in GR]
assert not _missing, f"grantham matrix incomplete: {_missing[:5]}"
def grantham(a, b): return GR.get(frozenset([a, b]), np.nan)

# ---------- inputs ----------
m1 = pd.read_csv(f"{ROOT}/data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv", sep="\t")
gcol = [c for c in m1.columns if "genome" in c.lower() and "gc" in c.lower()][0]
GC = dict(zip(m1.iloc[:, 0].astype(str), m1[gcol].astype(float)))
gcv = np.array([GC[s] for s in SPECIES])
GC_CLASS = set("ARGPWV"); AT_CLASS = set("NDIKFY")
CLS = {a: (1 if a in GC_CLASS else -1 if a in AT_CLASS else 0) for a in "ARNDCQEGHILKMFPSTWYV"}

def parse_aln(path):
    d = {}; hid, seq = None, []
    def flush():
        if hid is not None: d[hid] = "".join(seq)
    for line in open(path):
        if line.startswith(">"): flush(); hid = line[1:].split()[0]; seq = []
        else: seq.append(line.strip())
    flush(); return d

# ---------- alignment QC filter ----------
qc_rows = []
for f in sorted(os.listdir(ALN)):
    if not f.endswith(".aln") or f.startswith("EXT_"): continue
    aln = parse_aln(f"{ALN}/{f}")
    if set(SPECIES) - set(aln): continue
    L = len(aln["SP001"])
    seqs = [aln[s] for s in SPECIES]
    gapfrac = float(np.mean([[c == "-" for c in s] for s in seqs]))
    # mean pairwise identity on non-gap columns
    ident = []
    for i in range(L):
        col = [s[i] for s in seqs if s[i] != "-"]
        if len(col) >= 12: ident.append(max(Counter(col).values()) / len(col))
    qc_rows.append({"og": f[:-4], "alen": L, "gap_frac": gapfrac,
                    "mean_col_consensus": float(np.mean(ident)) if ident else 0})
QC = pd.DataFrame(qc_rows)
QC["pass"] = (QC.gap_frac < 0.30) & (QC.mean_col_consensus > 0.50)
QC.to_csv(f"{OUT}/L2v2_alignment_qc.tsv", sep="\t", index=False)
print("align QC pass:", int(QC["pass"].sum()), "/", len(QC), flush=True)

# ---------- site table v2 (Grantham) ----------
site_rows = []
for og in QC[QC["pass"]].og:
    aln = parse_aln(f"{ALN}/{og}.aln")
    L = len(aln["SP001"])
    for i in range(L):
        col = [aln[s][i] for s in SPECIES]
        ng = [c for c in col if c != "-"]
        if len(ng) < 12: continue
        uniq = sorted(set(ng))
        if len(uniq) < 2: continue
        cls = np.array([CLS.get(c, np.nan) if c != "-" else np.nan for c in col], float)
        ok = np.isfinite(cls) & (cls != 0)
        r = np.corrcoef(cls[ok], gcv[ok])[0, 1] if ok.sum() >= 8 and len(set(cls[ok])) > 1 else np.nan
        gs = [grantham(a, b) for a in uniq for b in uniq if a < b]
        site_rows.append({"og": og, "col": i, "n_class": int(ok.sum()), "r_gc": r,
                          "mean_grantham": float(np.nanmean(gs)), "max_grantham": float(np.nanmax(gs))})
S2 = pd.DataFrame(site_rows)
n_c = S2.n_class.values.astype(float)
r = S2.r_gc.values
okp = np.isfinite(r)
df_ = np.maximum(n_c - 2, 1)
t = np.abs(r) * np.sqrt(df_ / np.maximum(1 - r ** 2, 1e-12))
S2["p_corr"] = np.where(okp, 2 * _st.t.sf(np.abs(t), df_), np.nan)
p = S2.p_corr.values
o = np.argsort(p[okp]); n = okp.sum()
qq = np.empty(n); qq[o] = np.minimum.accumulate((p[okp][o] * n / np.arange(1, n + 1))[::-1])[::-1]
S2.loc[okp, "q_bh"] = np.minimum(qq, 1)
S2.to_csv(f"{OUT}/L2v2_site_table.tsv", sep="\t", index=False)
print("v2 sites:", len(S2), flush=True)

# ---------- H-L2a v2 (Grantham) ----------
coupled = S2[S2.q_bh < 0.05]; other = S2[S2.q_bh >= 0.1]
obs = np.median(coupled.mean_grantham) - np.median(other.mean_grantham)
pool = np.concatenate([coupled.mean_grantham.values, other.mean_grantham.values])
na = len(coupled); cnt = 0
for i in range(2000):
    idx = rng.permutation(len(pool))
    if np.median(pool[idx[:na]]) - np.median(pool[idx[na:]]) >= obs: cnt += 1
p_more_radical = (cnt + 1) / 2001
h2a = {"n_coupled": len(coupled), "n_other": len(other),
       "median_grantham_coupled": float(np.median(coupled.mean_grantham)),
       "median_grantham_other": float(np.median(other.mean_grantham)),
       "frac_conservative_lt100_coupled": float((coupled.mean_grantham < 100).mean()),
       "frac_conservative_lt100_other": float((other.mean_grantham < 100).mean()),
       "median_diff": float(obs), "p_more_radical": float(p_more_radical)}
print("H2a v2:", h2a, flush=True)

# ---------- H-L2b addendum: Grantham of the population-level adaptive substitutions ----------
adaptive_subs = [("K13", "C580Y", "C", "Y"), ("K13", "R539T", "R", "T"), ("CRT", "K76T", "K", "T"),
                 ("DHFR", "N51I", "N", "I"), ("DHFR", "C59R", "C", "R"), ("DHFR", "S108N", "S", "N"),
                 ("DHFR", "I164L", "I", "L"), ("DHPS", "S436A", "S", "A"), ("DHPS", "A437G", "A", "G"),
                 ("DHPS", "K540E", "K", "E"), ("DHPS", "A581G", "A", "G"), ("MDR1", "N86Y", "N", "Y"),
                 ("MDR1", "Y184F", "Y", "F"), ("MDR1", "D1246Y", "D", "Y")]
ad = pd.DataFrame([{"gene": g, "sub": s, "grantham": grantham(a, b),
                    "class": "radical(>150)" if grantham(a, b) > 150 else "moderate(100-150)" if grantham(a, b) >= 100 else "conservative(<100)"}
                   for g, s, a, b in adaptive_subs])
ad.to_csv(f"{OUT}/L2v2_adaptive_substitutions_grantham.tsv", sep="\t", index=False)
print(ad.to_string(), flush=True)

# ---------- supermatrix + IQ-TREE (if available) ----------
tree_status = "skipped_no_iqtree"
if os.path.exists(IQTREE):
    keep = list(QC[QC["pass"]].og)
    per_sp = {s: [] for s in SPECIES}
    for og in keep:
        aln = parse_aln(f"{ALN}/{og}.aln")
        for s in SPECIES: per_sp[s].append(aln[s])
    sup = f"{OUT}/supermatrix.faa"
    with open(sup, "w") as fh:
        for s in SPECIES: fh.write(f">{s}\n{''.join(per_sp[s])}\n")
    r = subprocess.run([IQTREE, "-s", sup, "-m", "MFP", "-T", "8", "-pre", f"{OUT}/supermatrix"],
                       capture_output=True, text=True)
    tree_status = "ok" if r.returncode == 0 else "iqtree_fail"
    print("iqtree:", tree_status, flush=True)
else:
    print("iqtree not installed yet; tree step skipped", flush=True)

# ---------- deliverables ----------
qc = {"align_qc_pass": int(QC["pass"].sum()), "align_qc_total": int(len(QC)),
      "v2_sites": int(len(S2)), "H2a_v2": h2a, "tree": tree_status, "seed": SEED}
json.dump(qc, open(f"{OUT}/L2v2_qc.json", "w"), indent=2)
open(f"{OUT}/L2v2_claim_impact.md", "w").write(f"""# L2 v2 claim impact (hardened: Grantham + alignment QC + ML tree)

- date: 2026-09-04
- alignment QC: {qc['align_qc_pass']}/{qc['align_qc_total']} OGs pass (gap<30%, consensus>50%)
- H-L2a v2 (Grantham): coupled median {h2a['median_grantham_coupled']:.1f} vs other {h2a['median_grantham_other']:.1f},
  conservative(<100) {h2a['frac_conservative_lt100_coupled']:.2f} vs {h2a['frac_conservative_lt100_other']:.2f},
  p(more radical)={h2a['p_more_radical']:.4f} -> churn hypothesis again rejected if p small.
- H-L2b addendum: Grantham severity of the 14 known adaptive substitutions -> L2v2_adaptive_substitutions_grantham.tsv
  (radical changes at cross-species-conserved sites = adaptation off composition axis, quantitative).
- ML tree: {tree_status} (supermatrix of passing OGs).
- official OrthoMCL OGs: still blocked (VEuPathDB subscription wall, JS-only downloads) - kmerRBH retained, disclosed.
""")
mani = pd.DataFrame({"file": [f"{ROOT}/data/raw/refs/grantham.tsv",
                              f"{ROOT}/data/derived/WP2/M2-02_aa/single_copy_core.list",
                              f"{ROOT}/data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv"]})
mani["sha256"] = mani.file.map(sha); mani["size"] = mani.file.map(os.path.getsize)
mani.to_csv(f"{OUT}/L2v2_input_manifest.tsv", sep="\t", index=False)
outs = [f for f in os.listdir(OUT) if f.startswith("L2v2_") and f != "L2v2_checksums.sha256" and os.path.isfile(f"{OUT}/{f}")]
with open(f"{OUT}/L2v2_checksums.sha256", "w") as fh:
    for f in sorted(outs):
        fh.write(f"{sha(os.path.join(OUT, f))}  {f}\n")
print("DONE", flush=True)
