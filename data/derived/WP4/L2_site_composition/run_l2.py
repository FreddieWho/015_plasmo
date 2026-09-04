#!/usr/bin/env python
# L2 v1: composition-aware site-level analysis (D-018/D-024)
# 181 single-copy core OGs x 16 spp; MAFFT MSA; per-site composition coupling vs genome GC;
# severity by property-class change; H-L2a near-neutral churn; H-L2b known adaptive sites off-axis test
import hashlib, json, os, re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
import numpy as np
import pandas as pd

ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4/L2_site_composition"
ALN = f"{OUT}/alignments"
os.makedirs(ALN, exist_ok=True)
SEED = 20260904
rng = np.random.default_rng(SEED)
NPERM = 500
MAFFT = os.path.expanduser("~/.conda/envs/sc/bin/mafft")
SPECIES = ["SP001","SP002","SP003","SP004","SP005","SP006","SP007","SP008",
           "SP011","SP012","SP013","SP014","SP015","SP016","SP017","SP018"]

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

# ---------- 0. inputs ----------
mani2 = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-02_aa/input_manifest.tsv", sep="\t", header=0)
faa = {}
for _, r in mani2.iterrows():
    m = re.match(r"protein_(SP\d+)_", str(r.iloc[0]))
    if m: faa[m.group(1)] = f"{ROOT}/{r.iloc[1]}"
assert set(SPECIES) <= set(faa), f"missing faa: {set(SPECIES)-set(faa)}"

m1 = pd.read_csv(f"{ROOT}/data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv", sep="\t")
spcol = m1.columns[0]
gcol = [c for c in m1.columns if "genome" in c.lower() and "gc" in c.lower()][0]
GC = dict(zip(m1[spcol].astype(str), m1[gcol].astype(float)))
gcv = np.array([GC[s] for s in SPECIES])
print("GC per species:", dict(zip(SPECIES, gcv.round(3))), flush=True)

OG = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-02_aa/orthogroups.tsv", sep="\t")
CORE = [l.strip() for l in open(f"{ROOT}/data/derived/WP2/M2-02_aa/single_copy_core.list") if l.strip()]
print("core OGs:", len(CORE), flush=True)

def parse_faa(path):
    d = {}
    with open(path) as fh:
        hid, seq = None, []
        def flush():
            if hid is not None: d[hid] = "".join(seq)
        for line in fh:
            if line.startswith(">"):
                flush(); hid = line[1:].split()[0]; seq = []
            else: seq.append(line.strip())
        flush()
    return d

PROT = {s: parse_faa(faa[s]) for s in SPECIES}
print("proteins loaded:", {s: len(v) for s, v in list(PROT.items())[:3]}, "...", flush=True)

# ---------- 1. build per-OG fasta + mafft ----------
def run_og(og):
    sub = OG[OG.orthogroup_id == og]
    recs = []
    for s in SPECIES:
        row = sub[sub.species_id == s]
        if len(row) != 1: return og, "skip_not_single"
        pid = row.iloc[0].protein_id
        seq = PROT[s].get(pid)
        if not seq: return og, "skip_no_seq"
        recs.append((s, seq))
    fa = f"{ALN}/{og}.faa"; aln = f"{ALN}/{og}.aln"
    with open(fa, "w") as fh:
        for s, sq in recs: fh.write(f">{s}\n{sq}\n")
    if not os.path.exists(aln):
        r = subprocess.run([MAFFT, "--auto", fa], capture_output=True, text=True)
        if r.returncode != 0: return og, "mafft_fail"
        open(aln, "w").write(r.stdout)
    return og, "ok"

with ThreadPoolExecutor(max_workers=8) as ex:
    res = list(ex.map(run_og, CORE))
stat = Counter(r[1] for r in res)
print("align stats:", dict(stat), flush=True)

def parse_aln(path):
    d = {}
    hid, seq = None, []
    def flush():
        if hid is not None: d[hid] = "".join(seq)
    for line in open(path):
        if line.startswith(">"): flush(); hid = line[1:].split()[0]; seq = []
        else: seq.append(line.strip())
    flush()
    return d

# ---------- 2. site stats ----------
GC_CLASS = set("ARGPWV"); AT_CLASS = set("NDIKFY")   # M2-02 GC-coupled classes
CLS = {a: (1 if a in GC_CLASS else -1 if a in AT_CLASS else 0) for a in "ARNDCQEGHILKMFPSTWYV"}
# property classes for severity
PROP = {}
for a in "AG": PROP[a] = "tiny"
for a in "VILM": PROP[a] = "hydrophobic"
for a in "FWY": PROP[a] = "aromatic"
for a in "KRH": PROP[a] = "positive"
for a in "DE": PROP[a] = "negative"
for a in "NQ": PROP[a] = "amide"
for a in "ST": PROP[a] = "hydroxyl"
PROP["C"] = "cys"; PROP["P"] = "proline"
ADJ = {frozenset(["tiny", "hydroxyl"]), frozenset(["tiny", "hydrophobic"]),
       frozenset(["amide", "hydroxyl"]), frozenset(["amide", "positive"]),
       frozenset(["negative", "amide"]), frozenset(["aromatic", "hydrophobic"]),
       frozenset(["hydroxyl", "hydrophobic"])}
def severity(a, b):
    if a == b: return 0
    if PROP.get(a) == PROP.get(b): return 0          # conservative (same class)
    if frozenset([PROP.get(a, "?"), PROP.get(b, "?")]) in ADJ: return 1  # moderate
    return 2                                          # radical

site_rows = []
aln_cache = {}
for og, st in res:
    if st != "ok": continue
    aln = parse_aln(f"{ALN}/{og}.aln")
    if set(SPECIES) - set(aln): continue
    aln_cache[og] = aln
    L = len(aln["SP001"])
    mat = np.array([[CLS.get(aln[s][i], np.nan) for s in SPECIES] for i in range(L)], dtype=float)
    aas = [[aln[s][i] for s in SPECIES] for i in range(L)]
    for i in range(L):
        col = mat[i]; ok = np.isfinite(col) & (col != 0)
        uniq = set(x for x in aas[i] if x != "-")
        if len(uniq) < 2 or ok.sum() < 10: continue
        r = np.corrcoef(col[ok], gcv[ok])[0, 1]
        sev = [severity(a, b) for a in uniq for b in uniq if a < b]
        site_rows.append({"og": og, "col": i, "n_class": int(ok.sum()), "r_gc": r,
                          "mean_severity": float(np.mean(sev)) if sev else 0.0,
                          "max_severity": max(sev) if sev else 0,
                          "aas": "".join(sorted(uniq))})
SITES = pd.DataFrame(site_rows)
print("variable classed sites:", len(SITES), flush=True)

# parametric p for correlation (t on r, df=n_class-2) — fast and standard at n=10..16
from scipy import stats as _st
n_c = SITES.n_class.values.astype(float)
df_ = np.maximum(n_c - 2, 1)
r_obs = SITES.r_gc.values
tstat = np.abs(r_obs) * np.sqrt(df_ / np.maximum(1 - r_obs ** 2, 1e-12))
SITES["p_perm"] = 2 * _st.t.sf(tstat, df_)  # column kept for schema compat; parametric
p = SITES.p_perm.values; okp = np.isfinite(p)
q = np.full(len(p), np.nan)
po = p[okp]; n = len(po); o = np.argsort(po)
qq = np.empty(n); qq[o] = np.minimum.accumulate((po[o] * n / np.arange(1, n + 1))[::-1])[::-1]
q[okp] = np.minimum(qq, 1)
SITES["q_perm"] = q
SITES.to_csv(f"{OUT}/L2_site_table.tsv", sep="\t", index=False)
print("site table written", flush=True)

# ---------- 3. H-L2a: coupled sites more conservative? ----------
coupled = SITES[SITES.q_perm < 0.05]
other = SITES[(SITES.q_perm >= 0.1)]
if len(coupled) >= 20 and len(other) >= 100:
    a = coupled.mean_severity.values; b = other.mean_severity.values
    obs = np.median(a) - np.median(b)
    pool = np.concatenate([a, b]); na = len(a); cnt = 0
    for i in range(1000):
        idx = rng.permutation(len(pool))
        if np.median(pool[idx[:na]]) - np.median(pool[idx[na:]]) <= obs: cnt += 1
    p_a = (cnt + 1) / 1001
else:
    obs, p_a = np.nan, np.nan
h2a = {"n_coupled": len(coupled), "n_other": len(other),
       "median_severity_coupled": float(np.median(coupled.mean_severity)) if len(coupled) else None,
       "median_severity_other": float(np.median(other.mean_severity)) if len(other) else None,
       "frac_conservative_coupled": float((coupled.mean_severity == 0).mean()) if len(coupled) else None,
       "frac_conservative_other": float((other.mean_severity == 0).mean()) if len(other) else None,
       "median_diff": None if np.isnan(obs) else float(obs), "p_perm": None if np.isnan(p_a) else float(p_a)}

# ---------- 4. H-L2b: known adaptive sites vs coupling ----------
KNOWN = [  # Pf 3D7 gene, positions, drug
    ("PF3D7_1343700", [580, 539], "K13 artemisinin"),
    ("PF3D7_0709000", [76], "CRT chloroquine"),
    ("PF3D7_0417200", [51, 59, 108, 164], "DHFR pyrimethamine"),
    ("PF3D7_0810800", [436, 437, 540, 581], "DHPS sulfadoxine"),
    ("PF3D7_0523000", [86, 184, 1246], "MDR1 multi"),
]
# map PF3D7 -> protein_id in SP001 faa headers
pf_hdr = {}
with open(faa["SP001"]) as fh:
    for line in fh:
        if line.startswith(">"):
            m = re.search(r"\[locus_tag=(PF3D7_\d+)\]", line)
            if m: pf_hdr[m.group(1)] = line[1:].split()[0]
kb = []
for g, positions, drug in KNOWN:
    pid = pf_hdr.get(g)
    og_rows = OG[(OG.species_id == "SP001") & (OG.protein_id == pid)] if pid else []
    if len(og_rows) == 0:
        kb.append({"gene": g, "drug": drug, "status": "not_in_kmerRBH_OG"}); continue
    og = og_rows.iloc[0].orthogroup_id
    if og not in aln_cache:
        kb.append({"gene": g, "drug": drug, "og": og, "status": "not_in_core181"}); continue
    aln = aln_cache[og]
    # map Pf sequence position -> alignment column
    seq = aln["SP001"]
    pos2col = {}
    c = 0
    for i, ch in enumerate(seq):
        if ch != "-": c += 1; pos2col[c] = i
    for pos in positions:
        col = pos2col.get(pos)
        if col is None: continue
        srow = SITES[(SITES.og == og) & (SITES.col == col)]
        kb.append({"gene": g, "drug": drug, "og": og, "pf_pos": pos, "aln_col": col,
                   "status": "mapped",
                   "r_gc": float(srow.r_gc.iloc[0]) if len(srow) else np.nan,
                   "q_perm": float(srow.q_perm.iloc[0]) if len(srow) else np.nan,
                   "variable": bool(len(srow))})
KB = pd.DataFrame(kb)
KB.to_csv(f"{OUT}/L2_known_sites.tsv", sep="\t", index=False)

# ---------- 5. deliverables ----------
qc = {"align_stats": dict(stat), "n_variable_sites": int(len(SITES)),
      "n_coupled_q05": int(len(coupled)), "H2a": h2a,
      "known_sites_mapped": int((KB.status == "mapped").sum()) if len(KB) else 0,
      "seed": SEED, "nperm": NPERM}
json.dump(qc, open(f"{OUT}/L2_qc.json", "w"), indent=2)
open(f"{OUT}/claim_impact.md", "w").write(f"""# L2 claim impact (v1, kmerRBH core OGs; exploratory D-019/D-022)

- date: 2026-09-04
- H-L2a (near-neutral churn): coupled n={h2a['n_coupled']} median severity {h2a['median_severity_coupled']} vs other {h2a['median_severity_other']}, permutation p={h2a['p_perm']}
- H-L2b (known adaptive sites): see L2_known_sites.tsv (mapped={qc['known_sites_mapped']})
- claim impact: none yet; exploratory. Official-orthogroup rerun planned when download lands (D-024).
- limitation: severity = property-class change (Grantham swap-in later); k-mer tree/OGs from M2 fallback; 16-spp n small per site.
""")
mani = pd.DataFrame({"file": [f"{ROOT}/data/derived/WP2/M2-02_aa/orthogroups.tsv",
                              f"{ROOT}/data/derived/WP2/M2-02_aa/single_copy_core.list",
                              f"{ROOT}/data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv"]})
mani["sha256"] = mani.file.map(sha); mani["size"] = mani.file.map(os.path.getsize)
mani.to_csv(f"{OUT}/input_manifest.tsv", sep="\t", index=False)
outs = [f for f in os.listdir(OUT) if f not in ("checksums.sha256", "alignments") and os.path.isfile(os.path.join(OUT, f))]
with open(f"{OUT}/checksums.sha256", "w") as fh:
    for f in sorted(outs):
        fh.write(f"{sha(os.path.join(OUT, f))}  {f}\n")
print("DONE", flush=True)
