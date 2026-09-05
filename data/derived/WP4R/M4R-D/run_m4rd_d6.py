#!/usr/bin/env python
"""M4R-D ROUND2 D6-formal: pre-registered transition contrasts on wide AP2/PUF/CHROM OGs.
Transitions (pre-registered sets; direction expectation from composition polarity):
  T1_AT_Lav:    A={SP001,SP002} (low-GC) vs B={SP003}              expect A>B (Asn)
  T2_GC_vivax:  A={SP004..007} (high-GC) vs B={SP008,SP011-013}    expect A<B
  T3_rodent:    A={SP011-013} vs B={SP008}                         expect per GC polarity (checked at runtime)
Support rule: >=2 transitions with composition-consistent direction (sign test per transition).
mafft timeout 300s, max 2 retries, then gene=UNRESOLVED. Writes M4RD_D6formal_*.tsv."""
import os, re, subprocess, sys
from collections import Counter
import numpy as np
import pandas as pd

ROOT = "/home/huyudi/015_plasmo"
L2 = f"{ROOT}/data/derived/WP4/L2_site_composition"
OUT = f"{ROOT}/data/derived/WP4R/M4R-D"
MAFFT = os.path.expanduser("~/.conda/envs/sc/bin/mafft")
TRANS = {"T1_AT_Lav": (["SP001", "SP002"], ["SP003"], 1),
         "T2_GC_vivax": (["SP004", "SP005", "SP006", "SP007"], ["SP008", "SP011", "SP012", "SP013"], -1),
         "T3_rodent": (["SP011", "SP012", "SP013"], ["SP008"], None)}

def parse_aln(path):
    d = {}; hid, seq = None, []
    def flush():
        if hid is not None: d[hid] = "".join(seq)
    for line in open(path):
        if line.startswith(">"): flush(); hid = line[1:].split()[0]; seq = []
        else: seq.append(line.strip())
    flush(); return d

_AA = 'ARNDCQEGHILKMFPSTWYV'
_AI = {a: i for i, a in enumerate(_AA)}
def _ent_counts(cnt, tot):
    if tot == 0: return 0.0
    if np.count_nonzero(cnt) >= 8: return 2.0
    p = cnt[cnt > 0] / tot
    return float(-np.sum(p * np.log2(p)))
def lcr_mask(seq, w=64, step=8, thr=1.5):
    L = len(seq); mask = np.zeros(L, bool)
    arr = [c if c != '-' else None for c in seq]
    def win_ent(win):
        cnt = np.zeros(20)
        for c in win:
            j = _AI.get(c)
            if j is not None: cnt[j] += 1
        return _ent_counts(cnt, len(win))
    if L < w:
        sub = [c for c in arr if c]
        if sub and win_ent(sub) < thr: mask[:] = True
        return mask
    for st in range(0, L - w + 1, step):
        win = arr[st:st + w]
        real = [c for c in win if c]
        if real and win_ent(real) < thr:
            mask[st:st + w] = True
    return mask

mani2 = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-02_aa/input_manifest.tsv", sep="\t", header=0)
faa = {}
for _, r in mani2.iterrows():
    m = re.match(r"protein_(SP\d+)_", str(r.iloc[0]))
    if m: faa[m.group(1)] = f"{ROOT}/{r.iloc[1]}"
PROT = {}
for s, fp in faa.items():
    if os.path.exists(fp): PROT[s] = parse_aln(fp)
    else: print(f"  skip missing proteome {s}: {fp}", flush=True)
print(f"loaded proteomes: {sorted(PROT)}", flush=True)
OG = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-02_aa/orthogroups.tsv", sep="\t")
PF = pd.read_csv(f"{OUT}/M4RD_D5_pf_features.tsv", sep="\t")
xp2pf = {}
gff = f"{ROOT}/data/raw/ncbi-datasets/GCF_000002765.6/ncbi_dataset/data/GCF_000002765.6/genomic.gff"
for line in open(gff):
    if "\tCDS\t" not in line: continue
    m1 = re.search(r"ID=cds-(XP_\d+\.\d+)", line); m2 = re.search(r"locus_tag=(PF3D7_\d+)", line)
    if m1 and m2: xp2pf[m1.group(1)] = m2.group(1)
pf2xp = {v: k for k, v in xp2pf.items()}
pid2og = dict(OG[OG.species_id == "SP001"].set_index("protein_id").orthogroup_id)
m1 = pd.read_csv(f"{ROOT}/data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv", sep="\t")
spcol = m1.columns[0]
gcol = [c for c in m1.columns if "genome" in c.lower() and "gc" in c.lower()][0]
GC = {s: float(m1[m1[spcol].astype(str) == s][gcol].iloc[0]) for s in
      ["SP001","SP002","SP003","SP004","SP005","SP006","SP007","SP008","SP011","SP012","SP013"]}
print("GC:", {k: round(v,3) for k,v in GC.items()}, flush=True)
# resolve T3 expectation from GC polarity
t3a = np.mean([GC[s] for s in TRANS["T3_rodent"][0]])
t3b = float(GC["SP008"])
TRANS["T3_rodent"] = (TRANS["T3_rodent"][0], TRANS["T3_rodent"][1], 1 if t3a < t3b else -1)
print(f"T3 polarity: rodent {t3a:.3f} vs malariae {t3b:.3f} -> expect sign {TRANS['T3_rodent'][2]}", flush=True)

cands = []
for name, col in [("AP2", "is_AP2"), ("PUF_RNA", "is_PUF_RNA"), ("CHROM", "is_CHROM")]:
    for g in PF[PF[col]].gene:
        og = pid2og.get(pf2xp.get(g))
        if og is None: continue
        nsp = OG[OG.orthogroup_id == og].species_id.nunique()
        cands.append((name, g, og, nsp))
cd = pd.DataFrame(cands, columns=["set", "gene", "og", "nsp"]).sort_values("nsp", ascending=False)
print(cd.groupby("set").nsp.describe().to_string(), flush=True)
picks = []
for name in ["AP2", "PUF_RNA"]:
    sub = cd[(cd["set"] == name) & (cd.nsp >= 8)].drop_duplicates("og")
    picks += sub.to_dict("records")
ch = cd[(cd["set"] == "CHROM") & (cd.nsp >= 8)].drop_duplicates("og").head(6)
picks += ch.to_dict("records")
print(f"picks: {len(picks)} OGs", flush=True)

def run_mafft(og):
    sub = OG[OG.orthogroup_id == og]
    faa_in = f"{OUT}/tmp_d6_{og}.faa"; aln_out = f"{OUT}/tmp_d6_{og}.aln"
    with open(faa_in, "w") as f:
        for r in sub.itertuples():
            s = r.species_id
            if s in PROT and r.protein_id in PROT[s]:
                f.write(f">{s}\n{PROT[s][r.protein_id]}\n")
    env = dict(os.environ); env.pop("LD_PRELOAD", None)
    for att in range(2):
        try:
            rr = subprocess.run([MAFFT, "--auto", "--quiet", faa_in],
                                capture_output=True, text=True, env=env, timeout=300)
            if rr.returncode == 0:
                open(aln_out, "w").write(rr.stdout)
                os.remove(faa_in)
                return parse_aln(aln_out)
        except subprocess.TimeoutExpired:
            print(f"  mafft timeout {og} att{att}", flush=True)
    try: os.remove(faa_in)
    except OSError: pass
    return None

rows = []
for p in picks:
    og, gene, lab = p["og"], p["gene"], p["set"]
    aln = run_mafft(og)
    if aln is None:
        for t in TRANS:
            rows.append({"gene": gene, "set": lab, "og": og, "transition": t,
                         "n_A": 0, "n_B": 0, "mean_Asn_A": np.nan, "mean_Asn_B": np.nan,
                         "diff_Asn": np.nan, "diff_LCR": np.nan, "status": "UNRESOLVED_mafft"})
        continue
    feats = {}
    for s, sq in aln.items():
        raw = sq.replace("-", "")
        if not raw: continue
        feats[s] = (float(lcr_mask(sq).mean()), raw.count("N") / len(raw))
    for t, (A, B, exp) in TRANS.items():
        a = [feats[s] for s in A if s in feats]; b = [feats[s] for s in B if s in feats]
        if len(a) < 2 or len(b) < 1:
            rows.append({"gene": gene, "set": lab, "og": og, "transition": t,
                         "n_A": len(a), "n_B": len(b), "mean_Asn_A": np.nan,
                         "mean_Asn_B": np.nan, "diff_Asn": np.nan, "diff_LCR": np.nan,
                         "status": "UNRESOLVED_coverage"})
            continue
        da = float(np.mean([x[1] for x in a]) - np.mean([x[1] for x in b]))
        dl = float(np.mean([x[0] for x in a]) - np.mean([x[0] for x in b]))
        rows.append({"gene": gene, "set": lab, "og": og, "transition": t,
                     "n_A": len(a), "n_B": len(b),
                     "mean_Asn_A": float(np.mean([x[1] for x in a])),
                     "mean_Asn_B": float(np.mean([x[1] for x in b])),
                     "diff_Asn": da, "diff_LCR": dl, "status": "OK"})
R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/M4RD_D6formal_pertransition.tsv", sep="\t", index=False)
# per-transition sign summary (Asn), sign test toward expected direction
from scipy import stats as st
summ = []
for (t, lab), grp in R[R.status == "OK"].groupby(["transition", "set"]):
    exp = TRANS[t][2]
    d = grp.diff_Asn.values
    npos = int((np.sign(d) == exp).sum()); n = len(d)
    p = float(st.binomtest(npos, n, 0.5, alternative="greater").pvalue) if n else np.nan
    summ.append({"transition": t, "set": lab, "n_genes": n,
                 "median_diff_Asn": float(np.median(d)) if n else np.nan,
                 "n_expected_sign": npos, "signtest_p": p,
                 "expected_sign": exp})
S = pd.DataFrame(summ)
S.to_csv(f"{OUT}/M4RD_D6formal_summary.tsv", sep="\t", index=False)
print(S.to_string(), flush=True)
print(f"D6 done: {len(R)} rows, {(R.status=='OK').sum()} OK", flush=True)
