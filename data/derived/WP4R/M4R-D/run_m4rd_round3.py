#!/usr/bin/env python
"""M4R-D ROUND1c: D6 natural-experiment pilot. AP2/PUF live outside core-181 (reported
limitation); salvage = mafft pilot on wide-coverage kmerRBH OGs + 5 CHROM core OGs."""
import os, re, subprocess
from collections import Counter
import numpy as np
import pandas as pd

ROOT = "/home/huyudi/015_plasmo"
L2 = f"{ROOT}/data/derived/WP4/L2_site_composition"
ALN = f"{L2}/alignments"
OUT = f"{ROOT}/data/derived/WP4R/M4R-D"
MAFFT = os.path.expanduser("~/.conda/envs/sc/bin/mafft")
SPECIES = ["SP001","SP002","SP003","SP004","SP005","SP006","SP007","SP008",
           "SP011","SP012","SP013","SP014","SP015","SP016","SP017","SP018"]
LOGC = ["SP001","SP002","SP003"]; HIGC = ["SP004","SP005","SP006","SP007"]
GC = {"SP001":0.193,"SP002":0.186,"SP003":0.182,"SP004":0.398,"SP005":0.40,"SP006":0.40,"SP007":0.40}

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
    for s in range(0, L - w + 1, step):
        win = [c for c in arr[s:s + w] if c]
        if len(win) < 16: continue
        if win_ent(win) < thr: mask[s:s + w] = True
    return mask

def parse_aln(path):
    d = {}; hid, seq = None, []
    def flush():
        if hid is not None: d[hid] = "".join(seq)
    for line in open(path):
        if line.startswith(">"): flush(); hid = line[1:].split()[0]; seq = []
        else: seq.append(line.strip())
    flush(); return d

import hashlib
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

mani2 = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-02_aa/input_manifest.tsv", sep="\t", header=0)
faa = {}
for _, r in mani2.iterrows():
    m = re.match(r"protein_(SP\d+)_", str(r.iloc[0]))
    if m: faa[m.group(1)] = f"{ROOT}/{r.iloc[1]}"
PROT = {s: parse_aln(faa[s]) for s in SPECIES}
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

# candidate OGs: wide coverage kmerRBH OGs holding AP2/PUF Pf genes
cands = []
for name, col in [("AP2", "is_AP2"), ("PUF_RNA", "is_PUF_RNA")]:
    for g in PF[PF[col]].gene:
        og = pid2og.get(pf2xp.get(g))
        if og is None: continue
        nsp = OG[OG.orthogroup_id == og].species_id.nunique()
        cands.append((name, g, og, nsp))
cd = pd.DataFrame(cands, columns=["set", "gene", "og", "nsp"]).sort_values("nsp", ascending=False)
print(cd.head(12).to_string(), flush=True)
picks = []
for name in ["AP2", "PUF_RNA"]:
    sub = cd[cd["set"] == name].drop_duplicates("og").head(2)
    picks += sub.to_dict("records")
print("picks:", picks, flush=True)

def analyze_og(og, gene, lab, aln):
    spp = [s for s in SPECIES if s in aln]
    L = len(aln[spp[0]])
    anchor = nn = 0
    for i in range(L):
        colaa = [aln[s][i] for s in spp if aln[s][i] != "-"]
        if len(colaa) >= max(8, int(0.6 * len(spp))):
            nn += 1
            if max(Counter(colaa).values()) / len(colaa) >= 0.8: anchor += 1
    feats = {}
    for s in spp:
        sq = aln[s].replace("-", "")
        if not sq: continue
        feats[s] = (float(lcr_mask(aln[s]).mean()), sq.count("N") / len(sq), len(sq))
    lo = [feats[s] for s in LOGC if s in feats]; hi = [feats[s] for s in HIGC if s in feats]
    return {"gene": gene, "set": lab, "og": og, "n_spp": len(spp), "alen": L,
            "frac_anchor_cols": (anchor / nn if nn else float("nan")),
            "lowGCpole_LCR": float(np.mean([x[0] for x in lo])) if lo else float("nan"),
            "highGCpole_LCR": float(np.mean([x[0] for x in hi])) if hi else float("nan"),
            "lowGCpole_Asn": float(np.mean([x[1] for x in lo])) if lo else float("nan"),
            "highGCpole_Asn": float(np.mean([x[1] for x in hi])) if hi else float("nan"),
            "pole_LCR_diff": (float(np.mean([x[0] for x in lo]) - np.mean([x[0] for x in hi])) if lo and hi else float("nan")),
            "pole_Asn_diff": (float(np.mean([x[1] for x in lo]) - np.mean([x[1] for x in hi])) if lo and hi else float("nan")),
            "aln_source": "mafft_pilot" if og.startswith("PILOT") else "core181"}

rows = []
for p in picks:
    og = p["og"]
    sub = OG[(OG.orthogroup_id == og)]
    faa_in = f"{OUT}/tmp_{og}.faa"; aln_out = f"{OUT}/tmp_{og}.aln"
    with open(faa_in, "w") as f:
        for r in sub.itertuples():
            s = r.species_id
            if s in PROT and r.protein_id in PROT[s]:
                f.write(f">{s}\n{PROT[s][r.protein_id]}\n")
    env = dict(os.environ); env.pop("LD_PRELOAD", None)
    rr = subprocess.run([MAFFT, "--auto", "--quiet", faa_in], capture_output=True, text=True, env=env, timeout=600)
    assert rr.returncode == 0, rr.stderr[:500]
    open(aln_out, "w").write(rr.stdout)
    aln = parse_aln(aln_out)
    rows.append(analyze_og("PILOT_" + og, p["gene"], p["set"], aln))
    os.remove(faa_in)
# 5 CHROM core OGs
for g, og in [('PF3D7_0320900', 'OG_003406'), ('PF3D7_0617800', 'OG_005217'),
              ('PF3D7_0714000', 'OG_001577'), ('PF3D7_0925700', 'OG_004073'),
              ('PF3D7_1105100', 'OG_000397')]:
    aln = parse_aln(f"{ALN}/{og}.aln")
    rows.append(analyze_og(og, g, "CHROM", aln))
NAT = pd.DataFrame(rows)
NAT.to_csv(f"{OUT}/M4RD_D6_natural_experiment.tsv", sep="\t", index=False)
print(NAT.to_string(), flush=True)
files = sorted([f for f in os.listdir(OUT) if f.startswith("M4RD_")])
pd.DataFrame({"file": files, "sha256": [sha(f"{OUT}/{f}") for f in files]}).to_csv(
    f"{OUT}/input_manifest.tsv", sep="\t", index=False)
print("D6 done:", len(NAT), flush=True)
