#!/usr/bin/env python
"""M4R-D ROUND1b: fix D6 OG mapping (GFF XP->PF3D7), D5 stratified counterfactual +
manual CMH, D4 Zhang-Pf constraint + Elsworth context. Fast follow-up; reuses round-1 files."""
import hashlib, os, re
from collections import Counter
import numpy as np
import pandas as pd
from scipy import stats as st

ROOT = "/home/huyudi/015_plasmo"
L2 = f"{ROOT}/data/derived/WP4/L2_site_composition"
ALN = f"{L2}/alignments"
OUT = f"{ROOT}/data/derived/WP4R/M4R-D"
SPECIES = ["SP001","SP002","SP003","SP004","SP005","SP006","SP007","SP008",
           "SP011","SP012","SP013","SP014","SP015","SP016","SP017","SP018"]
LOGC = ["SP001","SP002","SP003"]; HIGC = ["SP004","SP005","SP006","SP007"]

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

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
    for s in range(0, L - w + 1, step):
        win = [c for c in arr[s:s + w] if c]
        if len(win) < 16: continue
        if win_ent(win) < thr: mask[s:s + w] = True
    return mask

# ---------- XP -> PF3D7 map from RefSeq GFF ----------
xp2pf = {}
gff = f"{ROOT}/data/raw/ncbi-datasets/GCF_000002765.6/ncbi_dataset/data/GCF_000002765.6/genomic.gff"
with open(gff) as fh:
    for line in fh:
        if "\tCDS\t" not in line: continue
        m1 = re.search(r"ID=cds-(XP_\d+\.\d+)", line)
        m2 = re.search(r"locus_tag=(PF3D7_\d+)", line)
        if m1 and m2: xp2pf[m1.group(1)] = m2.group(1)
print("xp2pf:", len(xp2pf), flush=True)
OG = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-02_aa/orthogroups.tsv", sep="\t")
QC = pd.read_csv(f"{L2}/L2v2_alignment_qc.tsv", sep="\t")
PASS_OGS = set(o for o in QC[QC["pass"]].og if not str(o).startswith("EXT_"))
PF = pd.read_csv(f"{OUT}/M4RD_D5_pf_features.tsv", sep="\t")
pfset = {g: (bool(a), bool(p), bool(c)) for g, a, p, c in
         zip(PF.gene, PF.is_AP2, PF.is_PUF_RNA, PF.is_CHROM)}

# ---------- D6: natural experiment ----------
pid2og = dict(OG[OG.species_id == "SP001"].set_index("protein_id").orthogroup_id)
nat = []; mapped = {"AP2": 0, "PUF": 0}; total = {"AP2": 0, "PUF": 0}
for g, (a, p, c) in pfset.items():
    lab = "AP2" if a else ("PUF" if p else None)
    if lab is None: continue
    total[lab] += 1
    # find XP for this PF3D7
    xp = next((x for x, q in xp2pf.items() if q == g), None)
    og = pid2og.get(xp)
    if og is None or og not in PASS_OGS: continue
    if not os.path.exists(f"{ALN}/{og}.aln"): continue
    mapped[lab] += 1
    aln = parse_aln(f"{ALN}/{og}.aln")
    if set(SPECIES) - set(aln): continue
    L = len(aln[SPECIES[0]])
    anchor = 0; nn = 0
    for i in range(L):
        colaa = [aln[s][i] for s in SPECIES if aln[s][i] != "-"]
        if len(colaa) >= 12:
            nn += 1
            if max(Counter(colaa).values()) / len(colaa) >= 0.8: anchor += 1
    feats = {}
    for s in SPECIES:
        sq = aln[s].replace("-", "")
        if not sq: continue
        feats[s] = (float(lcr_mask(aln[s]).mean()), sq.count("N") / len(sq), len(sq))
    lo = [feats[s] for s in LOGC if s in feats]; hi = [feats[s] for s in HIGC if s in feats]
    if not lo or not hi: continue
    nat.append({"gene": g, "set": lab, "og": og, "alen": L,
                "frac_anchor_cols": (anchor / nn if nn else np.nan),
                "lowGCpole_LCR": float(np.mean([x[0] for x in lo])),
                "highGCpole_LCR": float(np.mean([x[0] for x in hi])),
                "lowGCpole_Asn": float(np.mean([x[1] for x in lo])),
                "highGCpole_Asn": float(np.mean([x[1] for x in hi])),
                "pole_LCR_diff": float(np.mean([x[0] for x in lo]) - np.mean([x[0] for x in hi])),
                "pole_Asn_diff": float(np.mean([x[1] for x in lo]) - np.mean([x[1] for x in hi]))})
NAT = pd.DataFrame(nat)
NAT.to_csv(f"{OUT}/M4RD_D6_natural_experiment.tsv", sep="\t", index=False)
print("D6 mapping:", total, mapped, "OG rows:", len(NAT), flush=True)
if len(NAT): print(NAT.to_string(), flush=True)

# ---------- D5b: stratified counterfactual (length tertile x polyN bin) + manual CMH ----------
PF["len_tert"] = pd.qcut(PF.length, 3, labels=["S", "M", "L"], duplicates="drop").astype(str)
PF["polyN_bin"] = pd.cut(PF.polyN_max, [-1, 0, 4, 1000], labels=["none", "short1-4", "tract5+"]).astype(str)
q90 = PF.asn_frac.quantile(0.9)
rows = []
for name, flag in [("AP2", PF.is_AP2), ("PUF_RNA", PF.is_PUF_RNA),
                   ("CHROM", PF.is_CHROM), ("regulator_any", PF.is_AP2 | PF.is_PUF_RNA | PF.is_CHROM)]:
    # within-stratum MW: regulator vs background Asn
    ds = []
    for (l, b), idx in PF.groupby(["len_tert", "polyN_bin"]).groups.items():
        g = PF.loc[idx]
        f = flag.loc[idx]
        if f.sum() >= 3 and ((~f).sum()) >= 10:
            u, p = st.mannwhitneyu(g[f].asn_frac, g[(~f)].asn_frac, alternative="two-sided")
            ds.append({"stratum": f"{l}x{b}", "n_reg": int(f.sum()), "n_bg": int((~f).sum()),
                       "median_reg": float(g[f].asn_frac.median()),
                       "median_bg": float(g[(~f)].asn_frac.median()),
                       "MW_p": float(p)})
    DS5 = pd.DataFrame(ds)
    DS5.to_csv(f"{OUT}/M4RD_D5_stratified_{name}.tsv", sep="\t", index=False)
    # manual CMH for asn_rich (top decile) across strata
    num = den = chisq_n = chisq_d = 0.0
    T = []
    for (l, b), idx in PF.groupby(["len_tert", "polyN_bin"]).groups.items():
        g = PF.loc[idx]; f = flag.loc[idx]
        a = int((f & (g.asn_frac >= q90)).sum()); b_ = int(((~f) & (g.asn_frac >= q90)).sum())
        c = int((f & (g.asn_frac < q90)).sum()); d = int(((~f) & (g.asn_frac < q90)).sum())
        n = a + b_ + c + d
        if n == 0 or (a + c) == 0 or (b_ + d) == 0: continue
        T.append([a, b_, c, d])
        num += a * d / n; den += b_ * c / n
    mh_or = num / den if den else np.nan
    # CMH chi-square (Mantel-Haenszel, continuity uncorrected)
    s_obs = s_exp = v = 0.0
    for a, b_, c, d in T:
        n = a + b_ + c + d
        n1 = a + c; m1 = a + b_
        e = n1 * m1 / n
        s_obs += a; s_exp += e
        v += n1 * (n - n1) * m1 * (n - m1) / (n * n * (n - 1)) if n > 1 else 0
    cmh_chi = (abs(s_obs - s_exp) - 0.5) ** 2 / v if v > 0 else np.nan
    cmh_p = float(st.chi2.sf(cmh_chi, 1)) if np.isfinite(cmh_chi) else np.nan
    sig = DS5[DS5.MW_p < 0.05] if len(DS5) else DS5
    rows.append({"set": name, "n": int(flag.sum()), "n_strata": len(DS5),
                 "MH_OR_asnrich": float(mh_or), "CMH_chi2": float(cmh_chi) if np.isfinite(cmh_chi) else np.nan,
                 "CMH_p": cmh_p,
                 "n_strata_regBG_diff_positive": int((DS5.median_reg > DS5.median_bg).sum()) if len(DS5) else 0,
                 "n_strata_MWsig": len(sig),
                 "median_withinstratum_diff": float((DS5.median_reg - DS5.median_bg).median()) if len(DS5) else np.nan})
E5b = pd.DataFrame(rows)
E5b.to_csv(f"{OUT}/M4RD_D5_stratified_summary.tsv", sep="\t", index=False)
print(E5b.to_string(), flush=True)

# ---------- D4: Zhang Pf constraint (Table S5 MIS/MFS) + Elsworth Pk categories ----------
import openpyxl
wz = openpyxl.load_workbook(f"{ROOT}/data/manual_inbox/M4R_ESSENTIALITY/zhang_S360_2018/NIHMS1004827-supplement-Table_S5.xlsx",
                            read_only=True, data_only=True)
ws = wz["Table S5"]; zrows = list(ws.iter_rows(values_only=True))
zhdr = [str(c) for c in zrows[1]]
print("zhang hdr:", zhdr, flush=True)
ZG = pd.DataFrame(zrows[2:], columns=zhdr)
ZG["Gene_ID"] = ZG["Gene_ID"].astype(str).str.strip()
ZG["MIS"] = pd.to_numeric(ZG["MIS"], errors="coerce")
ZG["MFS"] = pd.to_numeric(ZG["MFS"], errors="coerce")
ZG["n_ins"] = pd.to_numeric(ZG["# of insertions in CDS"], errors="coerce")
z = ZG.set_index("Gene_ID")
d4rows = [{"scope": "genome", "n": len(z), "frac_no_insertion": float((z.n_ins.fillna(0) == 0).mean()),
           "median_MIS": float(z.MIS.median()), "median_MFS": float(z.MFS.median())}]
for name, col in [("AP2", "is_AP2"), ("PUF_RNA", "is_PUF_RNA"), ("CHROM", "is_CHROM")]:
    genes = PF[PF[col]].gene.tolist()
    hit = [g for g in genes if g in z.index]
    if not hit: continue
    sub = z.loc[hit]
    bg = z.drop(index=hit, errors="ignore")
    u1, p1 = st.mannwhitneyu(sub.MIS.dropna(), bg.MIS.dropna(), alternative="two-sided")
    u2, p2 = st.mannwhitneyu(sub.MFS.dropna(), bg.MFS.dropna(), alternative="two-sided")
    d4rows.append({"scope": name, "n": len(hit),
                   "frac_no_insertion": float((sub.n_ins.fillna(0) == 0).mean()),
                   "median_MIS": float(sub.MIS.median()), "median_MFS": float(sub.MFS.median()),
                   "MW_MIS_p_vs_bg": float(p1), "MW_MFS_p_vs_bg": float(p2),
                   "bg_median_MIS": float(bg.MIS.median()), "bg_median_MFS": float(bg.MFS.median()),
                   "bg_frac_no_insertion": float((bg.n_ins.fillna(0) == 0).mean())})
D4 = pd.DataFrame(d4rows)
D4.to_csv(f"{OUT}/M4RD_D4_zhang_constraint.tsv", sep="\t", index=False)
we = openpyxl.load_workbook(f"{ROOT}/data/manual_inbox/M4R_ESSENTIALITY/elsworth_S387_2025/NIHMS2082619-supplement-Data_S3.xlsx",
                            read_only=True, data_only=True)
ws = we["S3a-Essentiality&Fitness.Score"]; erows = list(ws.iter_rows(values_only=True))
ehdr = [str(c) for c in erows[0]]
E = pd.DataFrame(erows[1:], columns=ehdr)
cat = {}
for c in ["Category1", "Category2", "Category3", "Category4"]:
    if c in E.columns: cat[c] = Counter(E[c].astype(str))
open(f"{OUT}/M4RD_D4_essentiality_note.txt", "w").write(
    f"Zhang S5 Pf genes={len(z)} (MIS/MFS/insertions; see M4RD_D4_zhang_constraint.tsv).\n"
    f"Elsworth S3a Pk genes={len(E)}; categories={dict(cat)}\n"
    f"NOTE: Elsworth IDs are Pk (PKNH_); Pf-regulator constraint uses Zhang Pf directly. "
    f"Pk-Pf ortholog-mapped constraint deferred to round 2.\n")
print(D4.to_string(), flush=True)
print("elsworth cats:", {k: dict(v) for k, v in cat.items()}, flush=True)

# refresh manifest
files = sorted([f for f in os.listdir(OUT) if f.startswith("M4RD_")])
mani = pd.DataFrame({"file": files, "sha256": [sha(f"{OUT}/{f}") for f in files]})
mani.to_csv(f"{OUT}/input_manifest.tsv", sep="\t", index=False)
print("WROTE2", files, flush=True)
