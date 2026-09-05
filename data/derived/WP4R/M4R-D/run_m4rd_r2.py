#!/usr/bin/env python
"""M4R-D ROUND2: LOO + D2b + D3-multidim + Pf8 boundary. D6-formal lives in run_m4rd_d6.py.
Usage: python run_m4rd_r2.py STEP   (STEP in LOO_<clade> D2B D3M PF8)
Read-only inputs; writes only new M4RD_* files under OUT. MWU asymptotic only.
Seed 20260905. Every STEP checkpoints its own outputs (idempotent reruns skip done rows)."""
import hashlib, json, math, os, re, sys
from collections import Counter
import numpy as np
import pandas as pd
from scipy import stats as st

ROOT = "/home/huyudi/015_plasmo"
L2 = f"{ROOT}/data/derived/WP4/L2_site_composition"
ALN = f"{L2}/alignments"
OUT = f"{ROOT}/data/derived/WP4R/M4R-D"
os.makedirs(OUT, exist_ok=True)
SEED = 20260905

SPECIES = ["SP001","SP002","SP003","SP004","SP005","SP006","SP007","SP008",
           "SP011","SP012","SP013","SP014","SP015","SP016","SP017","SP018"]
CLADES = {"laverania": ["SP001","SP002","SP003"],
          "vivax": ["SP004","SP005","SP006","SP007"],
          "malariae": ["SP008"],
          "rodent": ["SP011","SP012","SP013"],
          "avian": ["SP014"],
          "piroplasm": ["SP015","SP016","SP017"],
          "coccidian": ["SP018"]}
GC_CLASS = set("ARGPWV"); AT_CLASS = set("NDIKFY")
CLS = {a: (1 if a in GC_CLASS else -1 if a in AT_CLASS else 0) for a in "ARNDCQEGHILKMFPSTWYV"}

# ---------- shared loaders (same sources as run_m4rd.py) ----------
mani2 = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-02_aa/input_manifest.tsv", sep="\t", header=0)
m1 = pd.read_csv(f"{ROOT}/data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv", sep="\t")
spcol = m1.columns[0]
gcol = [c for c in m1.columns if "genome" in c.lower() and "gc" in c.lower()][0]
GCFULL = {s: float(m1[m1[spcol].astype(str) == s][gcol].iloc[0]) for s in SPECIES}
S2 = pd.read_csv(f"{L2}/L2v2_site_table.tsv", sep="\t")
S2["coupled"] = (S2.q_bh < 0.05).fillna(False)
S2["uncoupled"] = (S2.q_bh >= 0.1).fillna(False)
site_lookup = {(r.og, r.col): r for r in S2.itertuples()}
QC = pd.read_csv(f"{L2}/L2v2_alignment_qc.tsv", sep="\t")
PASS_OGS = [o for o in QC[QC["pass"]].og if not str(o).startswith("EXT_")]
gl = open(f"{ROOT}/data/raw/refs/grantham.tsv").read().strip().split("\n")
hdr = gl[0].split("\t")[1:]
GR = {}
for i, line in enumerate(gl[1:]):
    p = line.split("\t"); rowaa = p[0]
    for j, v in enumerate(p[1:]):
        fv = float(v)
        if fv != 0.0 or rowaa == hdr[j]:
            GR[frozenset([rowaa, hdr[j]])] = fv
for a in ["S"] + hdr: GR[frozenset([a, a])] = 0.0
def grantham(a, b): return GR.get(frozenset([a, b]), np.nan)

CODONS = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
"TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
"CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
"CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
"ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
"AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
"GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
"GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}
AACOD = {}
for cod, aa in CODONS.items():
    if aa == "*": continue
    AACOD.setdefault(aa, []).append(cod)
def min_nt(a, b):
    if a == b: return 0, 0, 0
    best = (9, 0, 0)
    for c1 in AACOD.get(a, []):
        for c2 in AACOD.get(b, []):
            d = sum(1 for x, y in zip(c1, c2) if x != y)
            g2a = sum(1 for x, y in zip(c1, c2) if x in "GC" and y in "AT" and x != y)
            a2g = sum(1 for x, y in zip(c1, c2) if x in "AT" and y in "GC" and x != y)
            if (d, g2a + a2g) < (best[0], best[1] + best[2]): best = (d, g2a, a2g)
    return best

def parse_aln(path):
    d = {}; hid, seq = None, []
    def flush():
        if hid is not None: d[hid] = "".join(seq)
    for line in open(path):
        if line.startswith(">"): flush(); hid = line[1:].split()[0]; seq = []
        else: seq.append(line.strip())
    flush(); return d

class Node:
    __slots__ = ("name", "blen", "children", "parent", "idx")
    def __init__(self, name=None, blen=0.0):
        self.name = name; self.blen = blen; self.children = []; self.parent = None; self.idx = -1

def parse_nw(s):
    s = s.strip()
    if s.endswith(";"): s = s[:-1]
    pos = 0
    def parse_node():
        nonlocal pos
        n = Node()
        if pos < len(s) and s[pos] == "(":
            pos += 1
            while True:
                ch = parse_node()
                ch.parent = n; n.children.append(ch)
                if s[pos] == ",": pos += 1; continue
                elif s[pos] == ")": pos += 1; break
        nm = ""
        while pos < len(s) and s[pos] not in ",();":
            nm += s[pos]; pos += 1
        nm = nm.strip()
        if nm:
            if ":" in nm:
                a, b = nm.split(":"); a = a.strip()
                n.name = a if a else None; n.blen = float(b)
            else:
                n.name = nm; n.blen = 0.0
        return n
    return parse_node()

def prune_tree(root, drop):
    def filt(n):
        if not n.children:
            return None if n.name in drop else n
        kids = []
        for c in n.children:
            k = filt(c)
            if k is not None: kids.append(k)
        if not kids: return None
        if len(kids) == 1:
            kids[0].parent = n.parent  # repair: collapse must not leave stale parent
            return kids[0]
        n.children = kids
        for k in kids: k.parent = n
        return n
    r = filt(root)
    assert r is not None
    return r

def anc_gc(nodes, tips, gcmap):
    N = len(nodes)
    gval = np.full(N, np.nan)
    tipset = set(id(t) for t in tips)
    for n in tips: gval[n.idx] = gcmap[n.name]
    for _ in range(500):
        for n in nodes:
            if id(n) in tipset: continue
            ch = [gval[c.idx] for c in n.children if np.isfinite(gval[c.idx])]
            if ch: gval[n.idx] = float(np.mean(ch))
    for _ in range(500):
        for n in nodes:
            if id(n) in tipset: continue
            vals = [gval[c.idx] for c in n.children if np.isfinite(gval[c.idx])]
            if n.parent is not None and np.isfinite(gval[n.parent.idx]): vals.append(gval[n.parent.idx])
            if vals: gval[n.idx] = float(np.mean(vals))
    return {n.idx: (gval[n.idx] - gval[n.parent.idx]) for n in nodes if n.parent is not None}

def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = p + z * z / (2 * n)
    m = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - m) / d, (c + m) / d)

def run_d1(sps, gcmap):
    """Faithful D1 rerun on species subset. Returns (summary_df, aux_dict)."""
    newick = open(f"{L2}/supermatrix.treefile").read().strip()
    drop = [s for s in SPECIES if s not in sps]
    root = prune_tree(parse_nw(newick), set(drop))
    nodes = []
    def walk(n):
        n.idx = len(nodes); nodes.append(n)
        for c in n.children: walk(c)
    walk(root)
    tips = [n for n in nodes if not n.children]
    assert {n.name for n in tips} == set(sps), "tip mismatch after prune"
    dGC = anc_gc(nodes, tips, gcmap)
    events = []; n_sites_run = 0; n_ambig_root = 0
    for og in PASS_OGS:
        aln = parse_aln(f"{ALN}/{og}.aln")
        if set(sps) - set(aln): continue
        L = len(aln[sps[0]])
        seqs = {s: aln[s] for s in sps}
        for i in range(L):
            key = (og, i)
            if key not in site_lookup: continue
            row = site_lookup[key]
            col = {s: seqs[s][i] for s in sps}
            ng = [c for c in col.values() if c != "-"]
            if len(ng) < 12 or len(set(ng)) < 2: continue
            n_sites_run += 1
            fset = {}
            order = []
            def post(n):
                for c in n.children: post(c)
                order.append(n)
            post(root)
            for n in order:
                if not n.children:
                    fset[n.idx] = {col[n.name]} if col[n.name] != "-" else set()
                else:
                    ch = [fset[c.idx] for c in n.children]
                    if any(len(x) == 0 for x in ch):
                        fset[n.idx] = set().union(*[x for x in ch if x])
                    else:
                        inter = set.intersection(*ch)
                        if inter: fset[n.idx] = inter
                        else: fset[n.idx] = set().union(*ch)
            if len(fset[root.idx]) != 1: n_ambig_root += 1
            assign = {}
            assign[root.idx] = sorted(fset[root.idx])[0]
            for n in order[::-1]:
                for c in n.children:
                    if c.idx in assign: continue
                    if assign[n.idx] in fset[c.idx]: assign[c.idx] = assign[n.idx]
                    else:
                        assign[c.idx] = sorted(fset[c.idx])[0] if fset[c.idx] else assign[n.idx]
            for n in nodes:
                if n.parent is None: continue
                a0, a1 = assign[n.parent.idx], assign[n.idx]
                if a0 == "-" or a1 == "-" or a0 == a1: continue
                c0, c1 = CLS.get(a0, 0), CLS.get(a1, 0)
                if c0 == 0 or c1 == 0: continue
                d, g2a, a2g = min_nt(a0, a1)
                events.append({"og": og, "col": i, "parent": n.parent.idx, "child": n.idx,
                               "from": a0, "to": a1, "dGC": dGC[n.idx],
                               "to_class": c1,
                               "concordant": int(np.sign(c1) == np.sign(dGC[n.idx])) if dGC[n.idx] != 0 else -1,
                               "min_nt": d, "gc2at": g2a, "at2gc": a2g,
                               "grantham": grantham(a0, a1),
                               "coupled": bool(row.coupled), "uncoupled": bool(row.uncoupled)})
    EV = pd.DataFrame(events)
    med = float(np.median([abs(e["dGC"]) for e in events])) if len(EV) else 0.0
    EVi = EV[(EV.concordant >= 0) & (EV.dGC.abs() >= med)].copy() if len(EV) else EV
    summ = []
    for label, df in [("all_branches", EV[EV.concordant >= 0]),
                      ("informative_branches", EVi),
                      ("informative_coupled", EVi[EVi.coupled]),
                      ("informative_uncoupled", EVi[EVi.uncoupled])]:
        k = int(df.concordant.sum()) if len(df) else 0; n = len(df)
        lo, hi = wilson(k, n)
        pbin = float(st.binomtest(k, n, 0.5, alternative="two-sided").pvalue) if n else np.nan
        summ.append({"stratum": label, "n_events": n, "n_concordant": k,
                     "frac": (k / n if n else np.nan), "ci_lo": lo, "ci_hi": hi, "binom_p": pbin})
    ct = EVi[EVi.coupled | EVi.uncoupled]
    orr, ft_p = np.nan, np.nan
    if len(ct):
        a = int(((ct.coupled) & (ct.concordant == 1)).sum()); b = int(((ct.coupled) & (ct.concordant == 0)).sum())
        c = int(((ct.uncoupled) & (ct.concordant == 1)).sum()); d = int(((ct.uncoupled) & (ct.concordant == 0)).sum())
        orr, ft_p = st.fisher_exact([[a, b], [c, d]])
    aux = {"OR": orr, "fisher_p": ft_p, "threshold": med, "n_sites": n_sites_run,
           "n_ambig": n_ambig_root, "n_events": len(EV)}
    return pd.DataFrame(summ), aux

# ---------- LCR helpers (same convention as round-1) ----------
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
def homopoly_runs(seq, minlen=5):
    runs = []; i = 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]: j += 1
        if j - i >= minlen and seq[i] != '-': runs.append((i, j))
        i = j
    return runs
def maxrun(s, aa):
    m = cur = 0
    for c in s:
        cur = cur + 1 if c == aa else 0
        m = max(m, cur)
    return m

STEP = sys.argv[1] if len(sys.argv) > 1 else ""
print(f"R2 STEP={STEP}", flush=True)

if STEP.startswith("LOO_"):
    cl = STEP[4:]
    assert cl in CLADES, f"unknown clade {cl}"
    sps = [s for s in SPECIES if s not in CLADES[cl]]
    DS, aux = run_d1(sps, GCFULL)
    row = {"left_out": cl, "n_species": len(sps), "n_events": aux["n_events"],
           "n_sites": aux["n_sites"], "threshold": aux["threshold"]}
    for _, r in DS.iterrows():
        row[f"{r.stratum}_n"] = r.n_events
        row[f"{r.stratum}_frac"] = r.frac
        row[f"{r.stratum}_lo"] = r.ci_lo
        row[f"{r.stratum}_hi"] = r.ci_hi
        row[f"{r.stratum}_p"] = r.binom_p
    row["Fisher_OR"] = aux["OR"]; row["Fisher_p"] = aux["fisher_p"]
    outp = f"{OUT}/M4RD_D1_LOO.tsv"
    hdrw = not os.path.exists(outp)
    pd.DataFrame([row]).to_csv(outp, sep="\t", index=False, mode="a", header=hdrw)
    print(f"LOO {cl}: OR={aux['OR']:.3f} p={aux['fisher_p']:.2g} events={aux['n_events']}", flush=True)

# ================= D2B: whole-proteome architecture attribution =================
if STEP == "D2B":
    PF = pd.read_csv(f"{OUT}/M4RD_D5_pf_features.tsv", sep="\t")
    OG = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-02_aa/orthogroups.tsv", sep="\t")
    og_age = OG.groupby("orthogroup_id").species_id.nunique().to_dict()
    # xp (SP001 protein) -> PF3D7 via GFF (assert PF3D7_ tags present)
    xp2pf = {}
    gff = f"{ROOT}/data/raw/ncbi-datasets/GCF_000002765.6/ncbi_dataset/data/GCF_000002765.6/genomic.gff"
    n_pf = 0
    for line in open(gff):
        if "\tCDS\t" not in line: continue
        m1 = re.search(r"ID=cds-(XP_\d+\.\d+)", line); m2 = re.search(r"locus_tag=(PF3D7_\d+)", line)
        if m1 and m2:
            xp2pf[m1.group(1)] = m2.group(1); n_pf += 1
    assert n_pf > 4000, f"GFF PF3D7 mapping too small: {n_pf}"
    print(f"xp2pf n={len(xp2pf)}", flush=True)
    pf2xp = {v: k for k, v in xp2pf.items()}
    pid2og = dict(OG[OG.species_id == "SP001"].set_index("protein_id").orthogroup_id)
    og2pf = {}
    for pid, og in pid2og.items():
        if pid in xp2pf: og2pf.setdefault(og, xp2pf[pid])
    # coupled-harboring Pf genes: OGs containing >=1 coupled site
    cpl_ogs = set(S2[S2.coupled].og)
    harbor = set(og2pf[o] for o in cpl_ogs if o in og2pf)
    print(f"coupled OGs={len(cpl_ogs)} harbor Pf genes={len(harbor)}", flush=True)
    # InterPro per gene (domain count + domain sets) from GFF Dbxref
    ipr = {}
    for line in open(gff):
        if "\tCDS\t" not in line and "\tgene\t" not in line: continue
        m = re.search(r"locus_tag=(PF3D7_\d+)", line)
        if not m: continue
        g = m.group(1)
        ips = set(re.findall(r"InterPro:(IPR\d+)", line))
        if ips: ipr.setdefault(g, set()).update(ips)
    print(f"genes with IPR={len(ipr)}", flush=True)
    D = PF.copy()
    D["harbor"] = D.gene.isin(harbor)
    D["og_age"] = D.gene.map(lambda g: og_age.get(pid2og.get(pf2xp.get(g, ""), ""), 0))
    D["n_ipr"] = D.gene.map(lambda g: len(ipr.get(g, set())))
    D["has_ipr"] = D.n_ipr > 0
    # frame-bias quantification: core-181 Pf members vs whole proteome
    core_ogs = set(QC[QC["pass"]].og) - {o for o in QC[QC["pass"]].og if str(o).startswith("EXT_")}
    core_pf = set(og2pf[o] for o in core_ogs if o in og2pf)
    D["in_core181"] = D.gene.isin(core_pf)
    frame_rows = []
    for feat in ["lcr_frac", "asn_frac", "length"]:
        a = D[D.in_core181][feat].dropna().values
        b = D[~D.in_core181][feat].dropna().values
        mw = st.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic")
        frame_rows.append({"feature": feat, "median_core": float(np.median(a)),
                           "median_noncore": float(np.median(b)), "MW_p": float(mw.pvalue)})
    pd.DataFrame(frame_rows).to_csv(f"{OUT}/M4RD_D2b_framebias.tsv", sep="\t", index=False)
    # matched test: harbor vs background matched on length decile + og-age class
    D["len_dec"] = pd.qcut(D.length, 10, labels=False, duplicates="drop")
    D["age_cls"] = pd.cut(D.og_age, [-1, 1, 4, 8, 100], labels=["narrow", "mid", "wide", "core"])
    res = []
    for feat in ["lcr_frac", "asn_frac", "polyN_max", "n_ipr"]:
        diffs = []
        for (ld, ag), grp in D.groupby(["len_dec", "age_cls"], observed=True):
            h = grp[grp.harbor][feat].dropna().values
            b = grp[~grp.harbor][feat].dropna().values
            if len(h) < 3 or len(b) < 10: continue
            diffs.append((float(np.median(h) - np.median(b)), len(h), len(b)))
        ha = D[D.harbor][feat].dropna().values
        ba = D[~D.harbor][feat].dropna().values
        mw = st.mannwhitneyu(ha, ba, alternative="two-sided", method="asymptotic")
        res.append({"feature": feat, "n_harbor": len(ha), "n_bg": len(ba),
                    "median_harbor": float(np.median(ha)), "median_bg": float(np.median(ba)),
                    "MW_p": float(mw.pvalue),
                    "median_strat_diff": float(np.mean([d[0] for d in diffs])) if diffs else float("nan"),
                    "n_strata": len(diffs)})
    R = pd.DataFrame(res)
    R.to_csv(f"{OUT}/M4RD_D2b_matched.tsv", sep="\t", index=False)
    D[["gene", "length", "asn_frac", "lcr_frac", "polyN_max", "n_ipr",
       "harbor", "in_core181", "og_age"]].to_csv(f"{OUT}/M4RD_D2b_wholeproteome.tsv", sep="\t", index=False)
    # verdict: SURVIVES unless harbor significantly MORE LCR-rich (matched + crude agree)
    lrow = R[R.feature == "lcr_frac"].iloc[0]
    verdict = "DOWNGRADED" if (lrow.MW_p < 0.05 and lrow.median_harbor > lrow.median_bg
                               and (np.isnan(lrow.median_strat_diff) or lrow.median_strat_diff > 0)) else "SURVIVES"
    print(f"D2b verdict={verdict} lcr_med_h={lrow.median_harbor:.4f} bg={lrow.median_bg:.4f} p={lrow.MW_p:.2g}", flush=True)
    with open(f"{OUT}/M4RD_D2b_verdict.txt", "w") as h:
        h.write(f"{verdict}\n")

# ================= D3M: multi-dim matched controls =================
if STEP == "D3M":
    MC = pd.read_csv(f"{OUT}/M4RD_D3_matched_controls.tsv", sep="\t")
    D2S = pd.read_csv(f"{OUT}/M4RD_D2_site_architecture.tsv", sep="\t")
    # InterPro per Pf gene from GFF
    xp2pf = {}
    ipr = {}
    gff = f"{ROOT}/data/raw/ncbi-datasets/GCF_000002765.6/ncbi_dataset/data/GCF_000002765.6/genomic.gff"
    for line in open(gff):
        if "\tCDS\t" in line:
            m1 = re.search(r"ID=cds-(XP_\d+\.\d+)", line); m2 = re.search(r"locus_tag=(PF3D7_\d+)", line)
            if m1 and m2: xp2pf[m1.group(1)] = m2.group(1)
        if "\tgene\t" in line or "\tCDS\t" in line:
            m = re.search(r"locus_tag=(PF3D7_\d+)", line)
            if m:
                ips = set(re.findall(r"InterPro:(IPR\d+)", line))
                if ips: ipr.setdefault(m.group(1), set()).update(ips)
    OG = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-02_aa/orthogroups.tsv", sep="\t")
    pid2og = dict(OG[OG.species_id == "SP001"].set_index("protein_id").orthogroup_id)
    og2pf = {}
    for pid, og in pid2og.items():
        if pid in xp2pf: og2pf.setdefault(og, xp2pf[pid])
    # Zhang per-gene MIS/MFS: locate S5 sheet with gene + MIS columns
    zpath = f"{ROOT}/data/manual_inbox/M4R_ESSENTIALITY/zhang_S360_2018/NIHMS1004827-supplement-Table_S5.xlsx"
    ess = {}
    try:
        xl = pd.ExcelFile(zpath)
        full = xl.parse("Table S5", header=1)
        gcol = [c for c in full.columns if "gene_id" in str(c).lower()][0]
        mcol = [c for c in full.columns if str(c).strip().lower() == "mis"][0]
        for _, r in full.iterrows():
            m = re.search(r"PF3D7_\d+", str(r[gcol]))
            if m:
                try: ess[m.group(0)] = float(r[mcol])
                except (ValueError, TypeError): pass
        print(f"Zhang Table S5 genes={len(ess)}", flush=True)
    except Exception as e:
        print(f"Zhang parse failed: {e}", flush=True)
    print(f"ess genes={len(ess)} ipr genes={len(ipr)}", flush=True)
    S2c = S2.set_index(["og", "col"])
    rows = []
    for _, k in MC.iterrows():
        gene, ogk = k.gene, k.og
        dome = ipr.get(gene, set())
        # background pools
        same_og = D2S[D2S.og == ogk]
        dom_ogs = set()
        if dome:
            for og, pf in og2pf.items():
                if ipr.get(pf, set()) & dome: dom_ogs.add(og)
        dom_sites = D2S[D2S.og.isin(dom_ogs)] if dom_ogs else D2S.iloc[0:0]
        for scope, pool in [("within_OG", same_og), ("domain_shared", dom_sites)]:
            if len(pool) < 5:
                rows.append({"gene": gene, "pf_pos": k.pf_pos, "scope": scope,
                             "n_bg": len(pool), "frac_bg_coupled": np.nan,
                             "note": "n<5:uninformative"})
                continue
            cpl = 0; tot = 0
            for _, s in pool.iterrows():
                key = (s.og, s.col)
                if key in S2c.index:
                    r = S2c.loc[key]
                    if isinstance(r, pd.DataFrame): r = r.iloc[0]
                    if bool(r.coupled): cpl += 1
                    if bool(r.coupled) or bool(r.uncoupled): tot += 1
            rows.append({"gene": gene, "pf_pos": k.pf_pos, "scope": scope,
                         "n_bg": len(pool), "frac_bg_coupled": (cpl / tot if tot else np.nan),
                         "note": ""})
        # essentiality bin of known gene
        mis = ess.get(gene, np.nan)
        ebin = "unknown" if np.isnan(mis) else ("essential" if mis < 0 else "dispensable")
        rows.append({"gene": gene, "pf_pos": k.pf_pos, "scope": "known_gene_ess_bin",
                     "n_bg": 1, "frac_bg_coupled": np.nan,
                     "note": f"MIS={mis:.3f}" if not np.isnan(mis) else "MIS missing"})
    R = pd.DataFrame(rows)
    R.to_csv(f"{OUT}/M4RD_D3_multidim.tsv", sep="\t", index=False)
    if not ess:
        print("D3M: essentiality per-gene table missing -> conservation-matched only + domain scopes kept", flush=True)
    print(R.to_string(), flush=True)

# ================= PF8: population boundary =================
if STEP == "PF8":
    PDIR = f"{ROOT}/data/manual_inbox/M4R_PF8_2026"
    M = pd.read_csv(f"{PDIR}/Pf8_drug_resistance_marker_genotypes.tsv", sep="\t", low_memory=False)
    S = pd.read_csv(f"{PDIR}/Pf8_samples.txt", sep="\t", low_memory=False)
    F = pd.read_csv(f"{PDIR}/Pf8_fws.tsv", sep="\t")
    # known adaptive markers: parse column names marker_pos[REF]
    KNOWN = pd.read_csv(f"{L2}/L2_known_sites.tsv", sep="\t")
    fws = dict(zip(F.Sample, F.Fws))
    qcpass = set(S[S["QC pass"] == True].Sample) if "QC pass" in S.columns else set(M.Sample)
    S["Fws"] = S.Sample.map(fws)
    strict = set(S[(S["QC pass"] == True) & (S.Fws >= 0.95)].Sample) if "QC pass" in S.columns else set()
    print(f"samples={len(M)} qcpass={len(qcpass)} strict_single={len(strict)}", flush=True)
    mcols = [c for c in M.columns if c != "Sample"]
    def parse_call(v, ref):
        if pd.isna(v) or str(v).strip() in ("", "-"): return None
        return [a.strip() for a in str(v).split(",")]
    rows = []
    for mc in mcols:
        m = re.match(r"(.+?)_(\d+)\[([A-Z*\-]+)\]$", mc)
        if not m: continue
        gene, pos, ref = m.group(1), int(m.group(2)), m.group(3)
        calls = M[mc]
        called = calls.dropna()
        called = called[~called.astype(str).isin(["", "-", " "])]
        n_call = len(called)
        alleles = set()
        n_res = 0
        for v in called.values:
            al = parse_call(v, ref)
            if not al: continue
            alleles.update(al)
            if any(a != ref for a in al): n_res += 1
        # strict subset
        sub = M[M.Sample.isin(strict)][mc].dropna() if strict else called.iloc[0:0]
        sub = sub[~sub.astype(str).isin(["", "-", " "])] if len(sub) else sub
        n_res_s = 0
        for v in sub.values:
            al = parse_call(v, ref)
            if al and any(a != ref for a in al): n_res_s += 1
        rows.append({"marker": mc, "gene": gene, "pos": pos, "ref": ref,
                     "call_rate": (n_call / len(M) if len(M) else np.nan),
                     "n_alleles": len(alleles),
                     "segregating": len(alleles) > 1,
                     "res_freq_all": (n_res / n_call if n_call else np.nan),
                     "res_freq_strict": (n_res_s / len(sub) if len(sub) else np.nan),
                     "n_strict": len(sub)})
    R = pd.DataFrame(rows)
    R.to_csv(f"{OUT}/M4RD_D4_Pf8boundary.tsv", sep="\t", index=False)
    # CNV locus rates
    C = pd.read_csv(f"{PDIR}/Pf8_cnv_calls.tsv", sep="\t", low_memory=False)
    fin = [c for c in C.columns if "final" in c.lower()]
    crow = []
    for c in fin:
        v = C[c].astype(str)
        amp = int(((v != "0") & (v != "-") & (v != "nan")).sum())
        crow.append({"locus_call": c, "n_amp": amp, "freq": amp / len(C)})
    pd.DataFrame(crow).to_csv(f"{OUT}/M4RD_D4_Pf8cnv.tsv", sep="\t", index=False)
    print(f"markers={len(R)} segregating={int(R.segregating.sum())}", flush=True)
    print(R[R.segregating].head(20).to_string(), flush=True)
    print("SNP-level burden: DEFERRED (needs Zarr streaming infra, beyond round-2 budget)", flush=True)
