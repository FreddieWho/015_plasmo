#!/usr/bin/env python
"""M4R-D ROUND1 first-gain analysis: branch-aware L2 + architecture partition +
matched adaptive-site controls + regulatory enrichment counterfactual + natural experiment.
Read-only inputs; writes only to data/derived/WP4R/M4R-D/."""
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
SKIP_D1 = ("--skip-d1" in sys.argv)
rng = np.random.default_rng(SEED)
SPECIES = ["SP001","SP002","SP003","SP004","SP005","SP006","SP007","SP008",
           "SP011","SP012","SP013","SP014","SP015","SP016","SP017","SP018"]
GC_CLASS = set("ARGPWV"); AT_CLASS = set("NDIKFY")
CLS = {a: (1 if a in GC_CLASS else -1 if a in AT_CLASS else 0) for a in "ARNDCQEGHILKMFPSTWYV"}

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

# ---------- inputs ----------
mani2 = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-02_aa/input_manifest.tsv", sep="\t", header=0)
faa = {}
for _, r in mani2.iterrows():
    m = re.match(r"protein_(SP\d+)_", str(r.iloc[0]))
    if m: faa[m.group(1)] = f"{ROOT}/{r.iloc[1]}"
m1 = pd.read_csv(f"{ROOT}/data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv", sep="\t")
spcol = m1.columns[0]
gcol = [c for c in m1.columns if "genome" in c.lower() and "gc" in c.lower()][0]
GC = {s: float(m1[m1[spcol].astype(str) == s][gcol].iloc[0]) for s in SPECIES}
S2 = pd.read_csv(f"{L2}/L2v2_site_table.tsv", sep="\t")
QC = pd.read_csv(f"{L2}/L2v2_alignment_qc.tsv", sep="\t")
PASS_OGS = [o for o in QC[QC["pass"]].og if not str(o).startswith("EXT_")]
KNOWN = pd.read_csv(f"{L2}/L2_known_sites.tsv", sep="\t")
ADAPT = pd.read_csv(f"{L2}/L2v2_adaptive_substitutions_grantham.tsv", sep="\t")
gl = open(f"{ROOT}/data/raw/refs/grantham.tsv").read().strip().split("\n")
hdr = gl[0].split("\t")[1:]
GR = {}
aas20 = ["S"] + hdr
for i, line in enumerate(gl[1:]):
    p = line.split("\t"); rowaa = p[0]
    for j, v in enumerate(p[1:]):
        fv = float(v)
        if fv != 0.0 or rowaa == hdr[j]:
            GR[frozenset([rowaa, hdr[j]])] = fv
for a in aas20: GR[frozenset([a, a])] = 0.0
def grantham(a, b): return GR.get(frozenset([a, b]), np.nan)

def parse_aln(path):
    d = {}; hid, seq = None, []
    def flush():
        if hid is not None: d[hid] = "".join(seq)
    for line in open(path):
        if line.startswith(">"): flush(); hid = line[1:].split()[0]; seq = []
        else: seq.append(line.strip())
    flush(); return d

def parse_faa(path):
    return parse_aln(path)

PROT = {s: parse_faa(faa[s]) for s in SPECIES}
OG = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-02_aa/orthogroups.tsv", sep="\t")
OG_SP1 = dict(OG[OG.species_id == "SP001"].set_index("orthogroup_id").protein_id)

# ---------- tree parse (tiny newick) ----------
newick = open(f"{L2}/supermatrix.treefile").read().strip()
class Node:
    __slots__ = ("name", "blen", "children", "parent", "idx")
    def __init__(self, name=None, blen=0.0):
        self.name = name; self.blen = blen; self.children = []; self.parent = None; self.idx = -1
def parse_newick(s):
    s = s.strip().rstrip(";")
    stack = []; cur = None; i = 0; root = None
    token = ""
    def finish_token(node):
        nonlocal token
        token = token.strip()
        if token:
            if ":" in token:
                nm, bl = token.split(":"); nm = nm.strip()
                node.name = nm if nm else None; node.blen = float(bl)
            else:
                node.name = token; node.blen = 0.0
        token = ""
    while i < len(s):
        c = s[i]
        if c == "(":
            n = Node(); n.parent = cur
            if cur is not None: cur.children.append(n)
            else: root = n
            stack.append(n); cur = n; i += 1
        elif c == ",":
            finish_token(Node())
            n = stack[-1]
            # token belonged to a tip child of current internal
            tip = Node(); tip.parent = cur
            nm_bl = token.strip()
            token = ""
            # re-parse: we consumed token already; reconstruct
            cur.children.append(tip)
            # fix: token was already finished; set from saved
            tip.name, tip.blen = None, 0.0
            # redo properly below
            cur.children.pop()
            i += 1
        elif c == ")":
            # last child before ) is a tip described by token
            if token.strip():
                tip = Node(); tip.parent = cur; cur.children.append(tip)
                finish_token(tip)
            i += 1
            # optional internal name/blen follows
            j = i
            tok = ""
            while j < len(s) and s[j] not in ",();":
                tok += s[j]; j += 1
            if tok.strip():
                token = tok; finish_token(cur); token = ""
            i = j
            stack.pop()
            cur = stack[-1] if stack else None
        else:
            token += c; i += 1
    # The comma branch above is broken; use simpler robust parser instead
    return None

# Simpler robust recursive-descent parser
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
    root = parse_node()
    return root
TROOT = parse_nw(newick)
NODES = []
def walk(n):
    n.idx = len(NODES); NODES.append(n)
    for c in n.children: walk(c)
walk(TROOT)
TIPS = [n for n in NODES if not n.children]
EDGES = [(n.parent.idx, n.idx) for n in NODES if n.parent is not None]
TIPSEQ = {n.name: n.idx for n in TIPS}
assert set(TIPSEQ) == set(SPECIES), f"tip mismatch {set(TIPSEQ)^set(SPECIES)}"
N = len(NODES)

# ---------- ancestral GC (squared-change parsimony, iterative mean) ----------
gval = np.full(N, np.nan)
for n in TIPS: gval[n.idx] = GC[n.name]
for _ in range(500):
    for n in NODES:
        if n in TIPS: continue
        ch = [gval[c.idx] for c in n.children if np.isfinite(gval[c.idx])]
        if ch: gval[n.idx] = float(np.mean(ch))
# one up-pass refine (parent info)
for _ in range(500):
    for n in NODES:
        if n in TIPS: continue
        vals = [gval[c.idx] for c in n.children if np.isfinite(gval[c.idx])]
        if n.parent is not None and np.isfinite(gval[n.parent.idx]): vals.append(gval[n.parent.idx])
        if vals: gval[n.idx] = float(np.mean(vals))
dGC = {n.idx: (gval[n.idx] - gval[n.parent.idx]) for n in NODES if n.parent is not None}

# ---------- genetic code accessibility ----------
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
    if a == b: return 0, 0, 0  # dist, gc_to_at, at_to_gc
    best = (9, 0, 0)
    for c1 in AACOD.get(a, []):
        for c2 in AACOD.get(b, []):
            d = sum(1 for x, y in zip(c1, c2) if x != y)
            g2a = sum(1 for x, y in zip(c1, c2) if x in "GC" and y in "AT" and x != y)
            a2g = sum(1 for x, y in zip(c1, c2) if x in "AT" and y in "GC" and x != y)
            if (d, g2a + a2g) < (best[0], best[1] + best[2]): best = (d, g2a, a2g)
    return best

# ---------- D1: Fitch parsimony per site, substitutions -> branches ----------
S2["coupled"] = (S2.q_bh < 0.05).fillna(False)
S2["uncoupled"] = (S2.q_bh >= 0.1).fillna(False)
site_lookup = {(r.og, r.col): r for r in S2.itertuples()}
events = []  # one row per branch substitution
n_sites_run = 0; n_ambig_root = 0; n_singleton_excluded = 0
for og in PASS_OGS:
    aln = parse_aln(f"{ALN}/{og}.aln")
    if set(SPECIES) - set(aln): continue
    L = len(aln[SPECIES[0]])
    seqs = {s: aln[s] for s in SPECIES}
    for i in range(L):
        key = (og, i)
        if key not in site_lookup: continue  # invariant / <12 ng / filtered
        row = site_lookup[key]
        col = {s: seqs[s][i] for s in SPECIES}
        ng = [c for c in col.values() if c != "-"]
        if len(ng) < 12 or len(set(ng)) < 2: continue
        n_sites_run += 1
        # Fitch post-order
        fset = {}
        order = []
        def post(n):
            for c in n.children: post(c)
            order.append(n)
        post(TROOT)
        steps = 0
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
                    else: fset[n.idx] = set().union(*ch); steps += 1
        if len(fset[TROOT.idx]) != 1: n_ambig_root += 1
        # pre-order assignment (parent-biased = DELTRAN-ish)
        assign = {}
        rootset = sorted(fset[TROOT.idx])
        assign[TROOT.idx] = rootset[0]
        for n in order[::-1]:
            for c in n.children:
                if c.idx in assign: continue
                if assign[n.idx] in fset[c.idx]: assign[c.idx] = assign[n.idx]
                else:
                    assign[c.idx] = sorted(fset[c.idx])[0] if fset[c.idx] else assign[n.idx]
        for n in NODES:
            if n.parent is None: continue
            a0, a1 = assign[n.parent.idx], assign[n.idx]
            if a0 == "-" or a1 == "-" or a0 == a1: continue
            c0, c1 = CLS.get(a0, 0), CLS.get(a1, 0)
            if c0 == 0 or c1 == 0: continue
            d, g2a, a2g = min_nt(a0, a1)
            events.append({"og": og, "col": i, "parent": n.parent.idx, "child": n.idx,
                           "from": a0, "to": a1, "dGC": dGC[n.idx],
                           "to_class": c1, "concordant": int(np.sign(c1) == np.sign(dGC[n.idx])) if dGC[n.idx] != 0 else -1,
                           "min_nt": d, "gc2at": g2a, "at2gc": a2g,
                           "grantham": grantham(a0, a1),
                           "coupled": bool(row.coupled), "uncoupled": bool(row.uncoupled),
                           "r_gc": row.r_gc})
EV = pd.DataFrame(events)
EV.to_csv(f"{OUT}/M4RD_D1_branch_substitutions.tsv", sep="\t", index=False)

def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = p + z * z / (2 * n)
    m = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - m) / d, (c + m) / d)

med_abs_dgc = float(np.median([abs(e["dGC"]) for e in events])) if len(EV) else 0.0
EVi = EV[(EV.concordant >= 0) & (EV.dGC.abs() >= med_abs_dgc)].copy() if len(EV) else EV
summ = []
for label, df in [("all_branches", EV[EV.concordant >= 0]),
                  ("informative_branches", EVi),
                  ("informative_coupled", EVi[EVi.coupled]),
                  ("informative_uncoupled", EVi[EVi.uncoupled]),
                  ("informative_1nt", EVi[EVi.min_nt == 1]),
                  ("informative_gen2nt", EVi[EVi.min_nt >= 2])]:
    k = int(df.concordant.sum()) if len(df) else 0; n = len(df)
    lo, hi = wilson(k, n)
    pbin = float(st.binomtest(k, n, 0.5, alternative="two-sided").pvalue) if n else np.nan
    summ.append({"stratum": label, "n_events": n, "n_concordant": k,
                 "frac": (k / n if n else np.nan), "ci_lo": lo, "ci_hi": hi, "binom_p": pbin})
DS = pd.DataFrame(summ)
# coupled vs uncoupled concordance (Fisher)
ct = EVi[EVi.coupled | EVi.uncoupled]
ft_p = np.nan; orr = np.nan
if len(ct):
    a = int(((ct.coupled) & (ct.concordant == 1)).sum()); b = int(((ct.coupled) & (ct.concordant == 0)).sum())
    c = int(((ct.uncoupled) & (ct.concordant == 1)).sum()); d = int(((ct.uncoupled) & (ct.concordant == 0)).sum())
    orr, ft_p = st.fisher_exact([[a, b], [c, d]])
DS2 = pd.DataFrame([{"comparison": "coupled_vs_uncoupled_concordance", "OR": orr, "fisher_p": ft_p,
                     "median_abs_dGC_threshold": med_abs_dgc, "n_sites_parsimony": n_sites_run,
                     "n_ambig_root": n_ambig_root}])
DS.to_csv(f"{OUT}/M4RD_D1_concordance_summary.tsv", sep="\t", index=False)
DS2.to_csv(f"{OUT}/M4RD_D1_concordance_aux.tsv", sep="\t", index=False)
print("D1 events:", len(EV), "sites:", n_sites_run, flush=True)
print(DS.to_string(), flush=True)

# ---------- D2: architecture partition (LCR via window64 entropy<1.5) ----------
# Project-convention LCR: exact M2 logic (window64/step8/entropy<1.5 with
# distinct>=8 -> 2.0 shortcut). Sensitivity: thr=2.5 (disclosed non-standard).
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

def lcr_mask_sens(seq, thr=2.5):
    return lcr_mask(seq, thr=thr)

def homopoly_runs(seq, minlen=5):
    runs = []; i = 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]: j += 1
        if j - i >= minlen and seq[i] != "-": runs.append((i, j, seq[i]))
        i = j
    return runs

site_arch = []
for og in PASS_OGS:
    aln = parse_aln(f"{ALN}/{og}.aln")
    if set(SPECIES) - set(aln): continue
    pf = aln["SP001"]
    mask = lcr_mask(pf)
    masks = lcr_mask_sens(pf)
    mask2 = lcr_mask_sens(pf, 2.0)
    runs = homopoly_runs(pf)
    runpos = {}
    for (a, b, aa) in runs:
        for k in range(a, b): runpos[k] = aa
    L = len(pf)
    for i in range(L):
        key = (og, i)
        if key not in site_lookup: continue
        row = site_lookup[key]
        colaa = [aln[s][i] for s in SPECIES if aln[s][i] != "-"]
        if len(colaa) < 12: continue
        cons = max(Counter(colaa).values()) / len(colaa)
        site_arch.append({"og": og, "col": i, "coupled": bool(row.coupled),
                          "uncoupled": bool(row.uncoupled), "r_gc": row.r_gc,
                          "mean_grantham": row.mean_grantham, "n_class": row.n_class,
                          "in_LCR": bool(mask[i]), "in_LCR_sens": bool(masks[i]),
                          "in_LCR_sens2": bool(mask2[i]),
                          "in_homorun": i in runpos,
                          "run_aa": runpos.get(i, ""),
                          "polyN": bool(runpos.get(i, "") == "N"),
                          "conservation": cons, "n_aa": len(set(colaa))})
AR = pd.DataFrame(site_arch)
AR.to_csv(f"{OUT}/M4RD_D2_site_architecture.tsv", sep="\t", index=False)
# overall + stratified stats
comp = AR[AR.coupled | AR.uncoupled].copy()
def grp(df, name, col="in_LCR"):
    cc = df[df.coupled]; uu = df[df.uncoupled]
    a = int(cc[col].sum()); b = int((~cc[col]).sum())
    c = int(uu[col].sum()); d = int((~uu[col]).sum())
    if min(a, b, c, d) < 0 or (a + b) == 0 or (c + d) == 0:
        orr, p = np.nan, np.nan
    else:
        try: orr, p = st.fisher_exact([[a, b], [c, d]])
        except Exception: orr, p = np.nan, np.nan
    return {"stratum": name, "n_coupled": len(cc), "n_uncoupled": len(uu),
            "frac_LCR_coupled": float(cc[col].mean()), "frac_LCR_uncoupled": float(uu[col].mean()),
            "OR_LCR": float(orr) if np.isfinite(orr) else np.nan, "fisher_p": float(p) if p == p else np.nan,
            "frac_polyN_coupled": float(cc.polyN.mean()), "frac_polyN_uncoupled": float(uu.polyN.mean()),
            "median_grantham_coupled": float(cc.mean_grantham.median()),
            "median_grantham_uncoupled": float(uu.mean_grantham.median())}
rows2 = [grp(comp, "overall_LCRconv"), grp(comp, "sens_LCRthr2.5", "in_LCR_sens"),
         grp(comp, "sens_LCRthr2.0", "in_LCR_sens2"),
         grp(comp, "homorun_ge5", "in_homorun")]
# within-OG permutation: is the coupled-LCR excess beyond OG structure?
obs = (comp[comp.coupled].in_LCR.mean() - comp[comp.uncoupled].in_LCR.mean())
perm = []
sub = comp[["og", "coupled", "in_LCR"]].copy()
for b in range(2000):
    sh = sub.copy()
    sh["coupled"] = sh.groupby("og")["coupled"].transform(lambda x: rng.permutation(x.values))
    perm.append(sh[sh.coupled].in_LCR.mean() - sh[sh.coupled == False].in_LCR.mean())
perm = np.array(perm)
p_perm = (1 + int((perm >= obs).sum())) / (1 + len(perm))
rows2.append({"stratum": "withinOG_permutation_p", "n_coupled": int(comp.coupled.sum()),
              "n_uncoupled": int(comp.uncoupled.sum()), "frac_LCR_coupled": float(comp[comp.coupled].in_LCR.mean()),
              "frac_LCR_uncoupled": float(comp[comp.uncoupled].in_LCR.mean()), "OR_LCR": np.nan,
              "fisher_p": float(p_perm), "frac_polyN_coupled": np.nan, "frac_polyN_uncoupled": np.nan,
              "median_grantham_coupled": np.nan, "median_grantham_uncoupled": np.nan,
              "note": f"obs_diff={obs:.4f}"})
# logistic: coupled ~ LCR + conservation + n_class (all sites with defined coupling)
from sklearn.linear_model import LogisticRegression
lz = comp.dropna(subset=["conservation", "n_class"]).copy()
X = np.column_stack([lz.in_LCR.astype(float), lz.conservation, lz.n_class.astype(float)])
X = (X - X.mean(0)) / np.where(X.std(0) == 0, 1, X.std(0))
lr = LogisticRegression(max_iter=2000).fit(X, lz.coupled.astype(int))
rows2.append({"stratum": "logistic_coef", "n_coupled": int(lz.coupled.sum()),
              "n_uncoupled": int(lz.uncoupled.sum()), "frac_LCR_coupled": float(lr.coef_[0][0]),
              "frac_LCR_uncoupled": float(lr.coef_[0][1]), "OR_LCR": float(lr.coef_[0][2]),
              "fisher_p": float(lr.score(X, lz.coupled.astype(int))),
              "frac_polyN_coupled": np.nan, "frac_polyN_uncoupled": np.nan,
              "median_grantham_coupled": np.nan, "median_grantham_uncoupled": np.nan,
              "note": "coef order: LCR, conservation, n_class; fisher_p col = accuracy"})
D2S = pd.DataFrame(rows2)
D2S.to_csv(f"{OUT}/M4RD_D2_partition_summary.tsv", sep="\t", index=False)
print(D2S.to_string(), flush=True)

# ---------- D3: matched adaptive-site controls ----------
# conservation per known site from alignments
mrows = []
for r in KNOWN.itertuples():
    if r.status != "mapped": continue
    alnpath = f"{ALN}/{r.og}.aln"
    if not os.path.exists(alnpath):
        # EXT drug-gene OGs live under EXT_ prefix
        cand = [f for f in os.listdir(ALN) if f.endswith(".aln") and r.og in f]
        alnpath = f"{ALN}/{cand[0]}" if cand else None
    if alnpath is None:
        mrows.append({"gene": r.gene, "drug": r.drug, "pf_pos": r.pf_pos, "og": r.og,
                      "cons_known": np.nan, "n_bg": 0, "note": "alignment_missing"})
        continue
    aln = parse_aln(alnpath)
    spp = [s for s in SPECIES if s in aln]
    i = int(r.aln_col)
    colaa = [aln[s][i] for s in spp if aln[s][i] != "-"]
    cons = max(Counter(colaa).values()) / len(colaa) if colaa else np.nan
    # background: same OG v2 sites within cons +- 0.1
    bg = AR[(AR.og == r.og) & (AR.conservation >= cons - 0.1) & (AR.conservation <= cons + 0.1)]
    bg_scope = "sameOG"
    if len(bg) == 0:
        bg = AR[(AR.conservation >= cons - 0.1) & (AR.conservation <= cons + 0.1)]
        bg_scope = "global_cons_matched"
    defined = bg.r_gc.notna()
    mrows.append({"gene": r.gene, "drug": r.drug, "pf_pos": r.pf_pos, "og": r.og,
                  "cons_known": cons, "n_bg": len(bg), "bg_scope": bg_scope,
                  "in_v2_OG": bool(r.og in PASS_OGS),
                  "frac_bg_r_defined": float(defined.mean()) if len(bg) else np.nan,
                  "median_abs_r_bg": float(bg.r_gc.abs().median()) if len(bg) else np.nan,
                  "frac_bg_coupled": float(bg[bg.r_gc.notna()].coupled.mean()) if defined.sum() else 0.0,
                  "median_grantham_bg": float(bg.mean_grantham.median()) if len(bg) else np.nan,
                  "median_n_aa_bg": float(bg.n_aa.median()) if len(bg) else np.nan,
                  "known_coupled": False, "known_r_defined": False})
M3 = pd.DataFrame(mrows)
M3.to_csv(f"{OUT}/M4RD_D3_matched_controls.tsv", sep="\t", index=False)
# adaptive-sub Grantham quantile vs coupled-site Grantham distribution
cg = AR[AR.coupled].mean_grantham.values
aq = []
for r in ADAPT.itertuples():
    q = float((cg <= r.grantham).mean()) if len(cg) else np.nan
    aq.append({"gene": r.gene, "sub": r.sub, "grantham": r.grantham, "class": r[3],
               "quantile_in_coupled_dist": q})
AQ = pd.DataFrame(aq)
AQ.to_csv(f"{OUT}/M4RD_D3_adaptive_grantham_quantiles.tsv", sep="\t", index=False)
print(M3.to_string(), flush=True)
print("adaptive grantham median:", float(ADAPT.grantham.median()),
      "coupled bg median:", float(np.median(cg)) if len(cg) else None, flush=True)

# ---------- D5: Pf-wide Asn/LCR/IDR split + AP2/PUF enrichment + counterfactual ----------
pf_path = f"{ROOT}/data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7_AnnotatedProteins.fasta"
genes, seqs, prods = [], [], []
hid = None; cur = []
def flush2():
    if hid:
        m = re.search(r"gene=(PF3D7_\S+)", hid)
        pm = re.search(r"gene_product=([^|]+)", hid)
        genes.append(m.group(1) if m else hid.split()[0])
        prods.append(pm.group(1).strip().lower() if pm else "")
        seqs.append("".join(cur))
for line in open(pf_path):
    if line.startswith(">"): flush2(); hid = line[1:].rstrip("\n"); cur = []
    else: cur.append(line.strip())
flush2()
PF = pd.DataFrame({"gene": genes, "prod_desc": prods, "seq": seqs})
PF["length"] = PF.seq.str.len()
PF["asn_frac"] = PF.seq.apply(lambda s: s.count("N") / len(s) if len(s) else 0)
PF["lcr_frac"] = PF.seq.apply(lambda s: float(lcr_mask(s).mean()) if len(s) else 0)
def maxrun(s, aa):
    best = cur = 0
    for c in s:
        cur = cur + 1 if c == aa else 0; best = max(best, cur)
    return best
PF["polyN_max"] = PF.seq.apply(lambda s: maxrun(s, "N"))
PF["is_AP2"] = PF.prod_desc.str.contains(r"\bap2|apiap2", regex=True)
PF["is_PUF_RNA"] = PF.prod_desc.str.contains(r"puf|rna[- ]binding|rrm", regex=True)
PF["is_CHROM"] = PF.prod_desc.str.contains(r"chromatin|histone|acetyltransferase|deacetylase|bromodomain|remodel", regex=True)
PF["is_regulator"] = PF.is_AP2 | PF.is_PUF_RNA | PF.is_CHROM
PF.to_csv(f"{OUT}/M4RD_D5_pf_features.tsv", sep="\t", index=False)
q90 = PF.asn_frac.quantile(0.9)
PF["asn_rich"] = PF.asn_frac >= q90
erows = []
for name, flag in [("AP2", PF.is_AP2), ("PUF_RNA", PF.is_PUF_RNA), ("CHROM", PF.is_CHROM), ("regulator_any", PF.is_regulator)]:
    a = int((flag & PF.asn_rich).sum()); b = int((~flag & PF.asn_rich).sum())
    c = int((flag & ~PF.asn_rich).sum()); d = int((~flag & ~PF.asn_rich).sum())
    orr, p = st.fisher_exact([[a, b], [c, d]])
    # length+LCR stratified CMH (tertiles of length x median split LCR)
    lbin = pd.qcut(PF.length, 3, labels=False, duplicates="drop")
    bbin = (PF.lcr_frac > PF.lcr_frac.median()).astype(int)
    cmh_tables = []
    for (l, bb), idx in PF.groupby([lbin, bbin]).groups.items():
        g = PF.loc[idx]
        cmh_tables.append([[int((flag.loc[idx] & g.asn_rich).sum()), int(((~flag.loc[idx]) & g.asn_rich).sum())],
                           [int((flag.loc[idx] & ~g.asn_rich).sum()), int(((~flag.loc[idx]) & ~g.asn_rich).sum())]])
    try:
        cmh = st.cochran_mantel_haenszel(cmh_tables) if hasattr(st, "cochran_mantel_haenszel") else None
        cmh_stat, cmh_p = (float(cmh.statistic), float(cmh.pvalue)) if cmh is not None else (np.nan, np.nan)
    except Exception:
        cmh_stat, cmh_p = np.nan, np.nan
    # counterfactual: within LCR-frac deciles, regulator vs background Asn (MW pooled z via stratified permutation)
    PF["lcr_dec"] = pd.qcut(PF.lcr_frac, 10, labels=False, duplicates="drop")
    diffs = []
    for dec, idx in PF.groupby("lcr_dec").groups.items():
        g = PF.loc[idx]
        if flag.loc[idx].sum() >= 3 and (~flag.loc[idx]).sum() >= 10:
            diffs.append(g[flag.loc[idx]].asn_frac.median() - g[(~flag.loc[idx])].asn_frac.median())
    erows.append({"set": name, "n": int(flag.sum()), "n_asn_rich": a,
                  "OR_crude": float(orr), "fisher_p": float(p),
                  "CMH_stat": cmh_stat, "CMH_p": cmh_p,
                  "median_withinLCRbin_diff": float(np.median(diffs)) if diffs else np.nan,
                  "n_bins_used": len(diffs)})
E5 = pd.DataFrame(erows)
E5.to_csv(f"{OUT}/M4RD_D5_regulator_enrichment.tsv", sep="\t", index=False)
print(E5.to_string(), flush=True)

# ---------- D6: natural experiment pilot (AP2/PUF OGs: domain anchor + LCR feature) ----------
# map AP2/PUF Pf genes -> OGs via SP001 protein headers
sp1_path = faa["SP001"]
pf2pid = {}
for line in open(sp1_path):
    if line.startswith(">"):
        m = re.search(r"\[locus_tag=(PF3D7_\d+)\]", line)
        if m: pf2pid[m.group(1)] = line[1:].split()[0]
pid2og = dict(OG[OG.species_id == "SP001"].set_index("protein_id").orthogroup_id)
ap2_genes = PF[PF.is_AP2].gene.tolist()
puf_genes = PF[PF.is_PUF_RNA].gene.tolist()
LOGC = ["SP001", "SP002", "SP003"]  # low-GC pole (Laverania)
HIGC = ["SP004", "SP005", "SP006", "SP007"]  # high-GC pole (vivax clade)
nat = []
for g in (ap2_genes + puf_genes):
    pid = pf2pid.get(g)
    og = pid2og.get(pid)
    if og is None or og not in PASS_OGS: continue
    aln = parse_aln(f"{ALN}/{og}.aln")
    if set(SPECIES) - set(aln): continue
    L = len(aln[SPECIES[0]])
    # conserved anchor: columns with consensus>=0.8 across all 16
    anchor = []
    for i in range(L):
        colaa = [aln[s][i] for s in SPECIES if aln[s][i] != "-"]
        if len(colaa) >= 12 and max(Counter(colaa).values()) / len(colaa) >= 0.8: anchor.append(i)
    # per-species LCR fraction + Asn fraction of the aligned sequence (ungapped)
    feats = {}
    for s in SPECIES:
        sq = aln[s].replace("-", "")
        if not sq: continue
        feats[s] = (float(lcr_mask(aln[s]).mean()), sq.count("N") / len(sq), len(sq))
    lo_lcr = np.mean([feats[s][0] for s in LOGC if s in feats])
    hi_lcr = np.mean([feats[s][0] for s in HIGC if s in feats])
    lo_asn = np.mean([feats[s][1] for s in LOGC if s in feats])
    hi_asn = np.mean([feats[s][1] for s in HIGC if s in feats])
    nat.append({"gene": g, "og": og, "alen": L, "frac_anchor_cols": len(anchor) / L,
                "lowGCpole_LCR": lo_lcr, "highGCpole_LCR": hi_lcr,
                "lowGCpole_Asn": lo_asn, "highGCpole_Asn": hi_asn,
                "pole_LCR_diff": lo_lcr - hi_lcr, "pole_Asn_diff": lo_asn - hi_asn})
NAT = pd.DataFrame(nat)
NAT.to_csv(f"{OUT}/M4RD_D6_natural_experiment.tsv", sep="\t", index=False)
print("D6 OGs:", len(NAT), flush=True)
if len(NAT): print(NAT.to_string(), flush=True)

# ---------- D4 light: essentiality constraint context ----------
try:
    import openpyxl
    w = openpyxl.load_workbook(f"{ROOT}/data/manual_inbox/M4R_ESSENTIALITY/elsworth_S387_2025/NIHMS2082619-supplement-Data_S3.xlsx",
                               read_only=True, data_only=True)
    ws = w["S3a-Essentiality&Fitness.Score"]
    rows = list(ws.iter_rows(values_only=True))
    hdr2 = [str(c) for c in rows[3]]
    gid = next(i for i, c in enumerate(hdr2) if "gene" in c.lower())
    # find essentiality-call column
    callcol = next((i for i, c in enumerate(hdr2) if "essent" in c.lower()), None)
    ess = {}
    for r in rows[4:]:
        if r[gid]:
            g = str(r[gid]).strip()
            ess[g] = str(r[callcol]).strip() if callcol is not None and r[callcol] else ""
    # map Pk genes? Elsworth is P. knowlesi; use only as background rate reference via orthologs — record counts only
    with open(f"{OUT}/M4RD_D4_essentiality_note.txt", "w") as f:
        f.write(f"Elsworth Data_S3 S3a rows={len(ess)} callcol={hdr2[callcol] if callcol is not None else None}\n")
        f.write("Pk quotes: " + str(Counter(ess.values()).most_common(10)) + "\n")
        f.write("NOTE: Pk gene IDs; Pf-regulator constraint mapping via orthology deferred to round 2.\n")
    print("D4 essentiality rows:", len(ess), Counter(ess.values()).most_common(5), flush=True)
except Exception as e:
    open(f"{OUT}/M4RD_D4_essentiality_note.txt", "w").write(f"FAILED: {e}\n")
    print("D4 essentiality FAILED:", e, flush=True)

# ---------- manifest / params ----------
files = sorted([f for f in os.listdir(OUT) if f.startswith("M4RD_")])
mani = pd.DataFrame({"file": files, "sha256": [sha(f"{OUT}/{f}") for f in files]})
mani.to_csv(f"{OUT}/input_manifest.tsv", sep="\t", index=False)
params = {"seed": SEED, "tree": f"{L2}/supermatrix.treefile",
          "ogs_pass": len(PASS_OGS), "sites_parsimony": int(n_sites_run),
          "lcr_rule": "window64 step8 entropy<1.5; short-protein whole-seq",
          "coupling": "q_bh<0.05 vs q>=0.1 (L2v2)", "branch_filter": f"|dGC|>=median({med_abs_dgc:.4f})"}
open(f"{OUT}/params.yaml", "w").write("\n".join(f"{k}: {v}" for k, v in params.items()) + "\n")
open(f"{OUT}/README.md", "w").write(
    "# M4R-D ROUND1\n\nBranch-aware L2 pilot (Fitch parsimony + squared-change GC) + "
    "architecture partition (window64 entropy LCR) + matched adaptive-site controls + "
    "Pf-wide regulator enrichment with LCR counterfactual + cross-pole natural experiment pilot + "
    "essentiality note. See M4R-D_ROUND1_REPORT.md and claim_impact.md.\n")
print("WROTE", files, flush=True)
