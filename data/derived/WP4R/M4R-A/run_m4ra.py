#!/usr/bin/env python
"""M4R-A ROUND1: formal enrichment + IDR counterfactual + lifecycle mapping + ChIP + GCN5 inventory."""
import gzip, hashlib, math, os, re, tarfile
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4R/M4R-A"
os.makedirs(OUT, exist_ok=True)
GFF = f"{ROOT}/data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7.gff"
FAA = f"{ROOT}/data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7_AnnotatedProteins.fasta"

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

def gkey(s):
    m = re.search(r"PF3D7_\d+", str(s))
    return m.group(0) if m else None

# ---- A. protein features (longest isoform per gene) ----
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
print("proteins:", len(prots), flush=True)

AA20 = "ACDEFGHIKLMNPQRSTVWY"
def feats(s):
    n = len(s)
    asn = s.count("N") / n if n else 0.0
    # max poly-Asn tract
    mt, cur = 0, 0
    for ch in s:
        cur = cur + 1 if ch == "N" else 0
        mt = max(mt, cur)
    # LCR proxy: window-20 Shannon entropy < 1.8
    w, th = 20, 1.8
    low = 0
    if n >= w:
        arr = np.frombuffer(s.encode(), dtype=np.uint8)
        cnt = np.zeros(256, int)
        for c in arr[:w]: cnt[c] += 1
        def ent():
            p = cnt[cnt > 0] / w
            return float(-(p * np.log2(p)).sum())
        e = ent(); low = int(e < th)
        for i in range(w, n):
            cnt[arr[i]] += 1; cnt[arr[i - w]] -= 1
            p = cnt[cnt > 0] / w
            if float(-(p * np.log2(p)).sum()) < th: low += 1
        lcr = low / (n - w + 1)
    else:
        lcr = 0.0
    # homopolymer-run fraction (>=4 same AA)
    run, i = 0, 0
    while i < n:
        j = i
        while j < n and s[j] == s[i]: j += 1
        if j - i >= 4: run += j - i
        i = j
    return n, asn, mt, lcr, run / n if n else 0.0

P = pd.DataFrame({g: feats(s) for g, s in prots.items()},
                 index=["length", "asn_frac", "polyn_max", "lcr_frac", "hprun_frac"]).T
P.index.name = "gene"
P = P.reset_index()
print(P.describe().to_string(), flush=True)

# ---- gene coords + products from GFF ----
genes, prod = {}, {}
for line in open(GFF):
    if line.startswith("#") or "\tprotein_coding_gene\t" not in line: continue
    f = line.rstrip("\n").split("\t")
    m = re.search(r"ID=(PF3D7_\d+)", f[8]); d = re.search(r"description=([^;]+)", f[8])
    if m:
        g = m.group(1)
        genes[g] = (f[0], int(f[3]), int(f[4]), f[6])
        prod[g] = (d.group(1) if d else "").replace("%2C", ", ").replace("%27", "'").lower()
D = P.merge(pd.DataFrame({"gene": list(genes), "product": [prod.get(g, "") for g in genes]}),
            on="gene", how="inner")
print("genes with coords+protein:", len(D), flush=True)
q90 = D.asn_frac.quantile(0.9)
D["asn_top10"] = D.asn_frac >= q90
print("asn q90:", q90, flush=True)

# ---- B. curated sets (description-derived; GFF carries no GO terms) ----
SETS = {
    "ApiAP2": r"\bap2|apiap2",
    "PUF": r"puf1|puf2|pumilio",
    "RNA_binding_broad": r"rna[- ]binding|rrm|helicase|zinc finger|alba|dozi|cith|caf\d|not\d|ccr4",
    "chromatin_reg": r"chromatin|histone|bromodomain|set domain|gcn5|sir2|hdac|acetyltransferase|methyltransferase|morne?c|hp1",
    "CCR4_NOT": r"ccr4|caf1|not1|not2|not5|pop2|caf40",
    "transcription_reg": r"transcription factor|transcriptional|rna polymerase|general transcription|mediator",
    "translation_reg": r"translation|elongation factor|initiation factor|eif|ribosome biogenesis",
    "sexual_gametocyte": r"gametocyte|gamete|oocyst|ookinete|zygote|fertili|p47|p48/45|ap2-g|ap2g|gdv1|mdv1",
    "sporozoite_liver": r"sporozoite|liver|uis|csp|trap|ssp|linup|lsa|exp1",
    "proteostasis": r"heat shock|chaperone|dnaj|hsp|proteasome|ubiquitin|unfolded",
}
for k, pat in SETS.items():
    D[k] = D["product"].str.contains(pat, regex=True)

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)

rows = []
for k in SETS:
    hit = D[k].values; top = D.asn_top10.values
    a = int((hit & top).sum()); b = int((~hit & top).sum())
    c = int((hit & ~top).sum()); d = int((~hit & ~top).sum())
    oddsr, p = stats.fisher_exact([[a, b], [c, d]])
    # Woolf 95% CI
    if min(a, b, c, d) > 0:
        se = math.sqrt(1/a + 1/b + 1/c + 1/d)
        lo, hi = math.exp(math.log(oddsr) - 1.96*se), math.exp(math.log(oddsr) + 1.96*se)
    else:
        lo, hi = (0.0, float("inf")) if oddsr in (0, float("inf")) else (float("nan"), float("nan"))
    rows.append({"set": k, "n_set": int(hit.sum()), "n_set_in_top10": a,
                 "odds_ratio": float(oddsr), "ci95_lo": lo, "ci95_hi": hi, "p_fisher": float(p)})
E = pd.DataFrame(rows)
E["q_bh"] = bh(E.p_fisher.values)
E = E.sort_values("p_fisher")
E.to_csv(f"{OUT}/M4RA_formal_enrichment.tsv", sep="\t", index=False)
print(E.to_string(), flush=True)

# ---- C. IDR/LCR counterfactual ----
D["loglen"] = np.log(D.length)
D["is_reg"] = (D.ApiAP2 | D.PUF | D.chromatin_reg | D.CCR4_NOT).astype(int)
X = sm.add_constant(D[["is_reg", "loglen", "lcr_frac"]])
mdl = sm.Logit(D.asn_top10.astype(int), X).fit(disp=0)
adj = pd.DataFrame({"coef": mdl.params, "se": mdl.bse,
                    "OR": np.exp(mdl.params),
                    "OR_lo": np.exp(mdl.params - 1.96*mdl.bse),
                    "OR_hi": np.exp(mdl.params + 1.96*mdl.bse),
                    "p": mdl.pvalues})
adj.to_csv(f"{OUT}/M4RA_counterfactual_logit.tsv", sep="\t")
print(adj.to_string(), flush=True)
# stratified MW within LCR tertiles: regulators vs rest on asn_frac
pos = D.lcr_frac[D.lcr_frac > 0]
cut = pos.median()
D["lcr_tert"] = pd.Categorical(
    np.where(D.lcr_frac == 0, "zero", np.where(D.lcr_frac <= cut, "lowpos", "high")),
    categories=["zero", "lowpos", "high"])
srows = []
for t, gdf in D.groupby("lcr_tert"):
    r = gdf[gdf.is_reg == 1].asn_frac; o = gdf[gdf.is_reg == 0].asn_frac
    u, p = stats.mannwhitneyu(r, o, alternative="two-sided")
    srows.append({"lcr_tertile": t, "n_reg": len(r), "n_rest": len(o),
                  "med_reg": float(r.median()), "med_rest": float(o.median()), "MW_p": float(p)})
S = pd.DataFrame(srows)
S.to_csv(f"{OUT}/M4RA_counterfactual_stratified.tsv", sep="\t", index=False)
print(S.to_string(), flush=True)
# Asn-rich vs poly-Asn distinction
D["has_polyn10"] = D.polyn_max >= 10
print(pd.crosstab(D.asn_top10, D.has_polyn10), flush=True)
print("frac Asn-top10 with polyN>=10:", float(((D.asn_top10) & (D.has_polyn10)).sum() / D.asn_top10.sum()), flush=True)

# ---- D. GSE75795 lifecycle mapping ----
tar = tarfile.open(f"{ROOT}/data/raw/geo/GSE75795/GSE75795_RAW.tar")
want = {"GSM1967965": "male", "GSM1967966": "female", "GSM1967967": "male2"}
counts = {}
for m in tar.getmembers():
    pre = m.name.split("_")[0]
    if pre not in want or "_9.1." not in m.name or ".raw." not in m.name: continue
    strand = "plus" if ".plus." in m.name else "minus"
    import io
    f = gzip.open(tar.extractfile(m), "rt")
    cov = {}
    for line in f:
        if line.startswith("track"): continue
        c = line.split("\t")
        if len(c) < 4: continue
        cov.setdefault(c[0], []).append((int(c[1]), int(c[2]), float(c[3])))
    f.close()
    counts[(pre, strand)] = cov
    print(pre, strand, "intervals:", sum(len(v) for v in cov.values()), flush=True)
# gene-level sums
G = pd.DataFrame({"gene": list(genes)})
for k, v in genes.items():
    pass
gdf = pd.DataFrame([{"gene": g, "chrom": v[0], "start": v[1], "end": v[2]} for g, v in genes.items()])
samp = {}
for (pre, strand), cov in counts.items():
    tot = np.zeros(len(gdf))
    per_chr = gdf.groupby("chrom").indices
    for ch, idx in per_chr.items():
        ivs = cov.get(ch)
        if not ivs: continue
        starts = np.array([a for a, b, v in ivs]); ends = np.array([b for a, b, v in ivs]); vals = np.array([v for a, b, v in ivs])
        for ii in idx:
            s, e = gdf.start.iloc[ii], gdf.end.iloc[ii]
            o = np.clip(np.minimum(ends, e) - np.maximum(starts, s), 0, None)
            tot[ii] = float((o * vals).sum())
    samp.setdefault(pre, []).append(tot)
gene_counts = pd.DataFrame({"gene": gdf.gene})
for pre, arrs in samp.items():
    gene_counts[pre] = arrs[0] + arrs[1]
lib = gene_counts[["GSM1967965", "GSM1967966", "GSM1967967"]].values
cpm = gene_counts[["GSM1967965", "GSM1967966", "GSM1967967"]].div(lib.sum(axis=0) / 1e6, axis=1)
gene_counts["male_cpm"] = (cpm.GSM1967965 + cpm.GSM1967967) / 2
gene_counts["female_cpm"] = cpm.GSM1967966
gene_counts["log2FC_FvsM"] = np.log2((gene_counts.female_cpm + 0.5) / (gene_counts.male_cpm + 0.5))
gene_counts.to_csv(f"{OUT}/M4RA_gametocyte_counts.tsv", sep="\t", index=False)
M = gene_counts.merge(D[["gene", "asn_frac", "asn_top10", "lcr_frac", "is_reg"] + list(SETS)], on="gene")
M["expressed"] = (M.male_cpm + M.female_cpm) > 2
Me = M[M.expressed].copy()
print("expressed genes:", len(Me), flush=True)
rho, prho = stats.spearmanr(Me.asn_frac, Me.log2FC_FvsM)
u, pmw = stats.mannwhitneyu(Me[Me.asn_top10].log2FC_FvsM, Me[~Me.asn_top10].log2FC_FvsM, alternative="two-sided")
lrows = [{"test": "spearman(asn_frac, F-vs-M logFC)", "stat": float(rho), "p": float(prho), "n": len(Me)},
         {"test": "MW Asn-top10 vs rest on F-vs-M logFC", "stat": float(u), "p": float(pmw),
          "n": len(Me), "med_top10": float(Me[Me.asn_top10].log2FC_FvsM.median()),
          "med_rest": float(Me[~Me.asn_top10].log2FC_FvsM.median())}]
for k in SETS:
    hit = Me[Me[k]].log2FC_FvsM; rst = Me[~Me[k]].log2FC_FvsM
    if len(hit) >= 3:
        u2, p2 = stats.mannwhitneyu(hit, rst, alternative="two-sided")
        lrows.append({"test": f"MW {k} vs rest on F-vs-M logFC", "stat": float(u2), "p": float(p2),
                      "n": len(Me), "med_top10": float(hit.median()), "med_rest": float(rst.median())})
L = pd.DataFrame(lrows)
L["q_bh"] = bh(L.p.values)
L.to_csv(f"{OUT}/M4RA_lifecycle_mapping.tsv", sep="\t", index=False)
print(L.to_string(), flush=True)

# ---- E. AP2-G ChIP trio ----
chip_files = {
    "GSE120448_AP2-G_S": ("GSE120448", ["GSM3400933_AP2-G_S_rep1.bed.gz", "GSM3400935_AP2-G_S_rep2.bed.gz"]),
    "GSE120448_AP2-G_R": ("GSE120448", ["GSM3400937_AP2-G_R_rep1.bed.gz", "GSM3400939_AP2-G_R_rep2.bed.gz"]),
    "GSE120448_AP2-G_G": ("GSE120448", ["GSM3400941_AP2-G_G_rep1.bed.gz", "GSM3400943_AP2-G_G_rep2.bed.gz"]),
    "GSE134268_NCC": ("GSE134268", ["GSM3940782_NCC_peaks_rep1.bed.gz", "GSM3940784_NCC_peaks_rep2.bed.gz"]),
    "GSE134268_SCC": ("GSE134268", ["GSM3940786_SCC_peaks_rep1.bed.gz", "GSM3940788_SCC_peaks_rep2.bed.gz"]),
    "GSE120488_AP2-G_schiz": ("GSE120488", ["GSM3401472_AP2-G_rep1.bed.gz", "GSM3401473_AP2-G_rep2.bed.gz"]),
    "GSE120488_AP2-I_schiz": ("GSE120488", ["GSM3401470_AP2-I_rep1.bed.gz", "GSM3401471_AP2-I_rep2.bed.gz"]),
}
def load_bed(gse, fn):
    out = []
    with gzip.open(f"{ROOT}/data/raw/geo/{gse}/{fn}", "rt") as f:
        for line in f:
            if line.startswith("track"): continue
            c = line.split("\t")
            out.append((c[0], int(c[1]), int(c[2])))
    return out

def consensus(p1, p2):
    # keep p1 peaks overlapping any p2 peak (same chrom)
    by = {}
    for c, s, e in p2: by.setdefault(c, []).append((s, e))
    keep = []
    for c, s, e in p1:
        for s2, e2 in by.get(c, []):
            if min(e, e2) > max(s, s2):
                keep.append((c, s, e)); break
    return keep

glist = gdf  # chrom,start,end,gene
bound = {}
crow = []
for cond, (gse, (f1, f2)) in chip_files.items():
    p1, p2 = load_bed(gse, f1), load_bed(gse, f2)
    con = consensus(p1, p2)
    # map to genes: peak within gene body or <=1.5kb upstream of start (strand-aware approx: just flank)
    bg = set()
    per_chr = glist.groupby("chrom").indices
    starts = {c: glist.start.iloc[i].values for c, i in per_chr.items()}
    ends = {c: glist.end.iloc[i].values for c, i in per_chr.items()}
    gnames = {c: glist.gene.iloc[i].values for c, i in per_chr.items()}
    for c, s, e in con:
        if c not in per_chr: continue
        ov = (starts[c] - 1500 <= e) & (ends[c] + 500 >= s)
        for g in gnames[c][ov]: bg.add(g)
    bound[cond] = bg
    crow.append({"condition": cond, "peaks_rep1": len(p1), "peaks_rep2": len(p2),
                 "consensus": len(con), "n_bound_genes": len(bg)})
C = pd.DataFrame(crow)
C.to_csv(f"{OUT}/M4RA_chip_consensus.tsv", sep="\t", index=False)
print(C.to_string(), flush=True)
# bound genes: Asn-rich regulator enrichment (union AP2-G conditions + AP2-I)
D["ap2g_bound"] = D.gene.isin(set().union(*[bound[k] for k in bound if "AP2-G" in k or "NCC" in k or "SCC" in k]))
D["ap2i_bound"] = D.gene.isin(bound["GSE120488_AP2-I_schiz"])
hrows = []
for col in ["ap2g_bound", "ap2i_bound"]:
    b = D[D[col]]
    hrows.append({"bound_set": col, "n_bound": len(b),
                  "med_asn_bound": float(b.asn_frac.median()),
                  "med_asn_unbound": float(D[~D[col]].asn_frac.median()),
                  "MW_p": float(stats.mannwhitneyu(b.asn_frac, D[~D[col]].asn_frac, alternative="two-sided")[1])})
    for k in ["ApiAP2", "PUF", "sexual_gametocyte", "asn_top10"]:
        kk = D.asn_top10 if k == "asn_top10" else D[k]
        a = int((D[col] & kk).sum()); bb = int(((~D[col]) & kk).sum())
        cc = int((D[col] & ~kk).sum()); dd = int(((~D[col]) & ~kk).sum())
        orr, pp = stats.fisher_exact([[a, bb], [cc, dd]])
        hrows.append({"bound_set": f"{col} x {k}", "n_bound": a, "med_asn_bound": float(orr),
                      "med_asn_unbound": float("nan"), "MW_p": float(pp)})
H = pd.DataFrame(hrows)
H.to_csv(f"{OUT}/M4RA_chip_overlap.tsv", sep="\t", index=False)
print(H.to_string(), flush=True)

# ---- F. GCN5 MOESM inventory ----
import openpyxl
grow = []
sup = f"{ROOT}/data/raw/pmc/M4R_GCN5_2026/suppl/41467_2026_74632_MOESM7_ESM.xlsx"
wb = openpyxl.load_workbook(sup, read_only=True, data_only=True)
for sn in wb.sheetnames:
    ws = wb[sn]
    vals = list(ws.values)
    grow.append({"sheet": sn, "nrows": ws.max_row, "ncols": ws.max_column,
                 "header": str(vals[0])[:300] if vals else ""})
wb2 = openpyxl.load_workbook(f"{ROOT}/data/raw/pmc/M4R_GCN5_2026/suppl/41467_2026_74632_MOESM3_ESM.xlsx",
                             read_only=True, data_only=True)
ws = wb2[wb2.sheetnames[0]]
vals = list(ws.values)
grow.append({"sheet": "MOESM3:" + wb2.sheetnames[0], "nrows": ws.max_row, "ncols": ws.max_column,
             "header": str(vals[0])[:300]})
wb3 = openpyxl.load_workbook(f"{ROOT}/data/raw/pmc/M4R_GCN5_2026/suppl/41467_2026_74632_MOESM4_ESM.xlsx",
                             read_only=True, data_only=True)
ws = wb3[wb3.sheetnames[0]]
vals = list(ws.values)
grow.append({"sheet": "MOESM4:" + wb3.sheetnames[0], "nrows": ws.max_row, "ncols": ws.max_column,
             "header": str(vals[0])[:500] if vals else ""})
G5 = pd.DataFrame(grow)
G5.to_csv(f"{OUT}/M4RA_gcn5_inventory.tsv", sep="\t", index=False)
print(G5.to_string(), flush=True)

# key GCN5 numbers: Figure 2 source data (repeat deletion phenotype)
wb = openpyxl.load_workbook(sup, read_only=False, data_only=True)
ws = wb["Figure 2"]
vals = list(ws.values)
with open(f"{OUT}/M4RA_gcn5_fig2_snapshot.tsv", "w") as f:
    for r in vals[:25]:
        f.write("\t".join("" if v is None else str(v)[:80] for v in r) + "\n")

# manifests / checksums
mani = pd.DataFrame({"file": [GFF, FAA,
    f"{ROOT}/data/derived/WP4/L1d_aat_function/L1d_keyword_enrichment.tsv",
    f"{ROOT}/data/raw/geo/GSE75795/GSE75795_RAW.tar",
    f"{ROOT}/data/raw/pmc/M4R_GCN5_2026/suppl/41467_2026_74632_MOESM7_ESM.xlsx"]})
mani["sha256"] = mani.file.map(sha); mani["size"] = mani.file.map(os.path.getsize)
mani.to_csv(f"{OUT}/input_manifest.tsv", sep="\t", index=False)
open(f"{OUT}/params.yaml", "w").write(
    "task: M4R-A-ROUND1\nasn_top10_rule: >=90th pct of asn_frac (PlasmoDB-71 longest isoform)\n"
    "lcr_proxy: window-20 Shannon entropy <1.8 fraction (heuristic, not InterPro)\n"
    "gene_sets: product-description regex (GFF has no GO terms); documented description-derived\n"
    "gametocyte: GSE75795 raw 9.1 bedGraphs summed over v71 gene bodies; n=1 per sex -> descriptive\n"
    "chip: rep-consensus (overlap) peaks; gene = body or <=1.5kb upstream; v71 coords vs 9.1 mapping approx\n"
    "stats: Fisher+WoolfCI+BH; Logit adjusted; MW stratified; spearman\n")
outs = sorted(f for f in os.listdir(OUT) if os.path.isfile(os.path.join(OUT, f)) and f != "checksums.sha256")
with open(f"{OUT}/checksums.sha256", "w") as fh:
    for f in outs:
        fh.write(f"{sha(os.path.join(OUT, f))}  {f}\n")
print("DONE", flush=True)
