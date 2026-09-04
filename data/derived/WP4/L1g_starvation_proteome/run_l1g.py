#!/usr/bin/env python
# L1g (ring 3b): protein-level test of the Asn axis using Li 2024 iScience starvation proteomics (mmc3)
# - same study as GSE226632 (our TE replication) -> true cross-layer test (mRNA/TE/protein in one experiment)
# - stage-delay confound controlled by excluding schizont-peak genes (from L1f GSE58402 stage profiles)
import hashlib, json, os, re
from collections import Counter
import numpy as np
import pandas as pd
import openpyxl
from scipy import stats

ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4/L1g_starvation_proteome"
os.makedirs(OUT, exist_ok=True)
LI = f"{ROOT}/data/raw/suppl/LI2024_ISCIENCE/mmc3.xlsx"

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

# ---------- features (same CDS rules as L1d/e) ----------
STD = dict(zip([a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
               "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
acc = "GCF_000002765.6"
CDS = f"{ROOT}/data/raw/ncbi-datasets/{acc}/ncbi_dataset/data/{acc}/cds_from_genomic.fna"
best = {}
with open(CDS) as fh:
    hid, seq = None, []
    def flush():
        if hid is None: return
        m = re.search(r"\[locus_tag=([^\]]+)\]", hid)
        if not m: return
        g = re.search(r"PF3D7_\d+", m.group(1))
        if not g: return
        g = g.group(0)
        s = "".join(seq).upper().replace("U", "T")
        if "[partial" in hid or len(s) < 150 or len(s) % 3 or not set(s) <= set("ATGC"): return
        cod = [s[i:i+3] for i in range(0, len(s), 3)]
        if "*" in [STD.get(c, "*") for c in cod][:-1]: return
        if len(s) > len(best.get(g, ("", ""))[1]): best[g] = (hid, s)
    for line in fh:
        if line.startswith(">"): flush(); hid, seq = line.rstrip(), []
        else: seq.append(line.strip())
    flush()
feat = {}
for g, (_, s) in best.items():
    cod = [s[i:i+3] for i in range(0, len(s), 3) if STD.get(s[i:i+3], "*") != "*"]
    if cod:
        cc = Counter(cod); n = len(cod)
        feat[g] = {"AAT_freq": cc["AAT"] / n, "Asn_content": (cc["AAT"] + cc["AAC"]) / n}
FE = pd.DataFrame(feat).T

# ---------- Li protein DE (AA starvation 6h vs CM) ----------
wb = openpyxl.load_workbook(LI, read_only=True)
ws = wb["differential expressed"]
rows = list(ws.iter_rows(min_row=3, values_only=True))
wb.close()
DE = pd.DataFrame([r[:6] for r in rows if r[0] and str(r[0]).startswith("PF3D7")],
                  columns=["gene", "logFC", "logCPM", "F", "p", "FDR"])
DE["gene"] = DE.gene.str.replace(r"\.\d+$", "", regex=True)
for c in ["logFC", "logCPM", "F", "p", "FDR"]: DE[c] = pd.to_numeric(DE[c], errors="coerce")
print("Li protein DE genes:", len(DE), flush=True)

# ---------- stage annotation from L1f (mRNA peak stage) ----------
STAGES = ["TE_2h_ring", "TE_10h_latering", "TE_19h_troph", "TE_31h_latetroph", "TE_46h_schizont"]
try:
    TEf = pd.read_csv(f"{ROOT}/data/derived/WP4/L1f_riboseq_te/L1f_te_by_stage.tsv", sep="\t", index_col=0)
    TEf.index.name = "gene"
    TEf["peak_stage"] = TEf[STAGES].idxmax(axis=1).str.replace("TE_", "")
    stage_map = TEf["peak_stage"]
except Exception as e:
    print("stage map failed:", e); stage_map = pd.Series(dtype=str)

# ---------- dTE from L1 replication (same study, translation layer) ----------
DTE = pd.read_csv(f"{ROOT}/data/derived/WP4/L1_codon_x_stress/L1_GSE226632_dTE.tsv", sep="\t").set_index("gene")

M = DE.set_index("gene").join(FE).join(DTE).join(stage_map.rename("peak_stage"))
M = M.dropna(subset=["logFC", "Asn_content"])
q90 = FE.Asn_content.quantile(0.9)  # genome-wide threshold (not detected-set re-derived)
M["AsnTop10"] = M.Asn_content >= q90
print("merged:", len(M), "AsnTop10 detected:", int(M.AsnTop10.sum()), flush=True)

def block(sub, tag):
    hi, lo = sub[sub.AsnTop10], sub[~sub.AsnTop10]
    sig = sub[sub.FDR < 0.05]
    sig_hi_up = ((sig.AsnTop10) & (sig.logFC > 0)).sum(); sig_lo_up = ((~sig.AsnTop10) & (sig.logFC > 0)).sum()
    sig_hi_dn = ((sig.AsnTop10) & (sig.logFC < 0)).sum(); sig_lo_dn = ((~sig.AsnTop10) & (sig.logFC < 0)).sum()
    r = {"set": tag, "n": len(sub), "n_AsnTop10": len(hi),
         "median_logFC_AsnTop10": float(hi.logFC.median()), "median_logFC_rest": float(lo.logFC.median()),
         "MW_p_two_sided": float(stats.mannwhitneyu(hi.logFC, lo.logFC).pvalue),
         "spearman_Asn_logFC": float(stats.spearmanr(sub.Asn_content, sub.logFC).statistic),
         "spearman_p": float(stats.spearmanr(sub.Asn_content, sub.logFC).pvalue),
         "frac_AsnTop10_in_sigUP": float(sig_hi_up / max(sig_hi_up + sig_lo_up, 1)),
         "frac_AsnTop10_in_sigDOWN": float(sig_hi_dn / max(sig_hi_dn + sig_lo_dn, 1)),
         "n_sig": len(sig)}
    return r

res = [block(M, "all_detected")]
if "peak_stage" in M:
    Mns = M[M.peak_stage != "46h_schizont"]
    res.append(block(Mns, "excl_schizont_peak"))
# cross-layer: protein logFC vs dTE (same study)
X = M.dropna(subset=["dTE"])
res.append({"set": "crosslayer_all", "n": len(X),
            "spearman_protLogFC_vs_dTE": float(stats.spearmanr(X.logFC, X.dTE).statistic),
            "spearman_p": float(stats.spearmanr(X.logFC, X.dTE).pvalue)})
Xh = X[X.AsnTop10]
res.append({"set": "crosslayer_AsnTop10", "n": len(Xh),
            "spearman_protLogFC_vs_dTE": float(stats.spearmanr(Xh.logFC, Xh.dTE).statistic) if len(Xh) > 5 else np.nan,
            "spearman_p": float(stats.spearmanr(Xh.logFC, Xh.dTE).pvalue) if len(Xh) > 5 else np.nan})
R = pd.DataFrame(res)
R.to_csv(f"{OUT}/L1g_tests.tsv", sep="\t", index=False)
M.reset_index().to_csv(f"{OUT}/L1g_merged_table.tsv", sep="\t", index=False)
print(R.to_string(), flush=True)

json.dump(res, open(f"{OUT}/L1g_qc.json", "w"), indent=2, default=str)
open(f"{OUT}/claim_impact.md", "w").write("""# L1g claim impact (ring 3b: starvation proteomics; exploratory)

- date: 2026-09-04
- data: Li 2024 iScience (PMC11544085) Table S2, NF54 6h AA starvation vs CM, 1148 proteins DE table (limma).
  SAME study as our TE replication dataset GSE226632 -> mRNA/TE/protein cross-layer test within one experiment.
- tests: Asn-top10% (genome-wide threshold) vs rest on protein logFC (MW, two-sided); spearman(Asn, logFC);
  sig-set direction fractions; stage-confound control = exclude schizont-peak genes (L1f annotation);
  cross-layer spearman(protein logFC, dTE).
- results: L1g_tests.tsv. Verdict after user review; honest reading: protein layer at 6h still integrates
  slow abundance dynamics; stage-delay confound documented by Li themselves (56/78 down = merozoite proteins).
- claim impact: direct protein-layer evidence for/against Asn axis under starvation; no claim change yet.
""")
mani = pd.DataFrame({"file": [LI, CDS, f"{ROOT}/data/derived/WP4/L1_codon_x_stress/L1_GSE226632_dTE.tsv",
                              f"{ROOT}/data/derived/WP4/L1f_riboseq_te/L1f_te_by_stage.tsv"]})
mani["sha256"] = mani.file.map(sha); mani["size"] = mani.file.map(os.path.getsize)
mani.to_csv(f"{OUT}/input_manifest.tsv", sep="\t", index=False)
outs = [f for f in os.listdir(OUT) if f != "checksums.sha256" and os.path.isfile(f"{OUT}/{f}")]
with open(f"{OUT}/checksums.sha256", "w") as fh:
    for f in sorted(outs):
        fh.write(f"{sha(os.path.join(OUT, f))}  {f}\n")
print("DONE", flush=True)
