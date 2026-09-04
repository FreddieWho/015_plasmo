#!/usr/bin/env python
# L1d (ring 2): are AAT-rich genes "a gang"? keyword enrichment on PlasmoDB product descriptions
import hashlib, json, os, re
from collections import Counter
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4/L1d_aat_function"
os.makedirs(OUT, exist_ok=True)

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

def gkey(s):
    m = re.search(r"PF3D7_\d+", str(s))
    return m.group(0) if m else None

# AAT freq per gene (same CDS rules as L1)
BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
STD = dict(zip(CODONS, AA))
acc = "GCF_000002765.6"
path = f"{ROOT}/data/raw/ncbi-datasets/{acc}/ncbi_dataset/data/{acc}/cds_from_genomic.fna"
best = {}
with open(path) as fh:
    hid, seq = None, []
    def flush():
        if hid is None: return
        m = re.search(r"\[locus_tag=([^\]]+)\]", hid)
        if not m: return
        g = gkey(m.group(1))
        if not g: return
        s = "".join(seq).upper().replace("U", "T")
        if "[partial" in hid or len(s) < 150 or len(s) % 3 or not set(s) <= set("ATGC"): return
        cod = [s[i:i+3] for i in range(0, len(s), 3)]
        if "*" in [STD.get(c, "*") for c in cod][:-1]: return
        if len(s) > len(best.get(g, ("", ""))[1]): best[g] = (hid, s)
    for line in fh:
        if line.startswith(">"): flush(); hid, seq = line.rstrip(), []
        else: seq.append(line.strip())
    flush()
aat = {}
for g, (_, s) in best.items():
    cod = [s[i:i+3] for i in range(0, len(s), 3) if STD.get(s[i:i+3], "*") != "*"]
    if cod: aat[g] = cod.count("AAT") / len(cod)
AAT = pd.Series(aat).rename("AAT_freq")
print("genes:", len(AAT), flush=True)

# product descriptions
prod = {}
for line in open(f"{ROOT}/data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7.gff"):
    if line.startswith("#") or "\tprotein_coding_gene\t" not in line: continue
    a = line.rstrip("\n").split("\t")[8]
    m = re.search(r"ID=(PF3D7_\d+)", a); d = re.search(r"description=([^;]+)", a)
    if m: prod[m.group(1)] = (d.group(1) if d else "").replace("%2C", ", ").replace("%27", "'").lower()
D = pd.DataFrame({"gene": AAT.index, "AAT_freq": AAT.values})
D["product"] = D.gene.map(prod).fillna("")
print("with products:", (D["product"] != "").sum(), flush=True)

# keyword sets
KW = {
    "AP2_transcription_factor": r"\bap2|apiap2",
    "RNA_binding_PUF_etc": r"puf|rna[- ]binding|zinc finger|helicase|rrm",
    "ribosomal_protein": r"ribosomal protein|ribosome",
    "translation_machine": r"translation|elongation factor|initiation factor|trna",
    "exported_surface_var_rif": r"exported|erythrocyte membrane|rifin|stevor|\bvar\b|pfemp",
    "kinase": r"kinase",
    "proteasome_ubiquitin": r"ubiquitin|proteasome",
    "heat_shock_chaperone": r"heat shock|chaperone|dnaj",
    "conserved_unknown": r"conserved.*unknown|unknown function|hypothetical",
    "asparagine_rich_antigen": r"asparagine.rich|asn.rich",
}
q90 = D.AAT_freq.quantile(0.9)
D["top10pct"] = D.AAT_freq >= q90
rows = []
for k, pat in KW.items():
    hit = D["product"].str.contains(pat, regex=True)
    a = int((hit & D.top10pct).sum()); b = int((~hit & D.top10pct).sum())
    c = int((hit & ~D.top10pct).sum()); d = int((~hit & ~D.top10pct).sum())
    or_, p = stats.fisher_exact([[a, b], [c, d]])
    rows.append({"keyword": k, "in_top10": a, "out_top10": c, "odds_ratio": or_, "p": p})
R = pd.DataFrame(rows).sort_values("p")
p = R.p.values; n = len(p); o = np.argsort(p)
q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
R["q_bh"] = np.minimum(q, 1)
R.to_csv(f"{OUT}/L1d_keyword_enrichment.tsv", sep="\t", index=False)
print(R.to_string(), flush=True)

# top-40 AAT genes with products (manual inspection aid)
top = D.nlargest(40, "AAT_freq")[["gene", "AAT_freq", "product"]]
top.to_csv(f"{OUT}/L1d_top40_AAT_genes.tsv", sep="\t", index=False)

open(f"{OUT}/claim_impact.md", "w").write(f"""# L1d claim impact (ring 2 of AAT mechanism chain; exploratory)

- date: 2026-09-04
- question: are AAT-rich genes a coherent functional group (a program) or a mixed bag?
- keyword Fisher enrichment (top-decile AAT vs rest), BH q in L1d_keyword_enrichment.tsv.
- top enriched: {R.iloc[0].keyword} (OR={R.iloc[0].odds_ratio:.2f}, q={R.iloc[0].q_bh:.3g})
- AP2 count in top10%: {int(R[R.keyword=='AP2_transcription_factor'].in_top10.iloc[0])}; interpretation after user review.
- claim impact: exploratory functional coherence check; no claim change.
""")
mani = pd.DataFrame({"file": [path, f"{ROOT}/data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7.gff"]})
mani["sha256"] = mani.file.map(sha); mani["size"] = mani.file.map(os.path.getsize)
mani.to_csv(f"{OUT}/input_manifest.tsv", sep="\t", index=False)
outs = [f for f in os.listdir(OUT) if f != "checksums.sha256" and os.path.isfile(os.path.join(OUT, f))]
with open(f"{OUT}/checksums.sha256", "w") as fh:
    for f in sorted(outs):
        fh.write(f"{sha(os.path.join(OUT, f))}  {f}\n")
print("DONE", flush=True)
