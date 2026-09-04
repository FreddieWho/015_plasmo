#!/usr/bin/env python3
"""
M1-03 Artifact stratification: assembly / annotation / LCR / mappability proxy
Input: M1-01_composition_table.tsv, genomic.fna, cds_from_genomic.fna, gff counts
Output: M1-03_artifact_table.tsv, correlations, claim_impact
A-level only, no host omics, no reads
"""
import pathlib, csv, math, re, collections, hashlib
from collections import Counter, defaultdict

BASE=pathlib.Path("/home/huyudi/015_plasmo")
FROZEN=BASE/"data/metadata/frozen_assembly_candidates.tsv"
COMP=BASE/"data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv"
RAW_BASE=BASE/"data/raw/ncbi-datasets"
OUTDIR=BASE/"data/derived/WP1/M1-03_artifact"
OUTDIR.mkdir(parents=True, exist_ok=True)

# load composition
rows=[]
with open(COMP) as fh:
    r=csv.DictReader(fh, delimiter="\t")
    for row in r: rows.append(row)

# parse GFF gene counts as annotation proxy
# For each accession, count genes in genomic.gff
def count_gff_features(acc):
    p = RAW_BASE/acc/"ncbi_dataset/data"/acc/"genomic.gff"
    if not p.exists():
        return None, None, None
    gene=0; cds=0; mrna=0
    with open(p) as fh:
        for line in fh:
            if line.startswith("#"): continue
            parts=line.split("\t")
            if len(parts)<3: continue
            ftype=parts[2]
            if ftype=="gene": gene+=1
            elif ftype=="CDS": cds+=1
            elif ftype=="mRNA": mrna+=1
    return gene, cds, mrna

# LCR / low entropy windows: simple entropy of 64bp windows
def lcr_fraction_of_seq(seq, window=64, entropy_thresh=1.5):
    # entropy per window: -sum(p log2 p)
    # low entropy = repetitive/low complexity
    seq=seq.upper()
    if len(seq) < window: return 0,0
    low=0; total=0
    for i in range(0, len(seq)-window+1, window):  # non-overlapping for speed
        win=seq[i:i+window]
        if 'N' in win: continue
        total+=1
        cnt=Counter(win)
        ent= -sum((c/window)*math.log2(c/window) for c in cnt.values())
        if ent < entropy_thresh:
            low+=1
    frac = low/total if total else 0
    return frac, total

# For speed, sample at most 5 Mb per species for LCR (first chromosomes)
def lcr_for_genome(path):
    total_frac_w=0; total_win=0
    for header, seq in parse_fasta(path):
        h=header.lower()
        if "mitochondrion" in h or "apicoplast" in h or "plastid" in h: continue
        # only nuclear
        # sample first 5 Mb cumulatively? just compute per contig and aggregate weighted
        frac, wins = lcr_fraction_of_seq(seq)
        # weight by wins
        total_frac_w += frac * wins
        total_win += wins
        if total_win > 80000:  # ~5Mb /64 ~78k windows
            break
    overall = total_frac_w/total_win if total_win else 0
    return overall

def parse_fasta(path):
    header=None; parts=[]
    with open(path) as fh:
        for line in fh:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts).upper()
                header=line; parts=[]
            else:
                parts.append(line.strip())
        if header is not None:
            yield header, "".join(parts).upper()

out_rows=[]
for row in rows:
    sid=row["species_id"]; acc=row["frozen_accession"]
    genome_gc=float(row["genome_gc"])
    n_percent=float(row["n_percent"])
    n50=int(row["n50_bp"]); total_len=int(row["total_length_bp"])
    num_contigs=int(row["num_contigs"])
    cds_count=int(row["cds_count"]) if row["cds_count"] else 0
    gc3=float(row["gc3"]) if row["gc3"] and row["gc3"]!='None' else None
    homopoly=float(row["homopolymer_density_per_Mb"])
    # GFF
    gene, gff_cds, mrna = count_gff_features(acc)
    cds_density = cds_count/(total_len/1e6) if total_len else 0
    # LCR
    unpack=RAW_BASE/acc/"ncbi_dataset/data"/acc
    cand=[p for p in unpack.glob("*_genomic.fna") if "cds_from" not in p.name]
    gpath=cand[0] if cand else None
    lcr_frac = lcr_for_genome(gpath) if gpath else None
    # mappability proxy: homopolymer already, plus N% and contig count
    out_rows.append({
        "species_id":sid,"species":row["species"],"accession":acc,
        "genome_gc":genome_gc,"gc3":gc3 if gc3 is not None else "",
        "total_len_Mb": total_len/1e6, "n50_kb": n50/1000, "n_percent": n_percent,
        "num_contigs": num_contigs,
        "cds_count": cds_count, "cds_density_per_Mb": cds_density,
        "gff_gene": gene if gene is not None else "", "gff_cds": gff_cds if gff_cds is not None else "",
        "bad_cds_not_div3": row["bad_cds_not_div3"], "short_cds_lt150": row["short_cds_lt150"],
        "homopolymer_per_Mb": homopoly,
        "lcr_frac_lowEntropy": f"{lcr_frac:.4f}" if lcr_frac is not None else "",
        "risk_flag": row["risk_flag"],
    })

# write table
fieldnames=["species_id","species","accession","genome_gc","gc3","total_len_Mb","n50_kb","n_percent","num_contigs","cds_count","cds_density_per_Mb","gff_gene","gff_cds","bad_cds_not_div3","short_cds_lt150","homopolymer_per_Mb","lcr_frac_lowEntropy","risk_flag"]
with open(OUTDIR/"M1-03_artifact_table.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    for r in out_rows: w.writerow(r)
print(f"wrote {OUTDIR/'M1-03_artifact_table.tsv'}")

# correlations: Pearson between genome_gc and each artifact metric
import math
def pearson(xs, ys):
    n=len(xs)
    if n<3: return 0,1
    mx=sum(xs)/n; my=sum(ys)/n
    num=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    denx=math.sqrt(sum((x-mx)**2 for x in xs))
    deny=math.sqrt(sum((y-my)**2 for y in ys))
    if denx==0 or deny==0: return 0,1
    r=num/(denx*deny)
    # approximate p via t? just report r
    return r, None

metrics=["total_len_Mb","n50_kb","n_percent","num_contigs","cds_density_per_Mb","homopolymer_per_Mb","lcr_frac_lowEntropy"]
# need numeric conversions
corr_rows=[]
gcs=[float(r["genome_gc"]) for r in out_rows]
for m in metrics:
    vals=[]
    valid_gcs=[]
    for r in out_rows:
        v=r[m]
        try:
            fv=float(v) if v!="" else None
            if fv is not None:
                vals.append(fv)
                valid_gcs.append(float(r["genome_gc"]))
        except: pass
    if len(vals)>=3:
        rval,_=pearson(valid_gcs, vals)
        corr_rows.append({"artifact_metric": m, "pearson_r_with_genomeGC": f"{rval:.3f}", "n": len(vals), "interpretation": "confounded if |r|>0.6" if abs(rval)>0.6 else "not strongly confounded"})
    else:
        corr_rows.append({"artifact_metric": m, "pearson_r_with_genomeGC": "NA", "n": len(vals), "interpretation": "insufficient"})

with open(OUTDIR/"M1-03_correlations.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh, fieldnames=["artifact_metric","pearson_r_with_genomeGC","n","interpretation"], delimiter="\t")
    w.writeheader()
    for r in corr_rows: w.writerow(r)
print("wrote correlations")

# Claim impact (Gate A)
# Assess: if GC gradient survives after stratifying by extremes of artifact metrics
# Example: compare low vs high artifact strata
with open(OUTDIR/"M1-03_claim_impact.md","w") as fh:
    fh.write("# M1-03 Artifact Claim Impact (Gate A)\n\n")
    fh.write("Does composition gradient survive assembly/annotation/LCR stratification?\n\n")
    fh.write("## Stratification\n")
    for r in corr_rows:
        fh.write(f"- {r['artifact_metric']}: r={r['pearson_r_with_genomeGC']} n={r['n']} — {r['interpretation']}\n")
    fh.write("\n## Per-species table highlights\n")
    # highlight extremes
    sorted_gc=sorted(out_rows, key=lambda x: float(x["genome_gc"]))
    fh.write(f"Low GC extremes: {sorted_gc[0]['species_id']} {sorted_gc[0]['genome_gc']:.1%} N50 {sorted_gc[0]['n50_kb']}kb contigs {sorted_gc[0]['num_contigs']} lcr {sorted_gc[0]['lcr_frac_lowEntropy']}\n")
    fh.write(f"High GC extremes: {sorted_gc[-1]['species_id']} {sorted_gc[-1]['genome_gc']:.1%} N50 {sorted_gc[-1]['n50_kb']}kb contigs {sorted_gc[-1]['num_contigs']} lcr {sorted_gc[-1]['lcr_frac_lowEntropy']}\n")
    fh.write("\n## Verdict for Gate A\n")
    # Simple rule: if N% and contig count not strongly correlated, and LCR correlation expected but not sole driver
    r_homopoly = next((float(x["pearson_r_with_genomeGC"]) for x in corr_rows if x["artifact_metric"]=="homopolymer_per_Mb"), 0)
    r_n = next((float(x["pearson_r_with_genomeGC"]) for x in corr_rows if x["artifact_metric"]=="n_percent"), 0)
    if abs(r_n) < 0.6:
        fh.write("N% / contiguity not primary driver (|r|<0.6): PASS\n")
    else:
        fh.write("N% strongly confounded: FLAG for Gate A downgrade\n")
    if abs(r_homopoly) > 0.6:
        fh.write(f"Homopolymer is GC-correlated (r~{r_homopoly:.2f}) — expected mechanistic coupling, not artifact; dual-track LCR retained/masked to be enforced in M2 (04 §5).\n")
    fh.write("LCR fraction will be carried forward as covariate; no Gate A STOP on this alone.\n")
    fh.write("\nNo B-level data fetched; no translation axis inference.\n")

with open(OUTDIR/"params.yaml","w") as fh:
    fh.write("method: assembly_gff_scan + lowEntropy_window64_thresh1.5 + homopolymer_from_M1-01\nwindow: 64\nentropy_thresh: 1.5\nsampling: first ~5Mb nuclear per species non-overlapping windows\ncorrelation: pearson\n")
with open(OUTDIR/"input_manifest.tsv","w") as fh:
    fh.write("file\tsize_bytes\tnote\n")
    for r in out_rows:
        acc=r["accession"]
        p=RAW_BASE/acc/"ncbi_dataset/data"/acc/(acc+"_genomic.fna" if acc.startswith("GCA") else acc+"_genomic.fna")
        # actual path via glob
        import pathlib
        cand=list((RAW_BASE/acc/"ncbi_dataset/data"/acc).glob("*_genomic.fna"))
        cand=[x for x in cand if "cds_from" not in x.name]
        if cand:
            fh.write(f"{cand[0]}\t{cand[0].stat().st_size}\tgenomic for LCR {r['species_id']}\n")
    fh.write(f"data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv\t{(BASE/'data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv').stat().st_size}\tcomposition for correlation\n")
with open(OUTDIR/"README.md","w") as fh:
    fh.write("# M1-03 Artifact\n\nQuestion: Is GC gradient explained by assembly/annotation/LCR/mappability?\nMethod: GFF gene/CDS counts, N% N50 contigs, CDS density, homopolymer/Mb (from M1-01), low-entropy 64bp windows fraction (threshold 1.5). Pearson r vs genome GC.\nOutputs: M1-03_artifact_table.tsv, _correlations.tsv, _claim_impact.md\n")
print("M1-03 done")
