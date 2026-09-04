#!/usr/bin/env python3
"""
M1-01 Composition map: genome / compartment / CDS / GC1-3 / 4D / dinucleotide / homopolymer
Input: frozen_assembly_candidates.tsv + NCBI Datasets unpacked genomic.fna / cds_from_genomic.fna
Output: M1-01_composition_table.tsv, compartment, dinucleotide, homopolymer + summary/docs
Contract: 04_ANALYSIS_CONTRACT §3 (A-level only), §10
Perf: uses str.count / slices / regex — no per-base Python loops over 400 Mb
"""
import pathlib, json, csv, re, collections, sys
from collections import Counter, defaultdict

BASE = pathlib.Path("/home/huyudi/015_plasmo")
FROZEN = BASE/"data/metadata/frozen_assembly_candidates.tsv"
RAW_BASE = BASE/"data/raw/ncbi-datasets"
OUTDIR = BASE/"data/derived/WP1/M1-01_composition"
OUTDIR.mkdir(parents=True, exist_ok=True)

FOURFOLD_CODONS = {"GCT","GCC","GCA","GCG","CGT","CGC","CGA","CGG","GGT","GGC","GGA","GGG",
                    "CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","TCT","TCC","TCA","TCG",
                    "ACT","ACC","ACA","ACG","GTT","GTC","GTA","GTG"}

# single compiled regex for homopolymer >=5 (ATGC only)
HP_REGEX = re.compile(r"A{5,}|T{5,}|G{5,}|C{5,}")

def parse_fasta(path):
    header=None; parts=[]
    with open(path) as fh:
        for line in fh:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts).upper()
                header=line
                parts=[]
            else:
                parts.append(line.strip())
        if header is not None:
            yield header, "".join(parts).upper()

def gc_stats(seq):
    """return gc_count, atgc, n"""
    # seq already upper
    g = seq.count("G"); c = seq.count("C")
    a = seq.count("A"); t = seq.count("T")
    n = seq.count("N")
    atgc = a+t+g+c
    gc = g+c
    return gc, atgc, n

def classify_contig(header):
    h = header.lower()
    if "apicoplast" in h or "plastid" in h or "chromosome: api" in h or "chromosome: pltd" in h:
        return "apicoplast"
    if "mitochondrion" in h or ("chromosome: mt" in h) or ("chromosome: mito" in h):
        return "mitochondrion"
    # NCBI mitochondrion header sometimes is "mitochondrion" alone without chromosome tag
    # check the word mito/apico before falling to nuclear
    if "mitochondrion" in h:
        return "mitochondrion"
    # chromosome detection must handle both "chromosome: 1" and "chromosome 1, whole genome shotgun sequence"
    # "chromosome Unknown" are unplaced scaffolds — treat as nuclear_other
    if "chromosome" in h:
        if "unknown" in h:
            return "nuclear_other"
        return "nuclear_chromosome"
    if "contig:" in h or "scaffold:" in h or "region:" in h or "contig " in h or "scaffold " in h:
        return "nuclear_contig"
    # fallback — treat as nuclear_other (will be summed)
    return "nuclear_other"

# frozen
species=[]
with open(FROZEN) as fh:
    header_line=fh.readline().strip().split("\t")
    idx={k:i for i,k in enumerate(header_line)}
    for line in fh:
        if not line.strip(): continue
        parts=line.rstrip("\n").split("\t")
        species.append((
            parts[idx["species_id"]], parts[idx["species"]], parts[idx["taxid"]],
            parts[idx["frozen_accession"]], parts[idx["frozen_is_annotated"]],
            parts[idx["plasmodb_counterpart"]], parts[idx["risk_flag"]], parts[idx["freeze_action"]]
        ))

results=[]; compartment_rows=[]; dinuc_rows=[]; homopoly_rows=[]

for sid, sp, taxid, acc, annotated, plasmodb, risk, action in species:
    print(f"Processing {sid} {sp} {acc} annotated={annotated}", flush=True)
    unpack_dir = RAW_BASE/acc/"ncbi_dataset/data"/acc
    # find genomic fna — exclude cds_from_genomic.fna
    candidates = [p for p in unpack_dir.glob("*_genomic.fna") if "cds_from" not in p.name]
    if not candidates:
        candidates = [p for p in unpack_dir.glob("*.fna") if "cds_from" not in p.name]
    if not candidates:
        print(f"  WARN no genomic fna for {acc}", flush=True)
        continue
    genomic_fna = candidates[0]
    cds_path = unpack_dir/"cds_from_genomic.fna"
    has_cds = cds_path.exists() and cds_path.stat().st_size>0

    comp_by_type=defaultdict(lambda: {"length":0,"gc_count":0,"atgc":0,"n_count":0,"count":0})
    total_len=0; total_gc=0; total_atgc=0; total_n=0
    total_dinuc_counter=Counter(); total_dinuc=0
    total_hp_counter=Counter(); total_hp_runs=0
    num_contigs=0
    lengths=[]

    for header, seq in parse_fasta(genomic_fna):
        num_contigs+=1
        comp = classify_contig(header)
        L = len(seq)
        lengths.append(L)
        gc_cnt, atgc, n_cnt = gc_stats(seq)
        total_len+=L; total_atgc+=atgc; total_gc+=gc_cnt; total_n+=n_cnt
        ct = comp_by_type[comp]
        ct["length"]+=L; ct["gc_count"]+=gc_cnt; ct["atgc"]+=atgc; ct["n_count"]+=n_cnt; ct["count"]+=1
        # dinucleotide — count over upper seq without crossing contig boundaries
        # use str.count overlapping? need per dinuc; loop in Python still heavy for 25Mb*18 but we can vectorize via counting 16 dinucs using str.count with trick
        # Safer: use collections.Counter on sliding window via Python loop — but we already did gc via count, dinuc still O(N*16). Use fast approach: iterate once counting via translate?
        # We'll do Counter on seq without N using fast method: build pairs string and count
        # Simple fast: count each of 16 dinucs using str.count on overlapping via regex? Instead do single pass in Python but with numba-like speed we approximate with slicing counts:
        # Use: for din in DINS: count = sum(1 for ...) — still loop.
        # Instead use a single Python loop over seq is still slow (25M * 18 = 450M iterations). So we need faster.
        # Use: seq_no_n = seq (with N removed? but we should skip N dinucs). Simplest fast: count dinucs via str.count for each of 16 using overlapping search.
        # Overlapping count via: seq.count(din) counts non-overlapping. For homopolymer AAA, AA overlapping would be 2 but non-overlapping gives 1. So wrong.
        # So we need overlapping. We can use regex lookahead or manual fast.
        # Use Python's efficient find loop per dinuc (16 * number_of_occurrences) — still heavy.
        # Alternative: use array('u') and numpy if available; else use a single pass with translation table in Python's C via `collections.Counter(zip(...))` which still loops in Python.
        # Pragmatic: use `re.finditer` with lookahead? Let's benchmark one contig: use `for i in range(len(seq)-1)` is slow.
        # Better: use `import numpy` if available to vectorize.
        # We'll try numpy path, fallback to optimized Python with memoryview.
        # For now attempt numpy:
        try:
            import numpy as np
            # map A->0, C->1, G->2, T->3, N->4
            # need fast mapping: create translation via str.translate
            # Use numpy from buffer
            # Quick: encode seq as bytes and map via lookup table
            arr = np.frombuffer(seq.encode(), dtype=np.uint8)  # ASCII
            # map: A=0,C=1,G=2,T=3,N=4, other=4
            lut = np.full(256, 4, dtype=np.uint8)
            lut[ord('A')] = 0; lut[ord('C')] = 1; lut[ord('G')] = 2; lut[ord('T')] = 3
            coded = lut[arr]
            # mask where either base is N (4)
            valid = (coded[:-1] != 4) & (coded[1:] != 4)
            # dinuc index = first*4+second (0..15)
            din_idx = coded[:-1]*4 + coded[1:]
            din_idx = din_idx[valid]
            if len(din_idx):
                counts = np.bincount(din_idx, minlength=16)
                DINS = ["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"]
                # Actually mapping A0 C1 G2 T3 -> idx 0*4+0=0 AA etc. Need order matching above
                # Our order: A0,C1,G2,T3 => AA(0), AC(1), AG(2), AT(3), CA(4), CC(5), CG(6), CT(7), GA(8), GC(9), GG(10), GT(11), TA(12), TC(13), TG(14), TT(15)
                for j, din in enumerate(DINS):
                    c = int(counts[j])
                    if c:
                        total_dinuc_counter[din]+=c
                        total_dinuc+=c
            else:
                pass
        except ImportError:
            # fallback single-pass Python with minimal overhead using memoryview + manual
            # still single pass but in Python loop — only if numpy missing
            for i in range(len(seq)-1):
                a=seq[i]; b=seq[i+1]
                if a=='N' or b=='N' or a not in "ACGT" or b not in "ACGT":
                    continue
                total_dinuc_counter[a+b]+=1
                total_dinuc+=1
        # homopolymer via regex (C-level)
        # HP_REGEX finds runs >=5 efficiently in C
        for m in HP_REGEX.finditer(seq):
            run = m.group(0)
            base = run[0]
            total_hp_counter[base]+=1
            total_hp_counter[f"{base}_bases"]+=len(run)
            total_hp_runs+=1

    genome_gc = total_gc/total_atgc if total_atgc else 0
    # N50
    lengths_sorted=sorted(lengths, reverse=True)
    half=total_len/2
    cum=0; n50=0
    for ll in lengths_sorted:
        cum+=ll
        if cum>=half:
            n50=ll; break

    # CDS composition — number of CDS ~5k per species, total ~500M bases overall worst but still OK; use count/slice
    cds_gc = cds_gc1 = cds_gc2 = cds_gc3 = cds_gc4d = None
    cds_count=0; cds_total_len=0
    cds_atgc=0; cds_gc_cnt=0
    gc1_gc=gc1_atgc=gc2_gc=gc2_atgc=gc3_gc=gc3_atgc=0
    gc4d_gc=gc4d_atgc=gc4d_sites=0
    bad_cds=0; short_cds=0
    if has_cds:
        for header, seq in parse_fasta(cds_path):
            cds_count+=1
            L=len(seq)
            cds_total_len+=L
            if L%3!=0: bad_cds+=1
            if L<150: short_cds+=1
            # overall CDS GC via counts
            cds_gc_cnt += seq.count("G")+seq.count("C")
            cds_atgc += seq.count("A")+seq.count("T")+seq.count("G")+seq.count("C")
            # GC1/2/3 via slicing (C-level)
            # pos0 = seq[0::3], pos1 seq[1::3], pos2 seq[2::3]
            s1=seq[0::3]; gc1_gc+=s1.count("G")+s1.count("C"); gc1_atgc+=s1.count("A")+s1.count("T")+s1.count("G")+s1.count("C")
            s2=seq[1::3]; gc2_gc+=s2.count("G")+s2.count("C"); gc2_atgc+=s2.count("A")+s2.count("T")+s2.count("G")+s2.count("C")
            s3=seq[2::3]; gc3_gc+=s3.count("G")+s3.count("C"); gc3_atgc+=s3.count("A")+s3.count("T")+s3.count("G")+s3.count("C")
            # 4D
            # iterate codons — at most ~500k codons per species, still OK in Python
            for j in range(0, L-2, 3):
                codon=seq[j:j+3]
                if len(codon)!=3 or 'N' in codon: continue
                if codon in FOURFOLD_CODONS:
                    gc4d_sites+=1
                    third=codon[2]
                    if third in "ATGC":
                        gc4d_atgc+=1
                        if third in "GC": gc4d_gc+=1
        cds_gc = cds_gc_cnt/cds_atgc if cds_atgc else None
        cds_gc1 = gc1_gc/gc1_atgc if gc1_atgc else None
        cds_gc2 = gc2_gc/gc2_atgc if gc2_atgc else None
        cds_gc3 = gc3_gc/gc3_atgc if gc3_atgc else None
        cds_gc4d = gc4d_gc/gc4d_atgc if gc4d_atgc else None

    # homopolymer density
    hp_density = total_hp_runs/(total_len/1e6) if total_len else 0

    row={
        "species_id": sid, "species": sp, "taxid": taxid, "frozen_accession": acc,
        "annotated": annotated, "risk_flag": risk, "freeze_action": action, "plasmodb_counterpart": plasmodb,
        "total_length_bp": total_len, "num_contigs": num_contigs, "n50_bp": n50,
        "n_content_bp": total_n, "n_percent": total_n/total_len if total_len else 0,
        "genome_gc": genome_gc, "genome_gc_count": total_gc, "genome_atgc": total_atgc,
        "cds_count": cds_count, "cds_total_len": cds_total_len,
        "cds_gc": cds_gc, "gc1": cds_gc1, "gc2": cds_gc2, "gc3": cds_gc3, "gc4d": cds_gc4d,
        "gc4d_sites": gc4d_sites, "gc4d_gc": gc4d_gc, "gc4d_atgc": gc4d_atgc,
        "bad_cds_not_div3": bad_cds, "short_cds_lt150": short_cds,
        "homopolymer_runs_ge5": total_hp_runs, "homopolymer_density_per_Mb": hp_density,
        "dinuc_total": total_dinuc,
    }
    for comp in ["nuclear_chromosome","nuclear_contig","nuclear_other","mitochondrion","apicoplast"]:
        ct=comp_by_type.get(comp, {"length":0,"gc_count":0,"atgc":0,"n_count":0,"count":0})
        row[f"{comp}_len"]=ct["length"]
        row[f"{comp}_count"]=ct["count"]
        row[f"{comp}_gc"]=ct["gc_count"]/ct["atgc"] if ct["atgc"] else 0
        row[f"{comp}_n_bp"]=ct["n_count"]
    results.append(row)

    for comp, ct in comp_by_type.items():
        compartment_rows.append({"species_id":sid,"species":sp,"accession":acc,"compartment":comp,"count":ct["count"],"length_bp":ct["length"],"gc":ct["gc_count"]/ct["atgc"] if ct["atgc"] else 0,"n_bp":ct["n_count"]})
    for din,cnt in total_dinuc_counter.items():
        dinuc_rows.append({"species_id":sid,"species":sp,"dinucleotide":din,"count":cnt,"frequency":cnt/total_dinuc if total_dinuc else 0,"genome_gc":genome_gc})
    for base in ["A","T","G","C"]:
        homopoly_rows.append({"species_id":sid,"species":sp,"base":base,"runs_ge5":total_hp_counter.get(base,0),"bases_in_runs":total_hp_counter.get(f"{base}_bases",0),"genome_len":total_len,"density_per_Mb":total_hp_counter.get(base,0)/(total_len/1e6) if total_len else 0})

# write — guard when no results (e.g. missing data) so len(results)==0 won't crash previous indexing pattern
import csv
out_main = OUTDIR/"M1-01_composition_table.tsv"
if results:
    fieldnames=list(results[0].keys())
    with open(out_main,"w",newline="") as fh:
        w=csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in results: w.writerow(r)
    print(f"Wrote {out_main} {len(results)} rows")
else:
    print("No results — check inputs")

out_comp = OUTDIR/"M1-01_compartment_breakdown.tsv"
with open(out_comp,"w",newline="") as fh:
    w=csv.DictWriter(fh, fieldnames=["species_id","species","accession","compartment","count","length_bp","gc","n_bp"], delimiter="\t")
    w.writeheader()
    for r in compartment_rows: w.writerow(r)
print(f"Wrote {out_comp}")

out_dinuc = OUTDIR/"M1-01_dinucleotide.tsv"
with open(out_dinuc,"w",newline="") as fh:
    w=csv.DictWriter(fh, fieldnames=["species_id","species","dinucleotide","count","frequency","genome_gc"], delimiter="\t")
    w.writeheader()
    for r in dinuc_rows: w.writerow(r)
print(f"Wrote {out_dinuc}")

out_hp = OUTDIR/"M1-01_homopolymer.tsv"
with open(out_hp,"w",newline="") as fh:
    w=csv.DictWriter(fh, fieldnames=["species_id","species","base","runs_ge5","bases_in_runs","genome_len","density_per_Mb"], delimiter="\t")
    w.writeheader()
    for r in homopoly_rows: w.writerow(r)
print(f"Wrote {out_hp}")

# summary
with open(OUTDIR/"M1-01_summary.md","w") as fh:
    fh.write("# M1-01 Composition Summary\n\n")
    if results:
        fh.write(f"Generated from {len(results)} species, frozen accessions {', '.join([r['frozen_accession'] for r in results])}\n\n")
        fh.write("| species_id | species | accession | len Mb | n50 kb | genome GC% | CDS GC% | GC1 | GC2 | GC3 | GC4D | CDS n | mito len | api len | homopoly/Mb |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
        for r in sorted(results, key=lambda x: x["genome_gc"]):
            def pct(v): return f"{v*100:.1f}" if v is not None else "NA"
            fh.write(f"| {r['species_id']} | {r['species'].replace('Plasmodium','P.').replace('Babesia','B.').replace('Theileria','T.').replace('Toxoplasma','T.')} | {r['frozen_accession']} | {r['total_length_bp']/1e6:.2f} | {r['n50_bp']/1000:.0f} | {r['genome_gc']*100:.1f} | {pct(r['cds_gc'])} | {pct(r['gc1'])} | {pct(r['gc2'])} | {pct(r['gc3'])} | {pct(r['gc4d'])} | {r['cds_count']} | {r['mitochondrion_len']} | {r['apicoplast_len']} | {r['homopolymer_density_per_Mb']:.1f} |\n")
        fh.write("\nNotes: mito/api lengths =0 indicate NCBI assembly did not include organellar contigs (e.g. Pf 3D7 lacks mito/api in GCF_000002765.6). Method limitation, not data error; do not fetch organelle separately until Gate B (04 §3). CDS GC* = NA where no cds_from_genomic.fna.\n")
    else:
        fh.write("No results.\n")
print("Wrote summary.md")
