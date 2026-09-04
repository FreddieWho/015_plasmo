#!/usr/bin/env python3
"""
M2-01 subfamily split check (post-Gate-B audit fix, 2026-09-03)
Motivation: primary null weights codons only by 3rd-position GC; for 6-fold families
(L/S/R) the 2-codon and 4-codon subfamilies differ at position 1 (TTR vs CTN etc.),
so position-1 composition effects leak into the "residual". Verified: CTN fraction of L
correlates with genome GC at r=0.990.
This script re-tests L/S/R with subfamily-aware nulls (TTR+CTN, AGY+TCN, AGR+CGN treated
as separate families), same per-gene GC3 rule, same QC filters, 200 pooled-multinomial
perms seed 42. Output: M2-01_subfamily_split_check.tsv + printed summary.
"""
import re, sys, math, collections
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/huyudi/015_plasmo")
OUTDIR = ROOT / "data/derived/WP2/M2-01_codon"
np.random.seed(42)
N_PERM = 200

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
STOPS = {'TAA','TAG','TGA'}
SUBFAMILIES = {
    'L_TTR': ['TTA','TTG'],
    'L_CTN': ['CTT','CTC','CTA','CTG'],
    'S_AGY': ['AGT','AGC'],
    'S_TCN': ['TCT','TCC','TCA','TCG'],
    'R_AGR': ['AGA','AGG'],
    'R_CGN': ['CGT','CGC','CGA','CGG'],
}

frozen = pd.read_csv(ROOT/"data/metadata/frozen_assembly_candidates.tsv", sep='\t')
acc_map = dict(zip(frozen['species_id'], frozen['frozen_accession']))
comp = pd.read_csv(ROOT/"data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv", sep='\t')
gc_by_sid = {r.species_id: r.genome_gc for _, r in comp.iterrows()}
SPECIES = ["SP001","SP002","SP003","SP004","SP005","SP006","SP007","SP008",
           "SP011","SP012","SP013","SP014","SP015","SP016","SP017","SP018"]

# original joint-family results for comparison
bg = pd.read_csv(OUTDIR/"M2-01_background_explained.tsv", sep='\t')
orig = bg[bg['aa'].isin(['L','S','R'])][['species_id','aa','prop_explained_GC3_masked','q_value_masked']]

def parse_cds_genes(path):
    """yield (locus, seq) representative = longest CDS per locus, with same QC filters as main run"""
    gene_dict = collections.defaultdict(list)
    header = None; chunks = []
    def flush():
        if header is None: return
        seq = ''.join(chunks).upper()
        m = re.search(r'\[locus_tag=([^\]]+)\]', header)
        locus = m.group(1) if m else header[:40]
        partial = 'partial' in header.lower()
        gene_dict[locus].append((seq, partial))
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            flush(); header = line; chunks = []
        else:
            chunks.append(line.strip())
    flush()
    for locus, seqs in gene_dict.items():
        seq, partial = max(seqs, key=lambda x: len(x[0]))
        L = len(seq)
        if L < 150 or L % 3 != 0 or partial or 'N' in seq: continue
        codons = [seq[i:i+3] for i in range(0, L-2, 3)]
        if any(c in STOPS for c in codons[:-1]): continue
        yield locus, codons

def gc3_of(codons):
    thirds = [c[2] for c in codons]
    return sum(1 for b in thirds if b in 'GC')/len(thirds)

rows = []
for sid in SPECIES:
    acc = acc_map[sid]
    path = ROOT/f"data/raw/ncbi-datasets/{acc}/ncbi_dataset/data/{acc}/cds_from_genomic.fna"
    obs = {sf: collections.Counter() for sf in SUBFAMILIES}
    exp = {sf: collections.Counter() for sf in SUBFAMILIES}
    for locus, codons in parse_cds_genes(path):
        gc3 = gc3_of(codons)
        cnt = collections.Counter(codons)
        for sf, cods in SUBFAMILIES.items():
            n = sum(cnt[c] for c in cods)
            if n == 0: continue
            w = {c: (gc3 if c[2] in 'GC' else 1-gc3) for c in cods}
            sw = sum(w.values())
            for c in cods:
                obs[sf][c] += cnt[c]
                exp[sf][c] += n * (w[c]/sw)
    for sf, cods in SUBFAMILIES.items():
        o = np.array([obs[sf][c] for c in cods], dtype=float)
        e = np.array([exp[sf][c] for c in cods], dtype=float)
        N = o.sum()
        if N < 20:
            continue
        chi2 = np.sum((o-e)**2/np.where(e>0, e, np.nan))
        eu = N/len(cods)
        chi2u = np.sum((o-eu)**2/eu)
        prop = 1 - chi2/chi2u if chi2u > 0 else np.nan
        # pooled-multinomial permutation p (same approximation as main run)
        pp = e/e.sum()
        sims = np.random.multinomial(int(N), pp, size=N_PERM)
        null = np.sum((sims - e)**2/np.where(e>0, e, 1), axis=1)
        p = (np.sum(null >= chi2) + 1)/(N_PERM+1)
        # GC-ending fraction residual (direction)
        gc_end_obs = sum(obs[sf][c] for c in cods if c[2] in 'GC')/N
        gc_end_exp = sum(exp[sf][c] for c in cods if c[2] in 'GC')/e.sum()
        rows.append({'species_id':sid, 'genome_gc':gc_by_sid[sid], 'family':sf[0], 'subfamily':sf,
                     'k':len(cods), 'N':int(N), 'prop_explained_GC3':prop,
                     'chi2':chi2, 'chi2_uniform':chi2u, 'perm_p':p,
                     'gcend_resid':gc_end_obs-gc_end_exp})

df = pd.DataFrame(rows)
df.to_csv(OUTDIR/"M2-01_subfamily_split_check.tsv", sep='\t', index=False)

print("=== Subfamily split null vs original joint null (prop_explained_GC3, unmasked track) ===\n")
for fam in ['L','S','R']:
    sub = df[df['family']==fam]
    o = orig[orig['aa']==fam]
    print(f"{fam}: joint-null median prop_explained = {o['prop_explained_GC3_masked'].median():.3f} (masked track), joint sig frac = {(o['q_value_masked']<0.05).mean():.2f}")
    for sf in sorted(sub['subfamily'].unique()):
        s = sub[sub['subfamily']==sf]
        gcresid_sign = np.sign(s['gcend_resid'])
        print(f"    {sf}: median prop_explained {s['prop_explained_GC3'].median():.3f}  "
              f"sig {int((s['perm_p']<0.05).sum())}/{len(s)} (perm p<0.05)  "
              f"gcend_resid +{int((gcresid_sign>0).sum())}/-{int((gcresid_sign<0).sum())}  "
              f"median |gcend_resid| {s['gcend_resid'].abs().median():.4f}")
    print()
print("Note: prop_explained here = 1 - chi2_GC3/chi2_uniform within subfamily, unmasked "
      "(LCR negligible: codon-level LCR ~1.6%). Compare vs joint null above and vs all-family median 0.845.")
