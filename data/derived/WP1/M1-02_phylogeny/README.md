# M1-02 Phylogeny + Ancestral GC + Independence

Scope: Gate A C1 evolution claim (02 §7). Tree is composition-independent (k-mer Jaccard, not GC).
Method: k=21 sketches 10k, blake2b min-hash, Jaccard distance, Bio.Phylo NJ rooted on Tg (SP018). Ancestral GC = mean of children (squared-change parsimony, equal branches). Leave-one-clade-out tests gradient robustness.
Outputs: M1-02_species_tree.nwk, _kmer_distance.tsv, _ancestral_gc.tsv, _leave_one_clade_out.tsv, _independence_grade.tsv

