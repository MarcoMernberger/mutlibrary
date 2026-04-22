def hamming_distance(codon, syn_codon):
    return sum(1 for a, b in zip(codon, syn_codon) if a != b)
