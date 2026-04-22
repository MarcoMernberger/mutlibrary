# Standard Codon Table
CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

BASES = ["A", "C", "G", "T"]

# ---------------------------
# Helper functions
# ---------------------------


def translate_codon(codon):
    return CODON_TABLE.get(codon, "?")


def hamming_distance(a, b):
    return sum(x != y for x, y in zip(a, b))


# -----------------------------
# 3. SNP GENERATION
# -----------------------------

BASES = ["A", "C", "G", "T"]


def generate_snps(seq):
    for i, base in enumerate(seq):
        for b in BASES:
            if b != base:
                yield {
                    "type": "SNP",
                    "pos": i,
                    "ref": base,
                    "alt": b,
                    "seq": seq[:i] + b + seq[i + 1 :],
                }


# ---------------------------
# 2. Amino acid substitutions
# ---------------------------


def generate_codon_variants(seq, exon_positions):
    """
    exon_positions: list of (start, end) (0-based, end exclusive)
    """
    variants = []

    for start, end in exon_positions:
        exon_seq = seq[start:end]

        for i in range(0, len(exon_seq), 3):
            codon = exon_seq[i : i + 3]
            if len(codon) < 3:
                continue

            original_aa = translate_codon(codon)

            # iterate over all codons
            for new_codon in CODON_TABLE:
                new_aa = CODON_TABLE[new_codon]

                if new_codon == codon:
                    continue

                # minimal edit rule
                dist = hamming_distance(codon, new_codon)

                # only take minimal changes for each AA
                if new_aa != original_aa:
                    mutated_exon = exon_seq[:i] + new_codon + exon_seq[i + 3 :]

                    mutated_seq = seq[:start] + mutated_exon + seq[end:]

                    variants.append(
                        {
                            "type": "AA_substitution",
                            "codon_pos": start + i,
                            "ref_codon": codon,
                            "alt_codon": new_codon,
                            "ref_aa": original_aa,
                            "alt_aa": new_aa,
                            "edit_distance": dist,
                            "sequence": mutated_seq,
                        }
                    )

    # optional: filter minimal edit per AA
    best = {}
    for v in variants:
        key = (v["codon_pos"], v["alt_aa"])
        if key not in best or v["edit_distance"] < best[key]["edit_distance"]:
            best[key] = v

    return list(best.values())


# ---------------------------
# 3. ClinVar mutations
# ---------------------------


def apply_clinvar_variants(seq, clinvar_list):
    """
    clinvar_list example:
    [
        {"type": "SNP", "pos": 100, "alt": "A"},
        {"type": "DEL", "start": 200, "end": 203},
        {"type": "INS", "pos": 300, "seq": "ATG"},
        {"type": "MNV", "start": 400, "seq": "GGA"}
    ]
    """
    variants = []

    for var in clinvar_list:
        if var["type"] == "SNP":
            mutated = seq[: var["pos"]] + var["alt"] + seq[var["pos"] + 1 :]

        elif var["type"] == "DEL":
            mutated = seq[: var["start"]] + seq[var["end"] :]

        elif var["type"] == "INS":
            mutated = seq[: var["pos"]] + var["seq"] + seq[var["pos"] :]

        elif var["type"] == "MNV":
            mutated = (
                seq[: var["start"]]
                + var["seq"]
                + seq[var["start"] + len(var["seq"]) :]
            )

        else:
            continue

        variants.append(
            {"type": "ClinVar", "variant": var, "sequence": mutated}
        )

    return variants


# ---------------------------
# MAIN PIPELINE
# ---------------------------


def generate_all_variants(seq, exon_positions, clinvar_list):
    snps = generate_snps(seq)
    aa_vars = generate_codon_variants(seq, exon_positions)
    clinvars = apply_clinvar_variants(seq, clinvar_list)

    return {"snps": snps, "aa_substitutions": aa_vars, "clinvar": clinvars}
