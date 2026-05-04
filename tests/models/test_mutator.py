import itertools

import pandas as pd  # type: ignore[import]
import pytest  # type: ignore[import]
from Bio.Data import CodonTable  # type: ignore[import]
from Bio.SeqIO import SeqRecord  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

from mutlibrary.models.annotator import (  # type: ignore[import]
    # Mutalyzer,
    HGVSMutationAnnotator,
    MutationAnnotator,
)
from mutlibrary.models.mutator import (  # type: ignore[import]
    MultiSeqMutatorHGVS,
    MutatorHGVS,
)


class MockAnnotator(MutationAnnotator):
    def annotate(self, mutation):
        return mutation


########################################################################################
# Fixtures
########################################################################################


@pytest.fixture
def mutator(sample, regions, cds, genomic, annotations):
    return MultiSeqMutatorHGVS(
        sample,
        genomic,
        cds,
        regions[sample.id],
        annotations,
        annotator=MockAnnotator(),
    )


@pytest.fixture
def real_mutator(
    real_sample, real_regions, real_cds, real_genomic, real_annotations
):
    return MultiSeqMutatorHGVS(
        real_sample,
        real_genomic,
        real_cds,
        real_regions[real_sample.id],
        real_annotations,
        annotator=HGVSMutationAnnotator(),
    )


@pytest.fixture
def real_mutations(real_mutator):
    return real_mutator.collect_mutations()


@pytest.fixture
def mutations(mutator):
    return mutator.generate_mutations()


pd.set_option("display.max_columns", None)


def get_mutable_length(sample, regions):
    mutable_len = 0
    for region in regions[sample.id]:
        if "INS" in region["mutations"]:
            mutable_len += region["length"]
    return mutable_len


########################################################################################
# Test cases
########################################################################################
def test_init():

    mutator = MutatorHGVS()
    assert mutator is not None


def test_records(test_fasta):
    major = MutatorHGVS()
    records = major.read_records(test_fasta)
    assert isinstance(records, dict)
    for item in records:
        assert isinstance(records[item], SeqRecord)


def test_regions(test_fasta, test_region):
    major = MutatorHGVS()
    records = major.read_records(test_fasta)
    regions = major.region_definitions(
        test_region,
        records,
        ["TP53-201_cds_protein_coding", "17", "Test_CDS", "Test_genomic"],
    )
    assert isinstance(regions, dict)
    for record_name in regions:
        region_list = regions[record_name]
        assert isinstance(
            region_list, list
        ), f"Expected list of regions for {record_name}, was {type(region_list)}"
        for item in region_list:
            assert "name" in item
            assert "type" in item
            assert "start" in item
            assert "end" in item
            assert "length" in item
            assert "mutations" in item


def test_annotation(real_annotations, real_regions):
    assert isinstance(real_annotations, dict)
    for record_name in real_regions:
        assert (
            record_name in real_annotations
        ), f"Missing annotation for {record_name}"
    annotation_for_gene = real_annotations["TP53-201_cds_protein_coding"]
    annotation_for_cds = real_annotations["TP53-201_cds_protein_coding"]
    annotations_for_sample = real_annotations["TP53-201_Ex9_realtest"]
    assert isinstance(annotation_for_cds, dict)
    assert isinstance(annotations_for_sample, dict)
    assert isinstance(annotation_for_gene, dict)
    assert "transcript_id" in annotation_for_cds
    assert "transcript_ac" in annotation_for_cds
    assert "chromosome_ac" in annotation_for_gene
    assert "protein_ac" in annotation_for_cds
    assert annotation_for_cds["transcript_id"] == "ENST00000269305"
    assert annotation_for_cds["transcript_ac"] == "NM_000546.6"
    assert annotation_for_gene["chromosome_ac"] == "NC_000017.11"
    assert annotation_for_cds["protein_ac"] == "NP_000537.3"


def test_hgvs_on_generated_sequences(
    real_mutations, real_mutator, p53_examples
):
    assert isinstance(
        real_mutator.annotator, HGVSMutationAnnotator
    ), f"annotator should be HGVSMutationAnnotator, was {type(real_mutator.annotator)}"
    unchecked = set(p53_examples.index.values)
    for mutation in real_mutations:
        assert mutation.hgvs_g, f"Missing HGVS g for {mutation}"
        if mutation.region_type == "exon":
            assert mutation.hgvs_c, f"Missing HGVS c for {mutation}"
        assert mutation.coding, f"Missing coding sequence for {mutation}"
        assert mutation.genomic, f"Missing genomic sequence for {mutation}"
        if "TAGCACTGCCCAACCTAAACACCAGCTCCTCT" in mutation.coding:
            print("Found expected sequence in coding:", mutation)
        if "TAGCACTGCCCAACGTTAACACCAGCTCCTCT" in mutation.coding:
            print("Found expected mutated sequence in coding:", mutation)
        if "CTAGCACTGCCCAACAAACACCAGCTCCTCTC" in mutation.genomic:
            print("Found expected sequence in genomic:", mutation)
        if mutation.hgvs_g in unchecked:
            print(unchecked)
            unchecked.remove(mutation.hgvs_g)
            expected_hgvs_c = p53_examples.loc[mutation.hgvs_g, "hgvsc"]
            expected_sequence = p53_examples.loc[mutation.hgvs_g, "sequence"]
            assert (
                mutation.hgvs_c == expected_hgvs_c
            ), f"Expected HGVS c {expected_hgvs_c} but got {mutation.hgvs_c} for {mutation}"
            assert (
                expected_sequence in mutation.coding
            ), f"Expected sequence {expected_sequence} not found in coding {mutation.coding} for {mutation}"
    assert not unchecked, f"HGVS g not tested for: {unchecked}"


def test_generate_mutations(mutations):

    assert isinstance(mutations, DataFrame)
    assert len(mutations) > 0


def test_all_mutations(mutations):

    assert isinstance(mutations, DataFrame)
    for idx, row in mutations.iterrows():
        assert "mutation_type" in row
        assert "hgvs_c" in row
        assert "hgvs_p" in row


def test_write_mutations(tmp_path, mutator, sample):

    output_file = tmp_path / f"{sample.id}_variants.tsv"
    mutator.write_mutations(output_file)

    assert (
        output_file.exists()
    ), "write_mutations must create the output TSV file"
    content = output_file.read_text(encoding="utf-8")
    assert "hgvs_c" in content
    assert "mutation_type" in content
    lines = [line for line in content.splitlines() if line.strip()]
    assert len(lines) > 1, "Expected header plus at least one mutation row"


def test_insertions(mutations, sample, regions, mutator):
    genomic_start = (
        mutator.genomic_start
    )  # annotations[sample.id]["genomic_start"]
    bases = ["A", "C", "G", "T"]

    # ---- Filter insertions ----
    insertions = mutations[mutations["mutation_type"].str.startswith("INS")]

    # ---- Expected counts ----
    mutable_len = get_mutable_length(sample, regions)
    expected_positions = mutable_len + 1
    expected_per_pos = sum(4**k for k in range(1, mutator.max_len + 1))
    expected_total = expected_positions * expected_per_pos

    assert (
        len(insertions) == expected_total
    ), f"Expected {expected_total}, got {len(insertions)}"

    # ---- All insertion sequences exist ----
    expected_insertions = {
        "".join(p)
        for ll in [1, 2, 3]
        for p in itertools.product(bases, repeat=ll)
    }

    observed_insertions = set(insertions["alt"].unique())

    missing = expected_insertions - observed_insertions
    assert not missing, f"Missing insertions: {missing}"

    # ---- Frameshift correctness ----
    for _, row in insertions.iterrows():
        inlength = len(row["alt"])
        # if inlength % 3 != 0 and row["region_type"] == "exon":
        #     assert "INS_fs" in row["mutation_type"], f"Expected INS_fs in {row}"
        # else:
        assert (
            f"INS_{inlength}" in row["mutation_type"]
        ), f"Expected INS_{inlength} in {row}"

    # ---- Sequence correctness ----
    for _, row in insertions.iterrows():
        pos = row["mutation_pos"]
        alt = row["alt"]
        seq = row["coding"]

        assert (
            seq[pos : pos + len(alt)] == alt
        ), f"Insertion {alt} not at correct position {pos}"
        assert (
            pos == row["genomic_pos"] - genomic_start
        ), f"Genomic position mismatch for {row}"


# def test_insertions_hgvs(mutations, annotations, sample):

#     # ---- HGVS checks ----
#     hgvs_c_pattern = re.compile(
#         r".+:(c|g)\.\d+([ACGT]+>[ACGT]+|del|ins[ACGT]+)"
#     )
#     hgvs_p_pattern = re.compile(r".*:p\.[A-Z\*]\d+[A-Z\*]?")

#     for _, row in mutations.iterrows():

#         hgvs_c = row["hgvs_c"]
#         hgvs_p = row["hgvs_p"]

#         # ---- c. vs g. logic ----
#         if row["region_type"] == "exon" and row["cds_pos"] is not None:
#             assert ":c." in hgvs_c, f"Expected c. HGVS for exon: {hgvs_c}"
#         else:
#             assert ":g." in hgvs_c, f"Expected g. HGVS for intron: {hgvs_c}"

#         # ---- basic syntax ----
#         assert hgvs_c_pattern.match(hgvs_c), f"Invalid HGVS c: {hgvs_c}"

#         # ---- protein HGVS only for exon ----
#         if row["region_type"] == "exon" and row["ref_aa"] and row["alt_aa"]:
#             assert hgvs_p, "Missing protein HGVS"
#             assert hgvs_p_pattern.match(hgvs_p), f"Invalid HGVS p: {hgvs_p}"

#         else:
#             assert hgvs_p == "" or pd.isna(hgvs_p)

#     # ---- No mutations in constant regions ----
#     constant_regions = mutations[
#         mutations["region_type"].str.contains("const", na=False)
#     ]

#     assert len(constant_regions) == 0, "Mutations found in constant regions"


def test_deletions(mutations, sample, mutator):
    genomic_start = mutator.genomic_start

    # ---- Filter deletions ----
    deletions = mutations[mutations["mutation_type"].str.startswith("DEL")]

    # ---- Expected counts ----
    mutable_len = 20  # same assumption as insertions test

    expected_positions = mutable_len  # deletions are length-based
    expected_per_pos = 3  # 1mer, 2mer, 3mer deletions
    not_deleted_at_end = (
        2 + 1
    )  # last 2 positions cannot have 3mer deletions and the last not 2mer

    expected_total = expected_positions * expected_per_pos - not_deleted_at_end

    assert (
        len(deletions) == expected_total
    ), f"Expected {expected_total}, got {len(deletions)}"

    # ---- All deletion sequences exist ----
    deletable = sample.seq[5:-5]
    expected_deletions = []
    for pos in range(len(deletable)):
        for length in range(1, 4):
            if pos + length > len(deletable):
                continue
            ref_seq = deletable[pos : pos + length]
            expected_deletions.append(str(ref_seq))

    observed_deletions = set(deletions["ref"].unique())

    missing = set(expected_deletions) - observed_deletions
    assert not missing, f"Missing deletions: {missing}"

    # ---- Frameshift correctness ----
    for _, row in deletions.iterrows():
        dellen = len(row["ref"])

        # if dellen % 3 != 0 and row["region_type"] == "exon":
        #     assert "DEL_fs" in row["mutation_type"], f"Expected DEL_fs in {row}"
        # else:
        assert (
            f"DEL_{dellen}" in row["mutation_type"]
        ), f"Expected DEL_{dellen} in {row}"

    # ---- Sequence correctness ----
    for _, row in deletions.iterrows():
        pos = row["mutation_pos"]
        ref = row["ref"]
        seq = row["coding"]
        # deleted segment must NOT be present anymore
        seq_stripped = sample.seq[:pos] + sample.seq[pos + len(ref) :]
        deleted = sample.seq[pos : pos + len(ref)]

        assert (
            seq == seq_stripped
        ), f"Mutated sequence should be {deleted}, was {seq}"
        assert (
            deleted == ref
        ), f"Deleted sequence should be {ref}, was {deleted}"
        # ensure sequence length consistency
        assert len(seq) == len(mutations["coding"].iloc[0]) - len(
            ref
        ), f"Length mismatch after deletion {ref}"

        # genomic mapping check
        assert (
            pos == row["genomic_pos"] - genomic_start
        ), f"Genomic position mismatch for {row}"


def test_snps(mutations, mutator):
    genomic_start = mutator.genomic_start
    bases = ["A", "C", "G", "T"]

    # ---- Filter SNPs ----
    snps = mutations[mutations["mutation_type"].str.startswith("SNP")]

    # ---- Expected counts ----
    mutable_len = 20
    expected_positions = mutable_len
    expected_per_pos = 3  # 3 alternatives per base
    expected_total = expected_positions * expected_per_pos

    assert (
        len(snps) == expected_total
    ), f"Expected {expected_total}, got {len(snps)}"

    # ---- All SNP alternatives exist ----
    expected_snps = set(bases)

    observed_snps = set(snps["alt"].unique())

    missing = expected_snps - observed_snps
    assert not missing, f"Missing SNP bases: {missing}"

    # ---- Mutation correctness ----
    for _, row in snps.iterrows():
        pos = row["mutation_pos"]
        ref = row["ref"]
        alt = row["alt"]
        seq = row["coding"]

        assert ref != alt, f"Ref and alt identical in {row}"

        assert seq[pos] == alt, f"SNP not applied correctly at position {pos}"

        assert (
            pos == row["genomic_pos"] - genomic_start
        ), f"Genomic position mismatch for {row}"


def test_amino_acid_mutations(mutations, mutator, sample):
    genomic_start = mutator.genomic_start

    # ---- Filter AA mutations ----
    aa = mutations[mutations["mutation_type"].isin(["missense", "nonsense"])]
    # ---- Expected structure ----
    # (This assumes ONLY exon region is tested)
    # expected_atg = {
    #     "ATA": "I",
    #     "AAG": "K",
    #     "AAC": "N",
    #     "AGC": "S",
    #     "ACG": "T",
    #     "AGG": "R",
    # }
    # expected_taa = {
    #     "GCA": "A",
    #     "GAA": "E",
    #     "GGA": "G",
    #     "ATA": "I",
    #     "AAA": "K",
    #     "TTA": "L",
    #     "CCA": "P",
    #     "CAA": "Q",
    #     "AGA": "R",
    #     "TCA": "S",
    #     "ACA": "T",
    #     "GTA": "V",
    # }

    mutable_full_codons = 2 * 20
    truncated_start = 6
    truncated_end = 12
    expected_total = mutable_full_codons + truncated_start + truncated_end

    assert len(aa) > 0, "No AA mutations generated"
    assert (
        len(aa) == expected_total
    ), f"Expected {expected_total} AA mutations, got {len(aa)}"

    # ---- Check mutation types ----
    for _, row in aa.iterrows():
        ref_aa = row["ref_aa"]
        alt_aa = row["alt_aa"]

        assert ref_aa is not None
        assert alt_aa is not None

        if alt_aa == "*":
            assert (
                row["mutation_type"] == "nonsense"
            ), f"Expected nonsense mutation: {row}"
        else:
            assert (
                row["mutation_type"] == "missense"
            ), f"Expected missense mutation: {row}"

    # ---- Codon consistency ----
    for _, row in aa.iterrows():
        ref_codon = row["ref"]
        alt_codon = row["alt"]

        assert len(ref_codon) == len(alt_codon)

        # ensure codon actually changes protein
        assert row["ref_aa"] != row["alt_aa"], f"No AA change detected: {row}"

    # --- group by codon position and check expected AA changes ---
    assert "mutation_pos" in aa.columns, "mutation_pos column missing"
    assert (
        aa["mutation_pos"].nunique() == 4
    ), "Expected exactly 4 reference codons"
    for codon_pos, group in aa.groupby("mutation_pos"):
        if codon_pos in [2, 5]:
            all_amino_acids = set(
                CodonTable.unambiguous_dna_by_name[
                    "Standard"
                ].forward_table.values()
            )
            all_amino_acids.remove(group["ref_aa"].values[0])

            for amino_acid in all_amino_acids:
                assert (
                    amino_acid in group["alt_aa"].unique()
                ), f"{amino_acid} codon expected at pos {codon_pos}"

    for seq in aa["coding"].unique():
        assert (
            seq[:10] == sample.seq[:10]
        ), f"First 10 bases should be unchanged in {seq}"
        assert (
            seq[-10:] == sample.seq[-10:]
        ), f"Last 10 bases should be unchanged in {seq}"

    # ---- Genomic consistency ----
    for _, row in aa.iterrows():
        pos = row["mutation_pos"]

        assert (
            pos == row["genomic_pos"] - genomic_start
        ), f"Genomic position mismatch for {row}"


def test_all_mutations_semantic(mutations):
    mutation_types = mutations["mutation_type"].unique()
    expected_types = {
        "SNP",
        "DEL_1",
        "DEL_2",
        "DEL_3",
        "INS_1",
        "INS_2",
        "INS_3",
        "missense",
        "nonsense",
    }
    missing_types = expected_types - set(mutation_types)
    assert not missing_types, f"Missing mutation types: {missing_types}"
    for idx, row in mutations.iterrows():
        assert "mutation_type" in row
        assert "hgvs_c" in row
        assert "hgvs_p" in row


def test_deduplication(mutations, mutator):
    mutations.to_csv("test.tsv", sep="\t")
    assert mutations.shape[0] == 1939
    dedup = mutator.deduplicate(mutations)
    assert dedup.shape[0] == 1471
