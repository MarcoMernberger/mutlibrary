import itertools

import pandas as pd
import pytest
from Bio.Data import CodonTable  # type: ignore[import]
from Bio.SeqIO import SeqRecord  # type: ignore[import]
from pandas import DataFrame

from mutlibrary.models.mutator import MultiSeqMutatorHGVS, MutatorHGVS

########################################################################################
# Fixtures
########################################################################################


@pytest.fixture
def records(test_fasta):
    major = MutatorHGVS()
    return major.read_records(test_fasta)


@pytest.fixture
def regions(test_region, records):
    major = MutatorHGVS()
    return major.region_definitions(test_region, records)


@pytest.fixture
def annotations(test_annotation, records):
    major = MutatorHGVS()
    return major.record_annotations(test_annotation, records)


@pytest.fixture
def sample(records):
    return records["Test_Exon_Truncation"]


@pytest.fixture
def cds(records):
    return records["Test_CDS"]


@pytest.fixture
def mutator(sample, regions, cds, annotations):
    return MultiSeqMutatorHGVS(
        sample, regions[sample.id], cds, annotations[sample.id]
    )


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
    regions = major.region_definitions(test_region, records)
    assert isinstance(regions, dict)
    for record_name in regions:
        region_list = regions[record_name]
        assert isinstance(region_list, list)
        for item in region_list:
            assert "name" in item
            assert "type" in item
            assert "start" in item
            assert "end" in item
            assert "length" in item


def test_annotation(annotations, regions):
    assert isinstance(annotations, dict)
    for record_name in regions:
        assert (
            record_name in annotations
        ), f"Missing annotation for {record_name}"
    annotation_for_sequence = annotations["Test_Exon_Truncation"]
    assert isinstance(annotation_for_sequence, dict)
    assert "genomic_start" in annotation_for_sequence
    assert annotation_for_sequence["chromosome"] == "17"
    assert annotation_for_sequence["strand"] == "+"
    assert annotation_for_sequence["genomic_start"] == 10000
    assert annotation_for_sequence["assembly"] == "GRCh38"
    assert annotation_for_sequence["transcript_ensembl_id"] == "ENST00000269305"
    assert annotation_for_sequence["genomic_id"] == "NM_000546.6"


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


def test_insertions(mutations, annotations, sample, regions, mutator):
    genomic_start = annotations[sample.id]["genomic_start"]
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
        for l in [1, 2, 3]
        for p in itertools.product(bases, repeat=l)
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
        seq = row["seq"]

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


def test_deletions(mutations, annotations, sample):
    genomic_start = annotations[sample.id]["genomic_start"]

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
        seq = row["seq"]
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
        assert len(seq) == len(mutations["seq"].iloc[0]) - len(
            ref
        ), f"Length mismatch after deletion {ref}"

        # genomic mapping check
        assert (
            pos == row["genomic_pos"] - genomic_start
        ), f"Genomic position mismatch for {row}"


def test_snps(mutations, annotations, sample):
    genomic_start = annotations[sample.id]["genomic_start"]
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
        seq = row["seq"]

        assert ref != alt, f"Ref and alt identical in {row}"

        assert seq[pos] == alt, f"SNP not applied correctly at position {pos}"

        assert (
            pos == row["genomic_pos"] - genomic_start
        ), f"Genomic position mismatch for {row}"


def test_amino_acid_mutations(mutations, annotations, sample):
    genomic_start = annotations[sample.id]["genomic_start"]

    # ---- Filter AA mutations ----
    aa = mutations[mutations["mutation_type"].isin(["missense", "nonsense"])]
    # ---- Expected structure ----
    # (This assumes ONLY exon region is tested)
    expected_atg = {
        "ATA": "I",
        "AAG": "K",
        "AAC": "N",
        "AGC": "S",
        "ACG": "T",
        "AGG": "R",
    }
    expected_taa = {
        "GCA": "A",
        "GAA": "E",
        "GGA": "G",
        "ATA": "I",
        "AAA": "K",
        "TTA": "L",
        "CCA": "P",
        "CAA": "Q",
        "AGA": "R",
        "TCA": "S",
        "ACA": "T",
        "GTA": "V",
    }
    aa.to_csv("test_aa.tsv", sep="\t")

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

    for seq in aa["seq"].unique():
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
        # "DEL_fs",
        "INS_1",
        "INS_2",
        "INS_3",
        # "INS_fs",
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
    dedup = mutator.consolidate(mutations)
    assert dedup.shape[0] == 1471
