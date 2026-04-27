import pytest

from mutlibrary.models.clinvar import (
    HGVSVariantGenerator,  # type: ignore[import]
)


@pytest.fixture
def hgvs_generator(tp53_201_cds, tp53_201_genomic):
    generator = HGVSVariantGenerator(
        genomic=tp53_201_genomic,
        cds=tp53_201_cds,
        # transcript_ac="NM_000546.6",
        # chrom_ac="NC_000017.11",
        cdot_json="/project/cdot-0.2.21.refseq.grch38_tp53.json",
        seq_id_prefix="var",
        genomic_flanks=(7673511, 7673636),  # checked in ensembl browser
    )
    return generator


@pytest.fixture
def variants(hgvs_generator, clinvar_cases):
    variants = {}
    for _, row in clinvar_cases.iterrows():
        hgvs_c = row["corrected_hgvs_c"]
        variant = hgvs_generator.from_hgvs(hgvs_c)
        variants[hgvs_c] = variant
    return variants


@pytest.fixture
def non_corrected_variants(hgvs_generator, clinvar_cases):
    variants = {}
    for _, row in clinvar_cases.iterrows():
        hgvs_c = row["hgvs_c"]
        variant = hgvs_generator.from_hgvs(hgvs_c)
        variants[hgvs_c] = variant
    return variants


def test_all_mutation_sequences(clinvar_cases, variants):
    for _, row in clinvar_cases.iterrows():
        hgvs_c = row["corrected_hgvs_c"]
        variant = variants[hgvs_c]
        expected_result = row["sequence_mutalyzer"]
        assert (
            variant.coding == expected_result
        ), f"Expected {expected_result} but got {variant.seq} for {variant}"
    assert len(variants) == len(
        clinvar_cases
    ), f"Expected {len(clinvar_cases)} variants but got {len(variants)}"


def test_all_mutation_sequences_with_non_canonical(
    clinvar_cases, non_corrected_variants
):
    print("reverse")
    for _, row in clinvar_cases.iterrows():
        print(reversed(row["sequence_mutalyzer"]))
        hgvs_c = row["hgvs_c"]
        variant = non_corrected_variants[hgvs_c]
        expected_result = row["sequence_mutalyzer"]
        assert (
            variant.coding == expected_result
        ), f"Expected {expected_result} but got {variant.coding} for {variant}"
        assert (
            variant.genomic == expected_result
        ), f"Expected {expected_result} but got {variant.coding} for {variant}"
    assert len(non_corrected_variants) == len(
        clinvar_cases
    ), f"Expected {len(clinvar_cases)} variants but got {len(non_corrected_variants)}"


def test_all_mutation_hgvs_g(clinvar_cases, variants):
    for _, row in clinvar_cases.iterrows():
        hgvs_c = row["corrected_hgvs_c"]
        variant = variants[hgvs_c]
        expected_result = row["corrected_hgvs_g"]
        assert (
            variant.hgvs_g == expected_result
        ), f"Expected {expected_result} but got {variant.hgvs_g} for {variant}"


def test_mutation_hgvs_p(clinvar_cases, variants):
    for _, row in clinvar_cases.iterrows():
        hgvs_c = row["corrected_hgvs_c"]
        variant = variants[hgvs_c]
        expected_result = row["predicted_hgvs_p"] or None
        assert (
            variant.hgvs_p == expected_result
        ), f"Expected {expected_result} but got {variant.hgvs_p} for {variant}"
        print(variant)


# def test_mutation_hgvs_r(clinvar_cases, variants):
#     for _, row in clinvar_cases.iterrows():
#         hgvs_c = row["corrected_hgvs_c"]
#         variant = variants[hgvs_c]
#         expected_result = row["predicted_hgvs_r"]
#         assert (
#             variant.hgvs_r == expected_result
#         ), f"Expected {expected_result} but got {variant.hgvs_r} for {hgvs_c}"


# def test_all_mutations(clinvar_cases, tp53_201_cds, tp53_201_genomic):
#     print(tp53_201_genomic, dir(tp53_201_genomic))
#     cds_start_in_genomic_0 = tp53_201_genomic.seq.find(tp53_201_cds.seq[:15])
#     generator = HGVSVariantGenerator(
#         genomic=tp53_201_genomic,
#         cds=tp53_201_cds,
#         # genomic_start=7668371,  # this is an 1-based chromosomal position, like in Ensembl
#         transcript_ac="NM_000546.6",
#         chrom_ac="NC_000017.11",
#         # chromosome="17",
#         cdot_json="/project/cdot-0.2.21.refseq.grch38.json.gz",
#         # assembly="GRCh38",
#         # strand=g_strand,
#         cds_start_in_genomic_0=cds_start_in_genomic_0,
#         flanking=10,
#         seq_id_prefix="var",
#         annotation=None,
#         genomic_flanks=(7673511, 7673636),
#     )
#     for _, row in clinvar_cases.iterrows():
#         hgvs_c = row["corrected_hgvs_c"]
#         expected_result = row["sequence_mutalyzer"]
#         variant = generator.from_hgvs_c(hgvs_c)
#         print(variant)
#         assert (
#             variant.seq == expected_result
#         ), f"Expected {expected_result} but got {variant.seq} for {hgvs_c}"
