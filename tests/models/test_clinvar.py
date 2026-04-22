from mutlibrary.models.clinvar import HGVSVariantGenerator


def test_all_mutations(clinvar_cases, tp53_201_cds, tp53_201_genomic):
    print(tp53_201_genomic, dir(tp53_201_genomic))
    cds_start_in_genomic_0 = tp53_201_genomic.seq.find(tp53_201_cds.seq[:15])
    generator = HGVSVariantGenerator(
        genomic=tp53_201_genomic,
        cds=tp53_201_cds,
        # genomic_start=7668371,  # this is an 1-based chromosomal position, like in Ensembl
        transcript_ac="NM_000546.6",
        chrom_ac="NC_000017.11",
        # chromosome="17",
        cdot_json="/project/cdot-0.2.21.refseq.grch38.json.gz",
        # assembly="GRCh38",
        # strand=g_strand,
        cds_start_in_genomic_0=cds_start_in_genomic_0,
        flanking=10,
        seq_id_prefix="var",
        annotation=None,
        genomic_flanks=(7673510, 7673635),
    )
    for _, row in clinvar_cases.iterrows():
        hgvs_c = row["corrected_hgvs_c"]
        expected_result = row["sequence_mutalyzer"]
        variant = generator.from_hgvs_c(hgvs_c)
        print(variant)
        assert (
            variant.seq == expected_result
        ), f"Expected {expected_result} but got {variant.seq} for {hgvs_c}"
