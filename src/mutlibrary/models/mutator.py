import itertools
import sys
from itertools import product
from pathlib import Path
from typing import Any, Literal

import pandas as pd  # type: ignore[import]
from Bio import Align, SeqIO  # type: ignore[import]
from Bio.Data import CodonTable  # type: ignore[import]
from mmalignments.models.elements import (  # type: ignore[import]
    Element,  # type: ignore[import]
    FileElement,  # type: ignore[import]
    # generate_element_key_name,
    element,  # type: ignore[import]
)
from mmalignments.models.parameters import Params  # type: ignore[import]
from mmalignments.models.tags import (  # type: ignore[import]
    ElementTag,  # type: ignore[import]
    # ExternalRunConfig,
    PartialElementTag,  # type: ignore[import]
)
from mmalignments.services.io import from_json  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

from mutlibrary.models.annotator import (  # type: ignore[import]
    HGVSMutationAnnotator,
    ManualMutationAnnotator,
    MaveMutationAnnotator,
    Mutation,
    MutationAnnotator,
)
from mutlibrary.models.clinvar import (
    HGVSVariantGenerator,  # type: ignore[import]
)

# HGVS standard requires three-letter amino acid codes for p. notation
_AA_1TO3: dict[str, str] = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
    "*": "Ter",
}

sys.path.append("mutlibrary/src")


# from mutlibrary.utils.distances import hamming_distance
def hamming_distance(codon, syn_codon):
    return sum(1 for a, b in zip(codon, syn_codon) if a != b)


class MultiSeqMutatorHGVS:
    def __init__(
        self,
        record: SeqIO.SeqRecord,
        regions: list[dict[str, Any]],
        genomic: SeqIO.SeqRecord | None = None,
        cds: SeqIO.SeqRecord | None = None,
        annotation: dict | None = None,
        annotator: MutationAnnotator | None = None,
    ):
        self.record = record
        self.regions = regions
        self.cds = cds
        self.genomic = genomic
        self.annotation = annotation
        self.annotator = annotator or ManualMutationAnnotator()
        self.genomic_start = (
            self.annotation["genomic_start"] if self.annotation else 0
        )
        self.genomic_id = annotation.get("genomic_id", "") if annotation else ""
        self.chromosome = annotation.get("chromosome", "") if annotation else ""
        self.__infer_from_genomic_description()
        self.bases = ["A", "C", "G", "T"]
        self.codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
        self.codon_map = dict(self.codon_table.forward_table)
        for stop in self.codon_table.stop_codons:
            self.codon_map[stop] = "*"
        self.combined_table = Path(f"cache/{self.record.id}_variants.tsv")
        self.cds_codons = self.build_cds_codon_mapping()
        self.max_len = 3

    def __infer_from_genomic_description(self):
        if self.genomic:
            desc = self.genomic.description
            splits = desc.split(":")
            strand = int(splits[-1])
            self.strand = int(strand)

    def map_exon_to_cds(self, exon_seq: str) -> int:
        if self.cds is None:
            raise ValueError(
                "CDS sequence is required for amino acid mutation generation"
            )
        cds_seq = str(self.cds.seq)
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -0.5
        best = next(iter(aligner.align(cds_seq, exon_seq)))
        return best.aligned[0][0][0]  # start of first aligned block in CDS

    def build_cds_codon_mapping(self):
        """
        Returns:
            mapping: list[int]  # cds_index -> genomic position
        """
        if self.cds is None:
            raise ValueError(
                "CDS sequence is required for amino acid mutation generation"
            )
        mapping = {}
        for index in range(0, len(self.cds.seq), 3):
            codon = str(self.cds.seq[index : index + 3])
            mapping[index] = codon
        return mapping

    def get_alternate_codons(
        self, codon: str, pos_in_codon: int, constant_front=True
    ) -> list[str]:
        mutations = []
        if constant_front:
            constant = codon[:pos_in_codon]
            variable = [
                "".join(p) for p in product(self.bases, repeat=3 - pos_in_codon)
            ]
            mutations = [constant + var for var in variable]

        else:
            constant = codon[pos_in_codon + 1 :]
            variable = [
                "".join(p) for p in product(self.bases, repeat=pos_in_codon + 1)
            ]
            mutations = [var + constant for var in variable]
        mutations = [m for m in mutations if m != codon]
        return mutations

    def get_minmimal_hamming_synonyms(
        self, codon: str, pos_in_codon: int, constant_front=True
    ) -> dict[str, tuple[str, int]]:
        mutations_with_synonyms = self.get_alternate_codons(
            codon, pos_in_codon, constant_front
        )

        minimal_hamming = {}
        for syn_codon in mutations_with_synonyms:
            syn = self.codon_map[syn_codon]
            hamming = hamming_distance(codon, syn_codon)
            if hamming == 0:
                continue
            if syn not in minimal_hamming or hamming < minimal_hamming[syn][1]:
                minimal_hamming[syn] = (syn_codon, hamming)
        return minimal_hamming

    def extract_mutable_exon_codons(self, region):
        """
        Returns list of tuples:
        (pos_in_exon, exon_subseq, codon, pos_in_codon, constant_front)
        """
        seq = str(self.record.seq)
        mapping = self.cds_codons
        mutatable_codon_positions = []
        start, end = region["start"], region["end"]
        exon = seq[start:end]
        cds_start = self.map_exon_to_cds(exon)
        # CDS Positionen, die in diesem Exon liegen
        cds_positions = [cds_start + i for i in range(len(exon))]
        if cds_start not in mapping:
            # handle truncated start codon
            for ii in [1, 2]:
                codon_start = cds_start - ii
                if codon_start in mapping:
                    codon = mapping[codon_start]
                    in_exon = exon[0 : 3 - ii]
                    # constant = codon[:ii]
                    mutatable_codon_positions.append(
                        (
                            codon_start,
                            0,  # position in exon
                            in_exon,
                            codon,
                            ii,  # position in codon if truncated
                            True,  # keep the front constant if truncated at the start
                        )
                    )
                    break

        last_covered_cds_pos = cds_positions[-1] + 1
        for exon_pos, cds_pos in enumerate(cds_positions):
            if cds_pos in mapping:  # this covers all full codons
                codon = mapping[cds_pos]
                cds_pos_end = cds_pos + 3
                if cds_pos_end > last_covered_cds_pos:
                    pos_in_codon = cds_pos_end - last_covered_cds_pos
                    in_exon = exon[exon_pos : cds_pos_end - pos_in_codon]
                    mutatable_codon_positions.append(
                        (
                            cds_pos,
                            exon_pos,  # position in exon
                            in_exon,  # truncated codon in exon
                            codon,  # full codon in CDS
                            pos_in_codon,  # position in codon if truncated
                            False,  # keep the end constant
                        )
                    )
                else:
                    # full codon in exo
                    in_exon = exon[exon_pos:cds_pos_end]
                    mutatable_codon_positions.append(
                        (
                            cds_pos,
                            exon_pos,  # position in exon
                            in_exon,  # codon in exon
                            codon,  # codon in frame, should be the same as above
                            0,  # position in codon if truncated, this is != 0
                            True,  # keep the front constant if truncated
                        )
                    )

        return mutatable_codon_positions

    def is_valid_window(self, pos, length, position_map):
        for i in range(pos, pos + length):
            if i not in position_map:
                return False
        return True

    ####################################################################################
    # The sequence generators
    ####################################################################################

    def generate_aa(self) -> list[Mutation]:
        # Amino Acid changes (exons only)

        variants = []
        seq = str(self.record.seq)
        seq_id = self.record.id

        for region in self.regions:
            if region["type"] != "exon" or "AA" not in region.get(
                "mutations", []
            ):
                continue
            if self.cds is None:
                raise ValueError(
                    "CDS sequence is required for amino acid mutation generation"
                )

            codon_blocks = self.extract_mutable_exon_codons(region)
            for (
                cds_pos,
                exon_pos,
                exon_subseq,
                codon,
                pos_in_codon,
                constant_front,
            ) in codon_blocks:

                current_start = region["start"] + exon_pos
                ref_aa = self.codon_map[codon]
                minimal_changes_codons = self.get_minmimal_hamming_synonyms(
                    codon, pos_in_codon, constant_front=constant_front
                )

                for alt_aa, (alt_codon, dist) in minimal_changes_codons.items():

                    if alt_aa == ref_aa:
                        continue

                    if pos_in_codon == 0:  # full_codon
                        mutated_seq = str(
                            seq[:current_start]
                            + alt_codon
                            + seq[current_start + 3 :]
                        )
                        alt = alt_codon
                        ref = codon
                    elif constant_front:
                        truncated_exon_length = len(exon_subseq)
                        alt = alt_codon[3 - truncated_exon_length :]
                        ref = codon[3 - truncated_exon_length :]
                        mutated_seq = str(
                            seq[:current_start]
                            + alt_codon[3 - truncated_exon_length :]
                            + seq[current_start + truncated_exon_length :]
                        )
                    elif not constant_front:
                        truncated_exon_length = len(exon_subseq)
                        alt = alt_codon[:truncated_exon_length]
                        ref = codon[:truncated_exon_length]
                        mutated_seq = str(
                            seq[:current_start]
                            + alt_codon[:truncated_exon_length]
                            + seq[current_start + truncated_exon_length :]
                        )
                    else:
                        raise ValueError("This should never happen")

                    # mutation type
                    mutation_type = "nonsense" if alt_aa == "*" else "missense"

                    mutation = Mutation(
                        seq_id=seq_id,
                        genomic_id=self.genomic_id,
                        chromosome=self.chromosome,
                        strand=self.strand,
                        mutation_pos=exon_pos,
                        genomic_pos=self.genomic_start + exon_pos,
                        region_type="exon",
                        ref=ref,
                        alt=alt,
                        mutation_type=mutation_type,
                        seq=mutated_seq,
                        cds_pos=cds_pos,
                        ref_codon=codon,
                        alt_codon=alt_codon,
                        ref_aa=ref_aa,
                        alt_aa=alt_aa,
                    )
                    mutation = self.annotator.annotate(mutation)
                    variants.append(mutation)

        return variants

    ####################################################################################
    # mutation generation
    ####################################################################################

    # -----------------------------
    # SNPs
    # -----------------------------
    def generate_snp(self) -> list[Mutation]:
        """
        Generates a list of all single SNP mutations for the record
        sequence. Each specified region in regions can have different mutation types
        specified, but the complete sequence is returned.
        The resulting mutations are annotated with HGVS nomenclature.

        Returns
        -------
        list[Mutation]
            A list of Mutation objects representing all generated SNP mutations.
        """
        variants = []
        seq = str(self.record.seq)
        seq_id = self.record.id
        # regions = self.annotations.get(
        #     seq_id, [{"start": 0, "end": len(seq), "type": "exon"}]
        # )
        position_map = {}
        for region in self.regions:
            start, end, muts = (
                region["start"],
                region["end"],
                region.get("mutations", []),
            )
            if "SNP" not in muts:
                continue
            for pos in range(
                start, end
            ):  # +1 because we want to insert after the last base as well
                position_map.setdefault(pos, []).append(region)
        mutable_positions = sorted(position_map.keys())
        for pos in mutable_positions:
            rtype = position_map[pos][0]["type"]  # first one in for annotation
            ref = seq[pos]
            for alt in self.bases:
                if alt != ref:
                    mutated_seq = seq[:pos] + alt + seq[pos + 1 :]
                    mutation = Mutation(
                        seq_id=seq_id,
                        genomic_id=self.genomic_id,
                        chromosome=self.chromosome,
                        strand=self.strand,
                        region_type=rtype,
                        mutation_pos=pos,
                        genomic_pos=self.genomic_start + pos,
                        cds_pos=None,
                        ref=ref,
                        alt=alt,
                        mutation_type="SNP",
                        seq=mutated_seq,
                        ref_codon=None,
                        alt_codon=None,
                        ref_aa=None,
                        alt_aa=None,
                    )
                    mutation = self.annotator.annotate(mutation)
                    variants.append(mutation)
        return variants

    # -----------------------------
    # InDels
    # -----------------------------
    def generate_insertion(self, max_len=3) -> list[Mutation]:
        variants = []
        seq = str(self.record.seq)
        seq_id = self.record.id
        position_map = {}
        for region in self.regions:
            start, end, muts = (
                region["start"],
                region["end"],
                region.get("mutations", []),
            )
            if "INS" not in muts:
                continue
            for pos in range(
                start, end + 1
            ):  # +1 because we want to insert after the last base as well
                position_map.setdefault(pos, []).append(region)
        mutable_positions = sorted(position_map.keys())
        for pos in mutable_positions:
            rtype = position_map[pos][0]["type"]  # first one in for annotation
            # Insertion
            for length in range(1, max_len + 1):
                for ins in itertools.product(self.bases, repeat=length):
                    ins_seq = "".join(ins)
                    mutated_seq = seq[:pos] + ins_seq + seq[pos:]
                    mutation_type = f"INS_{length}"
                    mutation = Mutation(
                        seq_id=seq_id,
                        genomic_id=self.genomic_id,
                        chromosome=self.chromosome,
                        strand=self.strand,
                        region_type=rtype,
                        mutation_pos=pos,
                        genomic_pos=self.genomic_start + pos,
                        cds_pos=None,
                        ref="",
                        alt=ins_seq,
                        mutation_type=mutation_type,
                        seq=mutated_seq,
                        ref_codon=None,
                        alt_codon=None,
                        ref_aa=None,
                        alt_aa=None,
                    )
                    mutation = self.annotator.annotate(mutation)
                    variants.append(mutation)
        return variants

    def generate_deletion(self, max_len=3) -> list[Mutation]:
        variants = []
        seq = str(self.record.seq)
        seq_id = self.record.id

        position_map = {}
        region_map = {}

        for region in self.regions:
            start, end, muts = (
                region["start"],
                region["end"],
                region.get("mutations", []),
            )

            if "DEL" not in muts:
                continue

            for pos in range(start, end):
                position_map.setdefault(pos, []).append(region)
                region_map[pos] = region  # overwrite ok or choose priority rule

        mutable_positions = sorted(position_map.keys())
        for pos in mutable_positions:

            for length in range(1, max_len + 1):

                if not self.is_valid_window(pos, length, position_map):
                    continue

                ref_seq = seq[pos : pos + length]
                mutated_seq = seq[:pos] + seq[pos + length :]

                rtype = region_map[pos]["type"]

                mutation_type = f"DEL_{length}"
                mutation = Mutation(
                    seq_id=seq_id,
                    genomic_id=self.genomic_id,
                    chromosome=self.chromosome,
                    strand=self.strand,
                    region_type=rtype,
                    mutation_pos=pos,
                    genomic_pos=self.genomic_start + pos,
                    ref=ref_seq,
                    alt="",
                    mutation_type=mutation_type,
                    seq=mutated_seq,
                    cds_pos=None,
                    ref_codon=None,
                    alt_codon=None,
                    ref_aa=None,
                    alt_aa=None,
                )
                mutation = self.annotator.annotate(mutation)
                variants.append(mutation)

        return variants

    def collect_mutations(self) -> list[Mutation]:
        """
        Collects all generated mutations (SNPs, AA changes, InDels) into a single list.
        """
        mutations = []
        mutations.extend(self.generate_snp())
        mutations.extend(self.generate_aa())
        mutations.extend(self.generate_insertion(self.max_len))
        mutations.extend(self.generate_deletion(self.max_len))
        return mutations

    # -----------------------------
    # Combine all variants
    # -----------------------------
    def combine_variant_tables(self, variants: list[Mutation]) -> DataFrame:
        combined = []
        for mutation in variants:
            combined.append(mutation.__dict__)
        return pd.DataFrame(combined)

    def _exon_codon_index(self, pos):
        exon_start = 0
        for region in self.regions:
            if (
                region["type"] == "exon"
                and region["start"] <= pos < region["end"]
            ):
                exon_start = region["start"]
                break
        return (pos - exon_start) // 3 + 1

    ########################################################################################
    # Elements
    ########################################################################################

    def mark_duplicates(
        self,
        mutations: DataFrame,
        order: list[str] | None = None,
        ascending: bool = False,
    ) -> DataFrame:
        order = order or [
            "missense",
            "nonsense",
            "INS_1",
            "SNP",
            "INS_2",
            "INS_3",
            "INS_fs",
            "DEL_1",
            "DEL_2",
            "DEL_3",
            "DEL_fs",
        ]
        mutations["mutation_type"] = pd.Categorical(
            mutations["mutation_type"], categories=order, ordered=True
        )
        mutations = mutations.sort_values(
            by=["mutation_type", "genomic_pos"],
            ascending=[True, ascending],
        )
        mutations["is_duplicate"] = mutations.duplicated(
            subset="seq", keep="first"
        )
        return mutations

    def generate_mutations(self) -> DataFrame:
        variants = self.collect_mutations()
        df = self.combine_variant_tables(variants)
        df = self.mark_duplicates(df)
        return df

    @staticmethod
    def deduplicate(mutations: DataFrame) -> DataFrame:
        mutations = mutations.drop_duplicates(subset="seq")
        return mutations

    def write_mutations(
        self, output_file: Path | str, deduplicate: bool = False
    ) -> None:
        # the subroutine
        df = self.generate_mutations()
        if deduplicate:
            df = self.consolidate(df)
        df.to_csv(output_file, sep="\t", index=False)


class MutatorHGVS:
    def __init__(self):
        self.bases = ["A", "C", "G", "T"]
        self.types = ["AA", "SNP", "IN", "DEL"]

    def read_records(self, fasta: Path | str) -> dict[str, SeqIO.SeqRecord]:
        return SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    def region_definitions(
        self, json_file: Path | str | None, records
    ) -> dict[str, list[dict[str, Any]]]:
        plan = from_json(Path(json_file)) if json_file else {}
        region_definitions = {}
        for seq_id, record in records.items():
            default = [
                {
                    "name": seq_id,
                    "genomic_start": 0,
                    "start": 0,
                    "end": len(record),
                    "type": "exon",
                    "mutations": self.types,
                    "strand": "+",
                    "chromosome": 0,
                    "length": len(record),
                }
            ]
            region_definitions[seq_id] = plan.get(seq_id, default)
        return region_definitions

    def record_annotations(self, json_file: Path | str | None, records) -> dict:
        annotations = from_json(Path(json_file)) if json_file else {}
        result = {}
        for seq_id in records:
            result[seq_id] = annotations.get(seq_id, None)
        return result

    def select_annotator(
        self, annotator_to_use: Literal["manual", "mave", "hgvs"]
    ) -> MutationAnnotator:
        if annotator_to_use == "manual":
            return ManualMutationAnnotator()
        elif annotator_to_use == "mave":
            return MaveMutationAnnotator()
        elif annotator_to_use == "hgvs":
            return HGVSMutationAnnotator()
        else:
            raise ValueError(f"Annotator {annotator_to_use} not implemented")

    def create_all_mutations(
        self,
        fasta: Path | str,
        region_json: Path | str,
        cds_id: str,
        genomic_id: str = "17",
        annotation_json: Path | str | None = None,
        annotator_to_use: Literal["manual", "mave", "hgvs"] = "manual",
        deduplicate: bool = False,
    ) -> DataFrame:
        records = self.read_records(fasta)
        regions = self.region_definitions(region_json, records)
        annotations = self.record_annotations(annotation_json, records)
        # annotator = self.select_annotator(annotator_to_use)
        to_concat = []
        cds = records[cds_id]
        genomic = records[genomic_id]
        for seq_id, record in records.items():
            if seq_id == cds_id:
                continue
            mutator = MultiSeqMutatorHGVS(
                record[seq_id],
                regions[seq_id],
                genomic,
                cds,
                annotations[seq_id],
                self.select_annotator(annotator_to_use),
            )
            mutations = mutator.generate_mutations()
            if deduplicate:
                mutations = mutator.deduplicate(mutations)
            to_concat.append(mutations)
        combined = pd.concat(to_concat, ignore_index=True)
        return combined

    def generate_mutations(
        self,
        output_file: str | Path,
        fasta: Path | str,
        region_json: Path | str,
        cds_id: str,
        genomic_id: str,
        annotation_json: Path | str | None = None,
        annotator_to_use: Literal["manual", "mave", "hgvs"] = "manual",
        deduplicate: bool = False,
    ) -> None:
        combined = self.create_all_mutations(
            fasta,
            region_json,
            cds_id,
            genomic_id,
            annotation_json,
            annotator_to_use,
            deduplicate,
        )
        combined.to_csv(output_file, sep="\t", index=False)

    def generate_from_clinvar(
        self,
        output_file: Path | str,
        fasta: Path | str,
        clinvar_cases_file: Path | str,
        cds_id: str,
        genomics_id: str = "17",
        # transcript_ac="NM_000546.6",
        # chrom_ac="NC_000017.11",
        cdot_json: Path | str = "/incoming/cdot-0.2.21.refseq.grch38_tp53.json",
        # annotation_json: Path | str | None = None,
        # seq_id_prefix="var",
        # genomic_flanks=(
        #     7673511,
        #     7673636,
        # ),  # checked in ensembl browser for exon 9 part of TP53
    ) -> None:
        clinvar_cases = pd.read_csv(clinvar_cases_file, sep="\t")
        records = self.read_records(fasta)
        genomic = records.get(genomics_id)
        cds = records.get(cds_id)
        if genomic is None:
            raise ValueError(
                f"Genomic sequence with id {genomics_id} not found in fasta"
            )
        if cds is None:
            raise ValueError(
                f"CDS sequence with id {cds_id} not found in fasta"
            )

        # annotations = self.record_annotations(annotation_json, records)
        generator = HGVSVariantGenerator(
            genomic=genomic,
            cds=cds,
            # transcript_ac=transcript_ac,
            # chrom_ac=chrom_ac,
            cdot_json=cdot_json,
            seq_id_prefix="predvar",
            # genomic_flanks=(7673511, 7673636),  # checked in ensembl browser
        )

        variants = []
        for _, row in clinvar_cases.iterrows():
            hgvs_c = row["corrected_hgvs_c"]
            genomic_start = row["genomic_start"]
            genomic_end = row["genomic_end"]
            variant = generator.from_hgvs(hgvs_c, (genomic_start, genomic_end))
            variants.append(variant)
        df = pd.DataFrame([v.__dict__ for v in variants])
        df.to_csv(output_file, sep="\t", index=False)

    @element
    def generate(
        self,
        fasta: FileElement | str | Path,
        regions: FileElement | str | Path | None = None,
        annotations: FileElement | str | Path | None = None,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
    ) -> Element:
        """
        High-level function to generate an output file of in silico mutations.

        Parameters
        ----------
        fasta : FileElement
            fasta file containing the reference sequences to be mutated.
        region_json : FileElement | str | Path | None, optional
            JSON file or path containing the region definitions, by default None. If
            None, default regions will be used.
        annotation_json : FileElement | str | Path | None, optional
            JSON file or path containing the sequence annotations, by default None. If
            None, default annotations will be used.
        tag : Tag | ElementTag | None, optional
            Tag for the mutation run, by default None
        outdir : Path | str | None, optional
            Output directory for the mutation results, by default None
        filename : Path | str | None, optional
            Filename for the mutation results, by default None
        params : Params | None, optional
            Parameters for the mutation run, by default None

        Returns
        -------
        Element
            An element that performs the mutation generation and produces an output file
            with the mutations.
        """
        # fasta_file = (
        #     Path(fasta) if not isinstance(fasta, FileElement) else fasta.fasta
        # )
        # regions_json_file = None
        # annotations_json_file = None
        # if regions:
        #     regions_json_file = (
        #         regions.json
        #         if isinstance(regions, FileElement)
        #         else Path(regions)
        #     )
        # if annotations:
        #     annotations_json_file = (
        #         annotations.json
        #         if isinstance(annotations, FileElement)
        #         else Path(annotations)
        #     )

        # outdir = Path(outdir or fasta_file.parent / "mutations")
        # if isinstance(fasta, Element):
        #     tag = from_prior(
        #         fasta.tag,
        #         tag,
        #         stage=Stage.MUTATE,
        #         method=Method.CUSTOM,
        #         state=State.GENERATED,
        #         ext="tsv",
        #     )
        # else:
        #     tag = ElementTag(
        #         stage=Stage.MUTATE,
        #         method=Method.CUSTOM,
        #         state=State.GENERATED,
        #         ext="tsv",
        #     )
        # outfile = outdir / (filename or tag.default_output).absolute()

        # runner = self.create_all_mutations(
        #     fasta=fasta_file,
        #     output_file=outfile,
        #     plan_json=plan_json,
        # )
        # determinants = ()
        # pres = (fasta,) if isinstance(fasta, Element) else ()
        # if isinstance(plan_json, Element):
        #     pres += (plan_json,)
        # key, name = generate_element_key_name(
        #     tag,
        #     "0.0.1",
        #     "MutatorHGVS",
        # )

        # return Element(
        #     key=key,
        #     name=name,
        #     run=runner,
        #     tag=tag,
        #     determinants=determinants,
        #     inputs=(fasta_file, plan_json) if plan_json else (fasta_file,),
        #     artifacts={
        #         "tsv": outfile,
        #     },
        #     pres=pres,
        # )


if __name__ == "__main__":
    major = MutatorHGVS()
    records = major.read_records("incoming/exons.fasta")
    regions = major.region_definitions("incoming/regions.json", records)
    print(records["TP53-201_cds_protein_coding"])
    cds = records["TP53-201_cds_protein_coding"]
    sample = records["TP53-201_Ex9_Flanks"]
    print(regions)
    mutator = MultiSeqMutatorHGVS(sample, regions["TP53-201_Ex9_Flanks"], cds)
    print(mutator.codon_table)
    print(
        mutator.get_alternate_codons(
            codon="ATG", pos_in_codon=0, constant_front=True
        )
    )
    print(
        mutator.get_minmimal_hamming_synonyms(
            codon="ATG", pos_in_codon=0, constant_front=True
        )
    )
    print("--------------------------------")
    print(
        mutator.get_alternate_codons(
            codon="ATG", pos_in_codon=1, constant_front=True
        )
    )
    print(
        mutator.get_minmimal_hamming_synonyms(
            codon="ATG", pos_in_codon=1, constant_front=True
        )
    )
    print("--------------------------------")
    print(
        mutator.get_alternate_codons(
            codon="ATG", pos_in_codon=1, constant_front=False
        )
    )
    print(
        mutator.get_minmimal_hamming_synonyms(
            codon="ATG", pos_in_codon=1, constant_front=False
        )
    )
