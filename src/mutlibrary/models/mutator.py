import itertools
import sys
from itertools import product
from pathlib import Path
from typing import Any, Literal

import pandas as pd  # type: ignore[import]
from Bio import Align, SeqIO  # type: ignore[import]
from Bio.Data import CodonTable  # type: ignore[import]
from Bio.Seq import Seq  # type: ignore[import]
from Bio.SeqIO import SeqRecord  # type: ignore[import]
from mmalignments.services.io import from_json, parents  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

from mutlibrary.models.annotator import (  # type: ignore[import]
    HGVSMutationAnnotator,
    # ManualMutationAnnotator,
    # Mutalyzer,
    Mutation,
    MutationAnnotator,
    get_info_from_description,
)
from mutlibrary.models.clinvar import (
    HGVSVariantGenerator,  # type: ignore[import]
)


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


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
        record: SeqRecord,
        genomic: SeqRecord,
        cds: SeqRecord,
        regions: list[dict[str, Any]],
        annotation: dict | Path | str,
        annotator: MutationAnnotator | None = None,
        cdot_json: Path | None = None,
    ):
        self.record = record
        info = get_info_from_description(record)
        self.genomic_start = int(info["start"])
        self.genomic_end = int(info["end"])
        self.chromosome = info["chromosome"]
        self.strand = int(info["strand"])
        self.assembly = info["assembly"]
        self.regions = regions
        self.cds = cds
        self.genomic_id = genomic.id if genomic else ""
        self.genomic = genomic
        if isinstance(annotation, (str, Path)):
            self.annotation = MultiSeqMutatorHGVS.annotations(annotation)
        elif isinstance(annotation, dict):
            self.annotation = annotation
        else:
            raise ValueError(
                "Annotation must be a dict or a path to a json file"
            )
        self.transcript_id = self.annotation[self.cds.id].get(
            "transcript_id", ""
        )
        self.transcript_ac = self.annotation[self.cds.id].get(
            "transcript_ac", ""
        )
        self.chromosome_ac = self.annotation[self.cds.id].get(
            "chromosome_ac", ""
        )
        self.protein_ac = self.annotation[self.cds.id].get("protein_ac", "")
        self.chromosome_ac = self.annotation[self.genomic.id].get(
            "chromosome_ac", ""
        )

        self.annotator = annotator or HGVSMutationAnnotator(
            assembly=self.assembly, cdot_json=cdot_json
        )  # Mutalyzer()
        self.bases = ["A", "C", "G", "T"]
        self.codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
        self.codon_map = dict(self.codon_table.forward_table)
        for stop in self.codon_table.stop_codons:
            self.codon_map[stop] = "*"
        self.combined_table = Path(f"cache/{self.record.id}_variants.tsv")
        self.cds_codons = self.build_cds_codon_mapping()
        self.max_len = 3

    ############################################################################
    # Helper
    ############################################################################
    @classmethod
    def annotations(
        cls, json_file: Path | str | None
    ) -> dict[str, dict[str, Any]]:
        return from_json(Path(json_file)) if json_file else {}

    def get_genomic_pos_from_local(self, local_pos: int) -> int:
        if self.strand == 1:
            return self.genomic_start + local_pos
        else:
            return self.genomic_end - local_pos  # - 1

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
        print("mapping", mapping)
        mutatable_codon_positions = []
        start, end = region["start"], region["end"]
        exon = seq[start:end]
        print("exon", exon)
        cds_start = self.map_exon_to_cds(exon)
        # CDS Positionen, die in diesem Exon liegen
        print("cds_start", cds_start)
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
                            True,  # keep front constant if truncated at start
                        )
                    )
                    break
        print("cds_positions", cds_positions)
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
                    print(
                        exon_pos,
                        cds_pos,
                        cds_pos_end,
                        last_covered_cds_pos,
                        in_exon,
                        codon,
                    )
                else:
                    # full codon in exo
                    in_exon = exon[exon_pos : min(exon_pos + 3, len(exon))]
                    mutatable_codon_positions.append(
                        (
                            cds_pos,
                            exon_pos,  # position in exon
                            in_exon,  # codon in exon
                            codon,  # codon in frame, should be the same
                            0,  # position in codon if truncated, this is != 0
                            True,  # keep the front constant if truncated
                        )
                    )
                    print(
                        exon_pos,
                        cds_pos,
                        cds_pos_end,
                        last_covered_cds_pos,
                        in_exon,
                        codon,
                    )
                # last_covered_cds_pos = cds_pos
        return mutatable_codon_positions

    def is_valid_window(self, pos, length, position_map):
        for i in range(pos, pos + length):
            if i not in position_map:
                return False
        return True

    ############################################################################
    # The sequence generators
    ############################################################################

    def generate_aa_mutseq(
        self,
        pos_in_codon: int,
        seq: str,
        current_start: int,
        alt_codon: str,
        codon: str,
        exon_subseq: str,
        constant_front: bool,
    ) -> tuple[str, str, str]:
        if pos_in_codon == 0:  # full_codon
            mutated_seq = str(
                seq[:current_start] + alt_codon + seq[current_start + 3 :]
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
            raise ValueError("No full codon and no truncated codon?")
        return mutated_seq, ref, alt

    def generate_aa(self) -> list[Mutation]:
        # Amino Acid changes (exons only)
        if self.cds is None:
            raise ValueError(
                "CDS sequence is required for amino acid mutation generation"  # noqa: E501
            )

        variants = []
        seq = str(self.record.seq)
        seq_id = self.record.id

        for region in self.regions:
            if region["type"] != "exon" or "AA" not in region.get(
                "mutations", []
            ):
                continue

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

                    mutated_seq, ref, alt = self.generate_aa_mutseq(
                        pos_in_codon,
                        seq,
                        current_start,
                        alt_codon,
                        codon,
                        exon_subseq,
                        constant_front,
                    )
                    print("the local index:", current_start)
                    print(
                        exon_pos,
                        cds_pos,
                        exon_subseq,
                        codon,
                        pos_in_codon,
                        constant_front,
                    )

                    # mutation type
                    mutation_type = "nonsense" if alt_aa == "*" else "missense"
                    genomic_seq = mutated_seq
                    if self.strand == -1:
                        genomic_seq = reverse_complement(mutated_seq)
                        ref = reverse_complement(ref)
                        alt = reverse_complement(alt)
                    mutation = Mutation(
                        seq_id=seq_id,
                        genomic_id=self.genomic_id,
                        coding=mutated_seq,
                        genomic=genomic_seq,
                        refseq=seq,
                        chromosome=self.chromosome,
                        genomic_start=self.genomic_start,
                        genomic_end=self.genomic_end,
                        strand=self.strand,
                        mutation_pos=current_start,
                        genomic_pos=self.get_genomic_pos_from_local(
                            current_start
                        ),
                        region_type="exon",
                        ref=ref,
                        alt=alt,
                        mutation_type=mutation_type,
                        cds_pos=cds_pos,
                        ref_codon=codon,
                        alt_codon=alt_codon,
                        ref_aa=ref_aa,
                        alt_aa=alt_aa,
                        chromosome_ac=self.chromosome_ac,
                        transcript_ac=self.transcript_ac,
                        protein_ac=self.protein_ac,
                    )
                    mutation = self.annotator.annotate(mutation)
                    variants.append(mutation)

        return variants

    # -----------------------------
    # SNPs
    # -----------------------------
    def generate_snp(self) -> list[Mutation]:
        """
        Generates a list of all single SNP mutations for the record
        sequence. Each specified region in regions can have different mutation
        types specified, but the complete sequence is returned.
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
            ref_base = seq[pos]
            for alt_base in self.bases:
                if alt_base != ref_base:
                    mutated_seq = seq[:pos] + alt_base + seq[pos + 1 :]
                    try:
                        assert (
                            mutated_seq != seq
                        ), f"Mutation did not change the sequence at position {pos}: {ref_base} -> {alt_base}"  # noqa: E501
                    except AssertionError:
                        print("position", pos, seq)
                        print(
                            "whaaat!",
                            seq[:pos],
                            ref_base,
                            alt_base,
                            seq[pos + 1 :],
                        )
                        print(
                            "wrong=", ref_base, alt_base, seq, mutated_seq, pos
                        )
                        print("at 0 is:", seq[0], "at pos is:", seq[pos])

                        raise
                    genomic_seq = mutated_seq
                    ref, alt = ref_base, alt_base
                    if self.strand == -1:
                        genomic_seq = reverse_complement(mutated_seq)
                        ref = reverse_complement(ref)
                        alt = reverse_complement(alt)

                    mutation = Mutation(
                        seq_id=seq_id,
                        genomic_id=self.genomic_id,
                        genomic=genomic_seq,
                        coding=mutated_seq,
                        refseq=seq,
                        chromosome=self.chromosome,
                        genomic_start=self.genomic_start,
                        genomic_end=self.genomic_end,
                        strand=self.strand,
                        region_type=rtype,
                        mutation_pos=pos,
                        genomic_pos=self.get_genomic_pos_from_local(pos),
                        cds_pos=None,
                        ref=ref,
                        alt=alt,
                        mutation_type="SNP",
                        ref_codon=None,
                        alt_codon=None,
                        ref_aa=None,
                        alt_aa=None,
                        chromosome_ac=self.chromosome_ac,
                        transcript_ac=self.transcript_ac,
                        protein_ac=self.protein_ac,
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
                    genomic_seq = mutated_seq
                    alt = ins_seq
                    if self.strand == -1:
                        genomic_seq = reverse_complement(mutated_seq)
                        alt = reverse_complement(alt)
                    mutation_type = f"INS_{length}"
                    try:
                        assert (
                            mutated_seq != seq
                        ), f"Mutation did not change the sequence at position {pos}: {"ref_base"} -> {alt}"  # noqa: E501
                    except AssertionError:
                        print("position", pos, seq)
                        print(
                            "whaaat!",
                            seq[:pos],
                            alt,
                            seq[pos + 1 :],
                        )
                        print("wrong=", alt, seq, mutated_seq, pos)
                        print("at 0 is:", seq[0], "at pos is:", seq[pos])

                        raise
                    print()
                    mutation = Mutation(
                        seq_id=seq_id,
                        genomic_id=self.genomic_id,
                        coding=mutated_seq,
                        genomic=genomic_seq,
                        refseq=seq,
                        chromosome=self.chromosome,
                        genomic_start=self.genomic_start,
                        genomic_end=self.genomic_end,
                        strand=self.strand,
                        region_type=rtype,
                        mutation_pos=pos,
                        genomic_pos=self.get_genomic_pos_from_local(pos),
                        cds_pos=pos,
                        ref="",
                        alt=alt,
                        mutation_type=mutation_type,
                        ref_codon=None,
                        alt_codon=None,
                        ref_aa=None,
                        alt_aa=None,
                        chromosome_ac=self.chromosome_ac,
                        transcript_ac=self.transcript_ac,
                        protein_ac=self.protein_ac,
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
                genomic_seq = (
                    mutated_seq
                    if self.strand == 1
                    else str(Seq(mutated_seq).reverse_complement())
                )
                genomic_seq = mutated_seq
                ref = ref_seq
                if self.strand == -1:
                    genomic_seq = reverse_complement(mutated_seq)
                    ref = reverse_complement(ref)
                mutation = Mutation(
                    seq_id=seq_id,
                    genomic_id=self.genomic_id,
                    coding=mutated_seq,
                    genomic=genomic_seq,
                    refseq=seq,
                    chromosome=self.chromosome,
                    genomic_start=self.genomic_start,
                    genomic_end=self.genomic_end,
                    strand=self.strand,
                    region_type=rtype,
                    mutation_pos=pos,
                    genomic_pos=self.get_genomic_pos_from_local(pos),
                    ref=ref,
                    alt="",
                    mutation_type=mutation_type,
                    cds_pos=pos,
                    ref_codon=None,
                    alt_codon=None,
                    ref_aa=None,
                    alt_aa=None,
                    chromosome_ac=self.chromosome_ac,
                    transcript_ac=self.transcript_ac,
                    protein_ac=self.protein_ac,
                )
                mutation = self.annotator.annotate(mutation)
                variants.append(mutation)

        return variants

    def collect_mutations(self) -> list[Mutation]:
        """
        Collects all generated mutations (SNPs, AA changes, InDels) into a
        single list.
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

        combined_df = pd.DataFrame(combined)
        combined_df = combined_df.sort_values(
            by=["seq_id", "mutation_type", "genomic_pos"],
            ascending=[True, True, True],
        )
        mut_ids = (
            combined_df["seq_id"] + "_mut_" + combined_df.index.astype(str)
        )
        combined_df.insert(0, "mut_id", mut_ids)
        return combined_df

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

    ############################################################################
    # Elements
    ############################################################################

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
            subset="coding", keep="first"
        )
        return mutations

    def generate_mutations(self) -> DataFrame:
        variants = self.collect_mutations()
        df = self.combine_variant_tables(variants)
        df = self.mark_duplicates(df)
        return df

    @staticmethod
    def deduplicate(mutations: DataFrame) -> DataFrame:
        mutations = DataFrame(mutations[~mutations["is_duplicate"]].copy())
        return mutations

    def write_mutations(
        self, output_file: Path | str, deduplicate: bool = False
    ) -> None:
        # the subroutine
        df = self.generate_mutations()
        if deduplicate:
            df = self.deduplicate(df)
        df.to_csv(output_file, sep="\t", index=False)


class MutatorHGVS:
    def __init__(self):
        self.bases = ["A", "C", "G", "T"]
        self.types = ["AA", "SNP", "IN", "DEL"]

    def read_records(self, fasta: Path | str) -> dict[str, SeqIO.SeqRecord]:
        return SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    def region_definitions(
        self,
        json_file: Path | str | None,
        records: dict[str, SeqIO.SeqRecord],
        exclude: list[str] | None = None,
    ) -> dict[str, list[dict[str, Any]]]:
        plan = from_json(Path(json_file)) if json_file else {}
        region_definitions = {}
        for seq_id, record in records.items():
            if seq_id in (exclude or []):
                continue
            default = [
                {
                    "chromosome": 0,
                    "strand": 1,
                    "genomic_start": 0,
                    "genomic_end": len(record),
                    "assembly": "",
                    "transcript_id": "",
                    "transcript_ac": "",
                    "chromosome_ac": "",
                    "protein_ac": "",
                }
            ]
            default = None
            region_definitions[seq_id] = plan.get(seq_id, default)
        return region_definitions

    @classmethod
    def record_annotations(
        cls, json_file: Path | str | None, records: dict | None = None
    ) -> dict:
        annotations = from_json(Path(json_file)) if json_file else {}
        result = {}
        ids = records.keys() if records else annotations.keys()
        for seq_id in ids:
            result[seq_id] = annotations.get(seq_id, None)
        return result

    def select_annotator(
        self, annotator_to_use: Literal["manual", "mave", "hgvs", "mutalyzer"]
    ) -> MutationAnnotator:
        if annotator_to_use == "hgvs":
            return HGVSMutationAnnotator()
        # elif annotator_to_use == "manual":
        #     return ManualMutationAnnotator()
        # elif annotator_to_use == "mave":
        #     return MaveMutationAnnotator()
        # elif annotator_to_use == "mutalyzer":
        #     return Mutalyzer()
        else:
            raise ValueError(f"Annotator {annotator_to_use} not implemented")

    def create_all_mutations(
        self,
        fasta: Path | str,
        region_json: Path | str,
        cds_id: str,
        genomic_id: str = "17",
        annotation_json: Path | str | None = None,
        annotator_to_use: Literal[
            "manual", "mave", "hgvs", "mutalyzer"
        ] = "manual",
        deduplicate: bool = False,
    ) -> DataFrame:
        records = self.read_records(fasta)
        regions = self.region_definitions(
            region_json, records, exclude=[cds_id, genomic_id]
        )
        annotations = self.record_annotations(annotation_json, records)
        annotations = (
            from_json(Path(annotation_json)) if annotation_json else {}
        )
        # annotator = self.select_annotator(annotator_to_use)
        to_concat = []
        cds = records[cds_id]
        genomic = records[genomic_id]
        for seq_id, record in records.items():
            if seq_id in [cds_id, genomic_id]:
                continue
            # print(
            #     "Creating mutator for",
            #     seq_id,
            #     record,
            #     regions[seq_id],
            #     genomic,
            #     cds,
            #     annotations[seq_id],
            # )
            # print("annotations", annotations)
            mutator = MultiSeqMutatorHGVS(
                record,
                genomic,
                cds,
                regions[seq_id],
                annotations,
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
        annotator_to_use: Literal[
            "manual", "mave", "hgvs", "mutalyzer"
        ] = "manual",
        deduplicate: bool = False,
    ) -> None:
        parents(output_file)
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
        clinvar_cases_file: Path | str,
        cdot_json: Path | str,
        fasta: Path | str,
        cds_id: str,
        genomic_id: str = "17",
        annotation_json: Path | str | None = None,
        clinvar_column: str = "Accession",
        input_column: str = "corrected_hgvs_c",
        genomic_start_column: str = "genomic_start",
        genomic_end_column: str = "genomic_end",
    ) -> None:
        parents(output_file)
        clinvar_cases = pd.read_csv(clinvar_cases_file, sep="\t")
        records = self.read_records(fasta)
        genomic = records.get(genomic_id)
        cds = records.get(cds_id)
        if genomic is None:
            raise ValueError(
                f"Genomic sequence with id {genomic_id} not found in fasta"
            )
        if cds is None:
            raise ValueError(
                f"CDS sequence with id {cds_id} not found in fasta"
            )

        generator = HGVSVariantGenerator(
            genomic=genomic,
            cds=cds,
            cdot_json=cdot_json,
            annotation_json=annotation_json,
        )
        variants = []
        # clinvar_ids = [
        # generator.vcv_to_variation_id(v)
        # for v in clinvar_cases[clinvar_column]
        # ]
        # spid = generator.get_spdi_from_clinvar("VCV000428872")
        # raise ValueError("Stop here for now")
        # spdis = [
        #     generator.get_spdi_from_clinvar(vcv)
        #     for vcv in clinvar_cases[clinvar_column]
        # ]
        # spdis = ["NC_000017.11:7669697:AGACAGA:AGA"]

        # clinvar_cases["spdi"] = spdis
        # vcvs = ["VCV001116590"]
        # hghvs = ["NM_000546.6:c.1101-11_1101-8del"]
        # for _, row in hghvs:
        # variant = generator.from_vcv(vcv, (7669584, 7669715))
        # variant = generator.from_spdi(spdi, (7669584, 7669715), "VCV000428872")
        # variant = generator.from_hgvs(
        #     vcv, (7669584, 7669715), "VCV000428872"
        # )
        # print(variant, variant.coding, variant.genomic)
        # variants.append(variant)
        for _, row in clinvar_cases.iterrows():
            hgvs_c = row[input_column]
            genomic_start = row[genomic_start_column]
            genomic_end = row[genomic_end_column]
            vcv = row[clinvar_column]
            variant = generator.from_hgvs(
                hgvs_c, (genomic_start, genomic_end), vcv
            )
            variants.append(variant)
        df = pd.DataFrame([v.__dict__ for v in variants])
        mut_ids = df["seq_id"] + "_mut_" + df.index.astype(str)
        df.insert(0, "mut_id", mut_ids)
        df.to_csv(output_file, sep="\t", index=False)

    def generate(
        self,
        output_combined: str | Path,
        output_file_generated: str | Path,
        output_file_clinvar: str | Path,
        fasta: Path | str,
        region_json: Path | str,
        clinvar_cases_file: Path | str,
        cdot_json: Path | str,
        cds_id: str,
        genomic_id: str,
        annotation_json: Path | str | None = None,
        annotator_to_use: Literal[
            "manual", "mave", "hgvs", "mutalyzer"
        ] = "manual",
        deduplicate: bool = False,
        input_column: str = "corrected_hgvs_c",
        genomic_start_column: str = "genomic_start",
        genomic_end_column: str = "genomic_end",
    ) -> None:
        self.generate_mutations(
            output_file_generated,
            fasta=fasta,
            region_json=region_json,
            cds_id=cds_id,
            genomic_id=genomic_id,
            annotation_json=annotation_json,
            annotator_to_use=annotator_to_use,
            deduplicate=deduplicate,
        )
        self.generate_from_clinvar(
            output_file_clinvar,
            clinvar_cases_file=clinvar_cases_file,
            cdot_json=cdot_json,
            fasta=fasta,
            cds_id=cds_id,
            genomic_id=genomic_id,
            annotation_json=annotation_json,
            input_column=input_column,
            genomic_start_column=genomic_start_column,
            genomic_end_column=genomic_end_column,
        )
        df_combined = pd.concat(
            [
                pd.read_csv(output_file_generated, sep="\t"),
                pd.read_csv(output_file_clinvar, sep="\t"),
            ],
            ignore_index=True,
        )
        df_combined.to_csv(output_combined, sep="\t", index=False)

    #     output_file: Path | str,
    #     fasta: Path | str,
    #     cds_id: str,
    #     genomic_id: str = "17",
    #     cdot_json: Path | str = "/incoming/cdot-0.2.21.refseq.grch38_tp53.json",
    #     annotation_json: Path | str | None = None,
    #     input_column: str = "corrected_hgvs_c",
    #     genomic_start_column: str = "genomic_start",
    #     genomic_end_column: str = "genomic_end",
    # )
    # @element
    # def generate(
    #     self,
    #     fasta: FileElement | str | Path,
    #     regions: FileElement | str | Path | None = None,
    #     annotations: FileElement | str | Path | None = None,
    #     *,
    #     tag: PartialElementTag | ElementTag | None = None,
    #     outdir: Path | str | None = None,
    #     filename: Path | str | None = None,
    #     params: Params | None = None,
    # ) -> Element:
    #     """
    #     High-level function to generate an output file of in silico mutations.

    #     Parameters
    #     ----------
    #     fasta : FileElement
    #         fasta file containing the reference sequences to be mutated.
    #     region_json : FileElement | str | Path | None, optional
    #         JSON file or path containing the region definitions, by default
    # None. If
    #         None, default regions will be used.
    #     annotation_json : FileElement | str | Path | None, optional
    #         JSON file or path containing the sequence annotations, by default
    # None. If
    #         None, default annotations will be used.
    #     tag : Tag | ElementTag | None, optional
    #         Tag for the mutation run, by default None
    #     outdir : Path | str | None, optional
    #         Output directory for the mutation results, by default None
    #     filename : Path | str | None, optional
    #         Filename for the mutation results, by default None
    #     params : Params | None, optional
    #         Parameters for the mutation run, by default None

    #     Returns
    #     -------
    #     Element
    #         An element that performs the mutation generation and produces an
    #  output file
    #         with the mutations.
    #     """
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


# if __name__ == "__main__":
#     major = MutatorHGVS()
#     records = major.read_records("incoming/exons.fasta")
#     regions = major.region_definitions("incoming/regions.json", records)
#     print(records["TP53-201_cds_protein_coding"])
#     cds = records["TP53-201_cds_protein_coding"]
#     sample = records["TP53-201_Ex9_Flanks"]
#     print(regions)
#     mutator = MultiSeqMutatorHGVS(sample, regions["TP53-201_Ex9_Flanks"], cds)
#     print(mutator.codon_table)
#     print(
#         mutator.get_alternate_codons(
#             codon="ATG", pos_in_codon=0, constant_front=True
#         )
#     )
#     print(
#         mutator.get_minmimal_hamming_synonyms(
#             codon="ATG", pos_in_codon=0, constant_front=True
#         )
#     )
#     print("--------------------------------")
#     print(
#         mutator.get_alternate_codons(
#             codon="ATG", pos_in_codon=1, constant_front=True
#         )
#     )
#     print(
#         mutator.get_minmimal_hamming_synonyms(
#             codon="ATG", pos_in_codon=1, constant_front=True
#         )
#     )
#     print("--------------------------------")
#     print(
#         mutator.get_alternate_codons(
#             codon="ATG", pos_in_codon=1, constant_front=False
#         )
#     )
#     print(
#         mutator.get_minmimal_hamming_synonyms(
#             codon="ATG", pos_in_codon=1, constant_front=False
#         )
#     )
