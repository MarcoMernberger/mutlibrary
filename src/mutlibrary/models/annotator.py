from abc import ABC, abstractmethod
from dataclasses import dataclass, replace
from pathlib import Path

import hgvs.assemblymapper  # type: ignore[import]
import hgvs.dataproviders.uta  # type: ignore[import]
import hgvs.edit  # type: ignore[import]
import hgvs.exceptions  # type: ignore[import]
import hgvs.location  # type: ignore[import]
import hgvs.parser  # type: ignore[import]
import hgvs.posedit  # type: ignore[import]
import hgvs.sequencevariant  # type: ignore[import]
from Bio.SeqIO import SeqRecord  # type: ignore[import]
from cdot.hgvs.dataproviders import (  # type: ignore[import]
    JSONDataProvider,
    RESTDataProvider,
)
from hgvs.exceptions import (  # type: ignore[import]
    HGVSDataNotAvailableError,
    HGVSParseError,
)
from hgvs.sequencevariant import SequenceVariant  # type: ignore[import]


def get_info_from_description(record: SeqRecord):
    desc = record.description
    splits = desc.split(":")
    strand = splits[-1]
    end = splits[-2]
    start = splits[-3]
    chromosome = splits[-4]
    assembly = splits[-5]
    return {
        "strand": strand,
        "end": end,
        "start": start,
        "chromosome": chromosome,
        "assembly": assembly,
    }


@dataclass
class Mutation:
    seq_id: str  # a name for the sequence
    # seq: str  # the mutated sequence on coding strand
    genomic: str  # the mutated sequence on + strand
    genomic_id: str  # the genomic_id of the reference record
    coding: str  # the mutated sequence on coding strand
    refseq: str  # the reference sequence with flanking
    chromosome: str
    genomic_start: int
    genomic_end: int
    strand: int  # +1 or -1
    region_type: str  # intron or exon
    mutation_pos: int  #   the position in the input_sequence
    genomic_pos: int  #   the genomic position of the mutation (1-based)
    cds_pos: int | None  # the position in the CDS (0-based), None if non-coding
    ref: str  # the reference allele (nucleotide(s)), as in vcf
    alt: str  # the alternate allele (nucleotide(s)), as in vcf
    mutation_type: str  # e.g. "SNV", "Del", "Ins", "MNV", "Delins"
    ref_codon: (
        str | None
    )  # the reference codon (3 nucleotides), None if non-coding
    alt_codon: (
        str | None
    )  # the alternate codon (3 nucleotides), None if non-coding
    ref_aa: str | None  # the reference amino acid, None if non-coding
    alt_aa: str | None  # the alternate amino acid, None if non-coding
    chromosome_ac: str  # the genomic_id for normalized hgvs, e.g. NC_000017.11
    transcript_ac: str
    protein_ac: str
    prot_pos: int | None = (
        None  # the position in the protein (1-based), None if non-coding
    )
    hgvs_g: str | None = (
        None  # the genomic HGVS annotation, e.g. "NC_000017.11:g.7674220A>T"
    )
    hgvs_c: str | None = (
        None  # the coding HGVS annotation, e.g. "NM_000546.6:c.215C>G"
    )
    hgvs_p: str | None = (
        None  # the protein HGVS annotation, e.g. "NP_000537.3:p.Pro72Arg"
    )
    hgvs_r: str | None = (
        None  # the RNA HGVS annotation, e.g. "NR_000537.3:r.215C>G"
    )

    def __repr__(self):
        return f"Mutation(seq_id={self.seq_id}, genomic_id={self.genomic_id}, chromosome={self.chromosome}, strand={self.strand}, region_type={self.region_type}, mutation_pos={self.mutation_pos}, genomic_pos={self.genomic_pos}, cds_pos={self.cds_pos}, ref={self.ref}, alt={self.alt}, mutation_type={self.mutation_type}, ref_codon={self.ref_codon}, alt_codon={self.alt_codon}, ref_aa={self.ref_aa}, alt_aa={self.alt_aa}, prot_pos={self.prot_pos}, hgvs_g={self.hgvs_g}, hgvs_c={self.hgvs_c}, hgvs_p={self.hgvs_p})"  # noqa: E501

    def __str__(self):
        return f"Mutation(\nseq_id={self.seq_id},\ngenomic_id={self.genomic_id},\nchromosome={self.chromosome},\nstrand={self.strand},\ngenomic_start={self.genomic_start},\ngenomic_end={self.genomic_end},\nregion_type={self.region_type},\nmutation_pos={self.mutation_pos},\ngenomic_pos={self.genomic_pos},\ncds_pos={self.cds_pos},\nref={self.ref},\nalt={self.alt},\nmutation_type={self.mutation_type},\nref_codon={self.ref_codon},\nalt_codon={self.alt_codon},\nref_aa={self.ref_aa},\nalt_aa={self.alt_aa},\nchromosome_ac={self.chromosome_ac},\ntranscript_ac={self.transcript_ac},\nprotein_ac={self.protein_ac},\nprot_pos={self.prot_pos},\nhgvs_g={self.hgvs_g},\nhgvs_c={self.hgvs_c},\nhgvs_p={self.hgvs_p})"  # noqa: E501


class MutationAnnotator(ABC):

    def __init__(self):
        pass

    @abstractmethod
    def annotate(
        self,
        mutation: Mutation,
    ) -> Mutation:
        raise NotImplementedError()

    def to_hgvs(self, var: SequenceVariant | None) -> str | None:
        if var is None:
            return None
        return str(var)


class HGVSMutationAnnotator(MutationAnnotator):

    # NC_ accessions for GRCh38 (hg38) and GRCh37 (hg19)
    _CHROM_TO_NC_GRCH38: dict[str, str] = {
        "1": "NC_000001.11",
        "2": "NC_000002.12",
        "3": "NC_000003.12",
        "4": "NC_000004.12",
        "5": "NC_000005.10",
        "6": "NC_000006.12",
        "7": "NC_000007.14",
        "8": "NC_000008.11",
        "9": "NC_000009.12",
        "10": "NC_000010.11",
        "11": "NC_000011.10",
        "12": "NC_000012.12",
        "13": "NC_000013.11",
        "14": "NC_000014.9",
        "15": "NC_000015.10",
        "16": "NC_000016.10",
        "17": "NC_000017.11",
        "18": "NC_000018.10",
        "19": "NC_000019.10",
        "20": "NC_000020.11",
        "21": "NC_000021.9",
        "22": "NC_000022.11",
        "X": "NC_000023.11",
        "Y": "NC_000024.10",
        "MT": "NC_012920.1",
    }
    _CHROM_TO_NC_GRCH37: dict[str, str] = {
        "1": "NC_000001.10",
        "2": "NC_000002.11",
        "3": "NC_000003.11",
        "4": "NC_000004.11",
        "5": "NC_000005.9",
        "6": "NC_000006.11",
        "7": "NC_000007.13",
        "8": "NC_000008.10",
        "9": "NC_000009.11",
        "10": "NC_000010.10",
        "11": "NC_000011.9",
        "12": "NC_000012.11",
        "13": "NC_000013.10",
        "14": "NC_000014.8",
        "15": "NC_000015.9",
        "16": "NC_000016.9",
        "17": "NC_000017.10",
        "18": "NC_000018.9",
        "19": "NC_000019.9",
        "20": "NC_000020.10",
        "21": "NC_000021.8",
        "22": "NC_000022.10",
        "X": "NC_000023.10",
        "Y": "NC_000024.9",
        "MT": "NC_012920.1",
    }

    def __init__(
        self,
        assembly: str = "GRCh38",
        cdot_json: str | Path | None = None,
    ):
        self.assembly = assembly

        self.hdp = (
            JSONDataProvider([str(cdot_json)])
            if cdot_json
            else RESTDataProvider()  # Uses API server at cdot.cc
        )
        self.am = hgvs.assemblymapper.AssemblyMapper(
            self.hdp,
            assembly_name=self.assembly,
            alt_aln_method="splign",
            replace_reference=True,
        )
        self.hp = hgvs.parser.Parser()
        self.CHROM_TO_NC = (
            self._CHROM_TO_NC_GRCH37
            if assembly in ("GRCh37", "hg19")
            else self._CHROM_TO_NC_GRCH38
        )
        self.hn = hgvs.normalizer.Normalizer(self.hdp)
        self.vm = hgvs.variantmapper.VariantMapper(self.hdp)

    @staticmethod
    def _build_hgvs_edit(ref: str, alt: str):
        """Baut das passende hgvs Edit-Objekt für SNV, Del, Ins, Delins."""
        ref_len = len(ref)
        alt_len = len(alt)

        if ref_len == 1 and alt_len == 1:
            # SNV
            return hgvs.edit.NARefAlt(ref=ref, alt=alt)
        elif ref_len == 0 and alt_len > 0:
            # pure Insertion
            return hgvs.edit.NARefAlt(ref=None, alt=alt)
        elif ref_len > 0 and alt_len == 0:
            # Deletion
            return hgvs.edit.NARefAlt(ref=ref, alt=None)
        else:
            # MNV / Delins
            return hgvs.edit.NARefAlt(ref=ref, alt=alt)

    @classmethod
    def _build_hgvs_input(cls, mutation: Mutation) -> str:
        """Builds a variant annotation to be mutalyzed."""
        prefix = f"{mutation.chromosome_ac}:g."
        # g_pos1 = genomic_pos + 1
        if mutation.mutation_type == "SNP":
            return (
                f"{prefix}{mutation.genomic_pos}{mutation.ref}>{mutation.alt}"
            )
        elif mutation.mutation_type.startswith("DEL"):
            return f"{prefix}{mutation.genomic_pos}_{mutation.genomic_pos + len(mutation.ref) - 1}del"  # noqa: E501
        elif mutation.mutation_type.startswith("INS"):
            return f"{prefix}{mutation.genomic_pos}_{mutation.genomic_pos + 1}ins{mutation.alt}"  # noqa: E501
        elif mutation.mutation_type in ("nonsense", "missense", "DELINS"):
            return f"{prefix}{mutation.genomic_pos}_{mutation.genomic_pos + len(mutation.ref) - 1}delins{mutation.alt}"  # noqa: E501
        return f"{prefix}{mutation.genomic_start-1}_{mutation.genomic_end}delins{mutation.genomic}"  # noqa: E501

    def build_g_variant(self, mutation: Mutation):
        """
        Build a genomic HGVS SequenceVariant from a Mutation object.
        Uses VCF-style (pos, ref, alt).
        """

        pos = mutation.genomic_pos  # 1-based
        ref = mutation.ref
        alt = mutation.alt

        # SNP
        if len(ref) == 1 and len(alt) == 1:
            loc = hgvs.location.SimplePosition(pos)
            edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)

            return hgvs.sequencevariant.SequenceVariant(
                ac=mutation.chromosome_ac,
                type="g",
                posedit=hgvs.posedit.PosEdit(loc, edit),
            )

        # DEL
        if len(ref) > 0 and alt == "":
            start = pos
            end = pos + len(ref) - 1

            loc = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(start),
                end=hgvs.location.SimplePosition(end),
            )
            edit = hgvs.edit.NARefAlt(ref=ref, alt=None)

            return hgvs.sequencevariant.SequenceVariant(
                ac=mutation.chromosome_ac,
                type="g",
                posedit=hgvs.posedit.PosEdit(loc, edit),
            )

        # INS
        if ref == "" and len(alt) > 0:
            # HGVS insertion: between bases → position is between pos-1 and pos
            start = pos - 1
            end = pos

            loc = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(start),
                end=hgvs.location.SimplePosition(end),
            )
            edit = hgvs.edit.NARefAlt(ref=None, alt=alt)

            return hgvs.sequencevariant.SequenceVariant(
                ac=mutation.chromosome_ac,
                type="g",
                posedit=hgvs.posedit.PosEdit(loc, edit),
            )

        # DELINS (catch-all: MNV, missense, nonsense, etc.)
        start = pos
        end = pos + len(ref) - 1

        loc = hgvs.location.Interval(
            start=hgvs.location.SimplePosition(start),
            end=hgvs.location.SimplePosition(end),
        )
        edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)

        return SequenceVariant(
            ac=mutation.chromosome_ac,
            type="g",
            posedit=hgvs.posedit.PosEdit(loc, edit),
        )

    ############################################################################
    # Create var_g
    ############################################################################

    @classmethod
    def spdi_to_var_g(cls, spdi: str) -> SequenceVariant:
        """
        Convert SPDI string to HGVS SequenceVariant (g.)
        """
        seq_id, pos, deleted, inserted = spdi.split(":")

        pos = int(pos) + 1  # SPDI (0-based) → HGVS (1-based)

        # SNP
        if len(deleted) == 1 and len(inserted) == 1:
            loc = hgvs.location.SimplePosition(pos)
            edit = hgvs.edit.NARefAlt(ref=deleted, alt=inserted)

        # deletion
        elif len(deleted) > 0 and inserted == "":
            start = pos
            end = pos + len(deleted) - 1

            loc = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(start),
                end=hgvs.location.SimplePosition(end),
            )
            edit = hgvs.edit.NARefAlt(ref=deleted, alt=None)

        # insertion
        elif deleted == "" and len(inserted) > 0:
            start = pos - 1
            end = pos

            loc = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(start),
                end=hgvs.location.SimplePosition(end),
            )
            edit = hgvs.edit.NARefAlt(ref=None, alt=inserted)

        # delins (catch-all)
        else:
            start = pos
            end = pos + len(deleted) - 1

            loc = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(start),
                end=hgvs.location.SimplePosition(end),
            )
            edit = hgvs.edit.NARefAlt(ref=deleted, alt=inserted)

        return hgvs.sequencevariant.SequenceVariant(
            ac=seq_id,
            type="g",
            posedit=hgvs.posedit.PosEdit(loc, edit),
        )

    @staticmethod
    def _build_g_variant(chrom_ac: str, genomic_pos: int, ref: str, alt: str):
        """
        Builds a SequenceVariant object for the genomic g. representation.
        genomic_pos is 0-based → will be internally converted to 1-based.
        """
        g_pos1 = genomic_pos  # + 1
        ref_len = len(ref)

        if ref_len <= 1:
            pos = hgvs.location.SimplePosition(g_pos1)
            interval = hgvs.location.Interval(start=pos, end=pos)
        else:
            start = hgvs.location.SimplePosition(g_pos1)
            end = hgvs.location.SimplePosition(g_pos1 + ref_len - 1)
            interval = hgvs.location.Interval(start=start, end=end)

        edit = HGVSMutationAnnotator._build_hgvs_edit(ref, alt)
        posedit = hgvs.posedit.PosEdit(pos=interval, edit=edit)

        return SequenceVariant(
            ac=chrom_ac,
            type="g",
            posedit=posedit,
        )

    ############################################################################
    # SequenceVariants
    ############################################################################

    def g_to_c(self, var_g: str, transcript_ac: str) -> SequenceVariant | None:
        try:
            var_c = self.am.g_to_c(var_g, transcript_ac)
        except (
            HGVSParseError,
            HGVSDataNotAvailableError,
        ) as e:
            print(
                f"Error parsing HGVS c variant: {e}\nCould not convert {var_g} to c. notation"  # noqa: E501
            )
            var_c = None
        return var_c

    def c_to_p(self, var_c: str | None, proc_ac: str) -> SequenceVariant | None:
        return ""  # this is currently not working, as our data provider does not know protein accession  # noqa_: E501
        if not var_c:
            return None
        try:
            print("var_c for P", var_c)
            var_p = self.am.c_to_p(var_c)
        except Exception as e:
            print(
                f"Error parsing HGVS p variant: {e}\nCould not convert {var_c} to p. notation"  # noqa: E501
            )
            var_p = None
        return var_p

    def c_to_g(self, var_c: str) -> SequenceVariant | None:
        try:
            var_g = self.am.c_to_g(var_c)
        except (HGVSParseError, HGVSDataNotAvailableError) as e:
            print(
                f"Error parsing HGVS c variant: {e}\nCould not convert {var_c} to g. notation"  # noqa: E501
            )
            var_g = None
        return var_g

    ############################################################################
    # HGVS
    ############################################################################

    def to_hgvs(self, var: SequenceVariant) -> str:
        if not var:
            return ""
        if var.type == "g":
            pass
        return self._normalize_hgvs(var)

    def hgvs_c_to_r(self, hgvs_c: str) -> str | None:
        if "-" in hgvs_c or "+" in hgvs_c:
            return None  # this is currently not working, as c notation can have outside positions # noqa_: E501  #noqa: E501
        if not isinstance(hgvs_c, str):
            print(
                f"Cannot convert {hgvs_c} of type {type(hgvs_c)} to r. notation"
            )  # noqa: E501
            return ""
        print("hgvs_c_to_r input: ", hgvs_c)
        return hgvs_c.replace("c.", "r.(") + ")"

    def _normalize_hgvs(
        self,
        var: SequenceVariant,
    ) -> str:
        if var is None:
            return ""
        try:
            normalized_hgvs = str(self.hn.normalize(var)) if var else str(var)
        except Exception:
            normalized_hgvs = str(var)
        return normalized_hgvs

    # def c_to_r(self, var_c: str) -> SequenceVariant:
    #     try:
    #         var_r = self.vm.c_to_r(var_c)
    #     except Exception as e:
    #         print(f"Could not convert {var_c} to r. notation: {e}")
    #         raise
    #     return var_r

    # def c_to_r(self, var_c: str) -> SequenceVariant:
    #     # # Setup (einmal global machen, nicht pro Call!)
    #     # # hp = hgvs.parser.Parser()
    #     # hdp = hgvs.dataproviders.uta.connect()
    #     # am = hgvs.assemblymapper.AssemblyMapper(
    #     #     hdp,
    #     #     assembly_name="GRCh38",
    #     #     alt_aln_method="splign",
    #     #     replace_reference=True,
    #     # )
    #     # Convert to r.
    #     try:
    #         var_r = self.am.c_to_r(var_c)
    #     except Exception as e:
    #         # sugggest = f"{chromosome_ac}({transcript_ac}):r.({str(var_c).split('.')[1]})"
    #         print(f"Could not convert {var_c} to r. notation. Error: {e}")
    #         raise
    #     return var_r

    ############################################################################
    # Annotate mutation
    ############################################################################

    def _annotate_variant(
        self,
        mutation: Mutation,
    ) -> tuple[str, str, str | None, str | None, int | None]:

        # var_g
        var_g = self._build_g_variant(
            mutation.chromosome_ac,
            mutation.genomic_pos,
            mutation.ref,
            mutation.alt,
        )
        try:
            var_g = self.am.normalize(var_g)
        except Exception:
            pass

        var_c = self.g_to_c(var_g, mutation.transcript_ac)
        var_p = self.c_to_p(var_c, mutation.protein_ac)

        hgvs_g = self.to_hgvs(var_g)
        hgvs_c = self.to_hgvs(var_c) if var_c else ""
        hgvs_p = self.to_hgvs(var_p) if var_p else ""
        hgvs_r = self.hgvs_c_to_r(hgvs_c) if hgvs_c else ""
        # # var_c
        # hgvs_c = ""
        # var_c = None
        # if mutation.cds_pos and mutation.transcript_ac:
        #     try:
        #         var_c = self.am.g_to_c(var_g, mutation.transcript_ac)
        #         hgvs_c = str(var_c)
        #     except hgvs.exceptions.HGVSError as e:
        #         # z.B. intronic / UTR variants die kein c. haben
        #         hgvs_c = f"?({e})"
        #         print(
        #             f"Could not convert {hgvs_g} to c. notation for transcript {mutation.transcript_ac}: {e}"
        #         )

        # # ── 3. c. → p. ────────────────────────────────────────────────────
        # hgvs_p = ""
        # if var_c is not None and mutation.region_type == "exon":
        #     try:
        #         var_p = self.c_to_p(var_c)
        #         hgvs_p = str(var_p)
        #     except hgvs.exceptions.HGVSError:
        #         hgvs_p = ""

        # hgvs_r = ""
        # if var_c is not None and mutation.region_type == "exon":
        #     try:
        #         var_r = self.c_to_r(var_c)
        #         hgvs_r = str(var_r)
        #     except hgvs.exceptions.HGVSError:
        #         hgvs_r = ""

        prot_pos = (mutation.cds_pos // 3 + 1) if mutation.cds_pos else None
        # print(hgvs_g, hgvs_c, hgvs_r, hgvs_p, prot_pos)
        return hgvs_g, hgvs_c, hgvs_r, hgvs_p, prot_pos

    def annotate(
        self,
        mutation: Mutation,
    ) -> Mutation:
        try:
            hgvs_g, hgvs_c, hgvs_r, hgvs_p, prot_pos = self._annotate_variant(
                mutation
            )
        except Exception as e:
            print(f"Error annotating mutation {mutation.seq_id}: {e}")
            print("mutation: ", mutation)
            raise
            hgvs_g, hgvs_c, hgvs_r, hgvs_p, prot_pos = "", "", "", "", None
        return replace(
            mutation,
            hgvs_g=hgvs_g,
            hgvs_c=hgvs_c,
            hgvs_r=hgvs_r,
            hgvs_p=hgvs_p,
            prot_pos=prot_pos,
        )


# class MaveMutationAnnotator(MutationAnnotator):

#     # NC_ accessions for GRCh38 (hg38) and GRCh37 (hg19)
#     _CHROM_TO_NC_GRCH38: dict[str, str] = {
#         "1": "NC_000001.11",
#         "2": "NC_000002.12",
#         "3": "NC_000003.12",
#         "4": "NC_000004.12",
#         "5": "NC_000005.10",
#         "6": "NC_000006.12",
#         "7": "NC_000007.14",
#         "8": "NC_000008.11",
#         "9": "NC_000009.12",
#         "10": "NC_000010.11",
#         "11": "NC_000011.10",
#         "12": "NC_000012.12",
#         "13": "NC_000013.11",
#         "14": "NC_000014.9",
#         "15": "NC_000015.10",
#         "16": "NC_000016.10",
#         "17": "NC_000017.11",
#         "18": "NC_000018.10",
#         "19": "NC_000019.10",
#         "20": "NC_000020.11",
#         "21": "NC_000021.9",
#         "22": "NC_000022.11",
#         "X": "NC_000023.11",
#         "Y": "NC_000024.10",
#         "MT": "NC_012920.1",
#     }
#     _CHROM_TO_NC_GRCH37: dict[str, str] = {
#         "1": "NC_000001.10",
#         "2": "NC_000002.11",
#         "3": "NC_000003.11",
#         "4": "NC_000004.11",
#         "5": "NC_000005.9",
#         "6": "NC_000006.11",
#         "7": "NC_000007.13",
#         "8": "NC_000008.10",
#         "9": "NC_000009.11",
#         "10": "NC_000010.10",
#         "11": "NC_000011.9",
#         "12": "NC_000012.11",
#         "13": "NC_000013.10",
#         "14": "NC_000014.8",
#         "15": "NC_000015.9",
#         "16": "NC_000016.9",
#         "17": "NC_000017.10",
#         "18": "NC_000018.9",
#         "19": "NC_000019.9",
#         "20": "NC_000020.10",
#         "21": "NC_000021.8",
#         "22": "NC_000022.10",
#         "X": "NC_000023.10",
#         "Y": "NC_000024.9",
#         "MT": "NC_012920.1",
#     }

#     # HGVS standard requires three-letter amino acid codes for p. notation
#     _AA_1TO3: dict[str, str] = {
#         "A": "Ala",
#         "C": "Cys",
#         "D": "Asp",
#         "E": "Glu",
#         "F": "Phe",
#         "G": "Gly",
#         "H": "His",
#         "I": "Ile",
#         "K": "Lys",
#         "L": "Leu",
#         "M": "Met",
#         "N": "Asn",
#         "P": "Pro",
#         "Q": "Gln",
#         "R": "Arg",
#         "S": "Ser",
#         "T": "Thr",
#         "V": "Val",
#         "W": "Trp",
#         "Y": "Tyr",
#         "*": "Ter",
#     }

#     def __init__(self, assembly: str = "GRCh38"):
#         # No UTA database connection required — mavehgvs is purely string-based.
#         self.CHROM_TO_NC = (
#             self._CHROM_TO_NC_GRCH37
#             if assembly in ("GRCh37", "hg19")
#             else self._CHROM_TO_NC_GRCH38
#         )

#     def _build_g_hgvs(
#         self, chrom_ac: str, genomic_pos: int, ref: str, alt: str
#     ) -> str:
#         """
#         Build and validate a genomic g. HGVS string.
#         genomic_pos is 0-based and is converted to 1-based internally.
#         Returns e.g. 'NC_000017.11:g.7674220A>T'
#         """
#         g1 = genomic_pos + 1
#         ref_len, alt_len = len(ref), len(alt)

#         if ref_len == 1 and alt_len == 1:
#             raw = f"g.{g1}{ref}>{alt}"
#         elif ref_len == 0:
#             # insertion between g1 and g1+1
#             raw = f"g.{g1}_{g1 + 1}ins{alt}"
#         elif alt_len == 0:
#             # deletion
#             raw = (
#                 f"g.{g1}del"
#                 if ref_len == 1
#                 else f"g.{g1}_{g1 + ref_len - 1}del"
#             )
#         else:
#             # MNV / delins
#             raw = (
#                 f"g.{g1}delins{alt}"
#                 if ref_len == 1
#                 else f"g.{g1}_{g1 + ref_len - 1}delins{alt}"
#             )

#         return f"{chrom_ac}:{HgvsVariant(raw)}"

#     def _annotate_variant(
#         self,
#         mutation: Mutation,
#     ) -> tuple[str, str, str, int | None]:
#         # use_cds = cds_pos is not None
#         chrom_ac = self.CHROM_TO_NC[str(mutation.chromosome)]
#         c_pos1 = (mutation.cds_pos + 1) if mutation.cds_pos else None  # type: ignore[operator]

#         # ── 1. g. ─────────────────────────────────────────────────────────
#         hgvs_g = self._build_g_hgvs(
#             chrom_ac, mutation.genomic_pos, mutation.ref, mutation.alt
#         )

#         # ── 2. c. ─────────────────────────────────────────────────────────
#         hgvs_c = ""
#         if mutation.cds_pos and mutation.genomic_id:
#             if not mutation.ref and mutation.alt:
#                 # insertion: 1-based adjacent positions
#                 raw_c = f"c.{c_pos1}_{c_pos1 + 1}ins{mutation.alt}"  # type: ignore[operator]
#             elif mutation.ref and not mutation.alt:
#                 if len(mutation.ref) == 1:
#                     raw_c = f"c.{c_pos1}del"
#                 else:
#                     end = c_pos1 + len(mutation.ref) - 1  # type: ignore[operator]
#                     raw_c = f"c.{c_pos1}_{end}del"
#             elif mutation.ref_codon and mutation.alt_codon:
#                 diffs = [
#                     i
#                     for i in range(3)
#                     if mutation.ref_codon[i] != mutation.alt_codon[i]
#                 ]
#                 if len(diffs) == 1:
#                     pos_nt = mutation.cds_pos + diffs[0] + 1  # type: ignore[operator]
#                     raw_c = (
#                         f"c.{pos_nt}"
#                         f"{mutation.ref_codon[diffs[0]]}>{mutation.alt_codon[diffs[0]]}"
#                     )
#                 else:
#                     start_nt = mutation.cds_pos + diffs[0] + 1  # type: ignore[operator]
#                     end_nt = mutation.cds_pos + diffs[-1] + 1  # type: ignore[operator]
#                     alt_sub = mutation.alt_codon[diffs[0] : diffs[-1] + 1]
#                     raw_c = f"c.{start_nt}_{end_nt}delins{alt_sub}"
#             else:
#                 raw_c = f"c.{c_pos1}{mutation.ref}>{mutation.alt}"
#             hgvs_c = f"{mutation.genomic_id}:{HgvsVariant(raw_c)}"

#         # ── 3. p. ─────────────────────────────────────────────────────────
#         hgvs_p = ""
#         if (
#             mutation.region_type == "exon"
#             and mutation.ref_aa
#             and mutation.alt_aa
#             and mutation.cds_pos
#         ):
#             codon_index = mutation.cds_pos // 3 + 1  # type: ignore[operator]
#             ref_aa_3 = self._AA_1TO3.get(mutation.ref_aa, mutation.ref_aa)
#             alt_aa_3 = self._AA_1TO3.get(mutation.alt_aa, mutation.alt_aa)
#             raw_p = f"p.{ref_aa_3}{codon_index}{alt_aa_3}"
#             hgvs_p = f"{mutation.genomic_id}:{HgvsVariant(raw_p)}"

#         prot_pos = (mutation.cds_pos // 3 + 1) if mutation.cds_pos else None  # type: ignore[operator]
#         return hgvs_g, hgvs_c, hgvs_p, prot_pos

#     def annotate(
#         self,
#         mutation: Mutation,
#     ) -> Mutation:
#         hgvs_g, hgvs_c, hgvs_p, prot_pos = self._annotate_variant(mutation)
#         return replace(
#             mutation,
#             hgvs_g=hgvs_g,
#             hgvs_c=hgvs_c,
#             hgvs_p=hgvs_p,
#             prot_pos=prot_pos,
#         )


# class ManualMutationAnnotator(MutationAnnotator):

#     def __generate_hgvs_g(self, mutation: Mutation) -> str:
#         # ============================================================
#         # 🧬 GENOMIC HGVS (g.)
#         # ============================================================
#         g_start = mutation.genomic_pos  # this should be 1-based
#         if not mutation.ref and mutation.alt:
#             # insertion
#             hgvs_g = f"{mutation.genomic_id}:g.{g_start}_{g_start+1}ins{mutation.alt}"  # noqa: E501

#         elif mutation.ref and not mutation.alt:
#             # deletion
#             if len(mutation.ref) == 1:
#                 hgvs_g = f"{mutation.genomic_id}:g.{g_start}del"
#             else:
#                 end = g_start + len(mutation.ref) - 1
#                 hgvs_g = f"{mutation.genomic_id}:g.{g_start}_{end}del"

#         elif (
#             mutation.ref
#             and mutation.alt
#             and len(mutation.ref) == len(mutation.alt) == 1
#         ):
#             # SNP
#             hgvs_g = f"{mutation.genomic_id}:g.{g_start}{mutation.ref}>{mutation.alt}"  # noqa: E501

#         else:
#             # complex replacement
#             end = g_start + len(mutation.ref) - 1
#             hgvs_g = (
#                 f"{mutation.genomic_id}:g.{g_start}_{end}delins{mutation.alt}"
#             )
#         return hgvs_g

#     def __generate_hgvs_c(self, mutation: Mutation) -> str:
#         if mutation.cds_pos:
#             c_start = (
#                 mutation.cds_pos + 1
#             )  # this may be 0-based? should be 1-based
#             # ---- insertion ----
#             if not mutation.ref and mutation.alt:
#                 hgvs_c = f"{mutation.seq_id}:c.{c_start}_{c_start+1}ins{mutation.alt}"  # noqa: E501

#             # ---- deletion ----
#             elif mutation.ref and not mutation.alt:
#                 if len(mutation.ref) == 1:
#                     hgvs_c = f"{mutation.seq_id}:c.{c_start}del"
#                 else:
#                     end = c_start + len(mutation.ref) - 1
#                     hgvs_c = f"{mutation.seq_id}:c.{c_start}_{end}del"

#             # ---- codon-aware substitution ----
#             elif mutation.ref_codon and mutation.alt_codon:

#                 diffs = [
#                     i
#                     for i in range(3)
#                     if mutation.ref_codon[i] != mutation.alt_codon[i]
#                 ]

#                 if len(diffs) == 1:
#                     pos_nt = mutation.cds_pos + diffs[0] + 1
#                     hgvs_c = (
#                         f"{mutation.seq_id}:c.{pos_nt}"
#                         f"{mutation.ref_codon[diffs[0]]}>{mutation.alt_codon[diffs[0]]}"  # noqa: E501
#                     )

#                 else:
#                     start_nt = mutation.cds_pos + diffs[0] + 1
#                     end_nt = mutation.cds_pos + diffs[-1] + 1
#                     alt_sub = mutation.alt_codon[diffs[0] : diffs[-1] + 1]

#                     hgvs_c = (
#                         f"{mutation.seq_id}:c.{start_nt}_{end_nt}"
#                         f"delins{alt_sub}"
#                     )

#             # ---- fallback SNP ----
#             elif (
#                 mutation.ref
#                 and mutation.alt
#                 and len(mutation.ref) == len(mutation.alt) == 1
#             ):
#                 hgvs_c = f"{mutation.seq_id}:c.{c_start}{mutation.ref}>{mutation.alt}"  # noqa: E501

#             else:
#                 end = c_start + len(mutation.ref) - 1
#                 hgvs_c = (
#                     f"{mutation.seq_id}:c.{c_start}_{end}delins{mutation.alt}"
#                 )

#         else:
#             hgvs_c = ""
#         return hgvs_c

#     def __generate_hgvs_p(self, mutation: Mutation) -> tuple[str, int | None]:
#         if (
#             mutation.region_type == "exon"
#             and mutation.ref_aa
#             and mutation.alt_aa
#             and mutation.cds_pos
#         ):

#             prot_pos = mutation.cds_pos // 3 + 1

#             if mutation.alt_aa == "*":
#                 hgvs_p = f"{mutation.seq_id}:p.{mutation.ref_aa}{prot_pos}Ter"

#             elif mutation.ref_aa == mutation.alt_aa:
#                 hgvs_p = f"{mutation.seq_id}:p.{mutation.ref_aa}{prot_pos}="

#             else:
#                 hgvs_p = f"{mutation.seq_id}:p.{mutation.ref_aa}{prot_pos}{mutation.alt_aa}"  # noqa: E501

#         else:
#             hgvs_p = ""
#             prot_pos = None
#         return hgvs_p, prot_pos

#     def _annotate_variant(
#         self,
#         mutation: Mutation,
#     ) -> tuple[str, str, str, int | None]:
#         hgvs_g = self.__generate_hgvs_g(mutation)
#         hgvs_c = self.__generate_hgvs_c(mutation)
#         hgvs_p, prot_pos = self.__generate_hgvs_p(mutation)

#         return hgvs_g, hgvs_c, hgvs_p, prot_pos

#     def annotate(
#         self,
#         mutation: Mutation,
#     ) -> Mutation:

#         hgvs_g, hgvs_c, hgvs_p, prot_pos = self._annotate_variant(mutation)
#         return replace(
#             mutation,
#             hgvs_g=hgvs_g,
#             hgvs_c=hgvs_c,
#             hgvs_p=hgvs_p,
#             prot_pos=prot_pos,
#         )


# class Mutalyzer(MutationAnnotator):

#     BASE_URL = "https://mutalyzer.nl/api/normalize/"

#     @classmethod
#     def _normalize(cls, hgvs: str) -> dict:
#         encoded = quote(hgvs, safe="")
#         url = f"{cls.BASE_URL}{encoded}"

#         r = requests.get(
#             url,
#             params={"only_variants": "false"},
#             headers={"accept": "application/json"},
#             timeout=30,
#         )
#         if r.status_code != 200:
#             print(hgvs, url)
#             raise ValueError(f"Mutalyzer error {r.status_code}: {r.text[:200]}")

#         if not r.text.strip():
#             raise ValueError(f"Empty Mutalyzer response for {hgvs}")

#         data = r.json()
#         # print("json data", json.dumps(data, indent=2))
#         # fallback-safe extraction
#         return cls._extract_fields(data)

#     @staticmethod
#     def _extract_fields(data: dict) -> dict:
#         """
#         Robust Mutalyzer normalization extractor.
#         Works with /api/normalize response structure.
#         """

#         # ----------------------------
#         # 1. canonical description (MOST IMPORTANT)
#         # ----------------------------
#         canonical = data.get("normalized_description")
#         transcript_ac = (
#             data.get("normalized_model", {})
#             .get("reference", {})
#             .get("selector", {})
#             .get("id", None)
#         )
#         chrom_ac = (
#             data.get("normalized_model", {})
#             .get("reference", {})
#             .get("reference", {})
#             .get("id", None)
#         )

#         if not canonical:
#             # fallback: sometimes only corrected_description exists
#             canonical = data.get("corrected_description")

#         # ----------------------------
#         # 2. genomic equivalents
#         # ----------------------------
#         eq = data.get("equivalent_descriptions", {})

#         g_list = eq.get("g", [])
#         c_list = eq.get("c", [])

#         hgvs_g = g_list[0]["description"] if g_list else None

#         # try to pick matching transcript (prefer MANE / NM_000546.6 if present)
#         hgvs_c = None
#         for item in c_list:
#             desc = item.get("description")
#             if "NM_000546.6" in desc:
#                 hgvs_c = desc
#                 break

#         # fallback: first c
#         if hgvs_c is None and c_list:
#             hgvs_c = c_list[0]["description"]

#         # ----------------------------
#         # 3. protein + rna (optional, often missing)
#         # ----------------------------
#         protein = data.get("protein", {})
#         rna = data.get("rna", {})

#         hgvs_p = None
#         hgvs_r = None

#         # Mutalyzer does NOT reliably return p/r normalized strings
#         # so we only extract if explicitly present (future-proof)
#         protein_ac = None
#         if isinstance(protein, dict):
#             hgvs_p = protein.get("description") or protein.get("hgvs_p")
#             if hgvs_p:
#                 m = re.search(r"\((NP_[0-9]+\.[0-9]+)\)", hgvs_p)
#                 protein_ac = m.group(1) if m else None
#             else:
#                 print("protein was", protein, hgvs_p)

#         if isinstance(rna, dict):
#             hgvs_r = rna.get("description") or rna.get("hgvs_r")

#         # ----------------------------
#         # 4. fallback consistency rule
#         # ----------------------------
#         if canonical and not hgvs_c:
#             hgvs_c = canonical
#         protein = "notImplemented"
#         return {
#             "hgvs_c": hgvs_c,
#             "hgvs_g": hgvs_g,
#             "hgvs_p": hgvs_p,
#             "hgvs_r": hgvs_r,
#             "reference_ac": chrom_ac,
#             "transcript_ac": transcript_ac,
#             "protein_ac": protein_ac,
#             "raw": data,
#         }

#     @classmethod
#     def normalize_hgvs_c(cls, hgvs_c: str) -> dict:
#         if "c." not in hgvs_c and "NM_" not in hgvs_c:
#             raise ValueError("Invalid hgvs_c")
#         return cls._normalize(hgvs_c)

#     @classmethod
#     def normalize_hgvs_g(cls, hgvs_g: str) -> dict:
#         if ":g." not in hgvs_g:
#             raise ValueError("Invalid hgvs_g")
#         return cls._normalize(hgvs_g)

#     @classmethod
#     def normalize_hgvs_r(cls, hgvs_r: str) -> dict:
#         if ":r." not in hgvs_r and not hgvs_r.startswith("r."):
#             raise ValueError("Invalid hgvs_r")
#         return cls._normalize(hgvs_r)

#     @classmethod
#     def normalize_hgvs_p(cls, hgvs_p: str) -> dict:
#         try:
#             return cls._normalize(hgvs_p)
#         except Exception:
#             return {
#                 "hgvs_p": hgvs_p,
#                 "hgvs_c": None,
#                 "hgvs_g": None,
#                 "hgvs_r": None,
#             }

#     @classmethod
#     def normalize(cls, hgvs: str) -> dict:
#         return cls._normalize(hgvs)

#     @classmethod
#     def _normalize_mutalyzer_to_hgvs(
#         cls, hgvs: str | None, id: str | None = None
#     ) -> str | None:
#         ret = f"{id}:{hgvs.split(':', 1)[-1]}" if hgvs and id else hgvs
#         if ret and "[" in ret:
#             ret = cls.fix_mutalyzer_repeat(ret)
#         return ret

#     @classmethod
#     def fix_mutalyzer_repeat(cls, hgvs: str) -> str:
#         """
#         Converts Mutalyzer repeat notation (A[3]) into HGVS-compatible ins
#         notation.
#         Works for c. and g. coordinates.
#         """
#         pattern = re.compile(r"([cg]\.)(\d+)([ACGT])\[(\d+)\]")

#         def repl(match):
#             coord_type = match.group(1)  # c. or g.
#             pos = int(match.group(2))
#             base = match.group(3)
#             n = int(match.group(4))

#             if n <= 1:
#                 return match.group(0)

#             inserted = base * (n - 1)
#             return f"{coord_type}{pos}_{pos+1}ins{inserted}"

#         return pattern.sub(repl, hgvs)

#     @classmethod
#     def mutalyzer_to_hgvs(
#         cls,
#         mutalyzer_out: dict[str, str],
#     ) -> dict[str, str | None]:
#         """
#         Extracts the most important fields from Mutalyzer output for the
#         mutation instance.
#         Returns a dict with hgvs_c, hgvs_g, hgvs_p, hgvs_r, chrom_ac,
#         transcript_ac, protein_ac.
#         """
#         try:
#             data = {}
#             hgvs_c = mutalyzer_out.get("hgvs_c")
#             hgvs_g = mutalyzer_out.get("hgvs_g")
#             hgvs_p = mutalyzer_out.get("hgvs_p")
#             hgvs_r = mutalyzer_out.get("hgvs_r")
#             chrom_ac = mutalyzer_out.get("chrom_ac")
#             transcript_ac = mutalyzer_out.get("transcript_ac")
#             protein_ac = mutalyzer_out.get("protein_ac")
#             data["hgvs_c"] = cls._normalize_mutalyzer_to_hgvs(
#                 hgvs_c, transcript_ac
#             )
#             data["hgvs_g"] = cls._normalize_mutalyzer_to_hgvs(hgvs_g, chrom_ac)
#             data["hgvs_p"] = cls._normalize_mutalyzer_to_hgvs(
#                 hgvs_p, protein_ac
#             )
#             data["hgvs_r"] = cls._normalize_mutalyzer_to_hgvs(
#                 hgvs_r, transcript_ac
#             )
#             return data
#         except Exception:
#             print(mutalyzer_out["raw"])
#             raise

#     @classmethod
#     def mutalyzer_to_hgvs_compliant(
#         cls, hgvs: str, accession: str
#     ) -> str | None:
#         return cls._normalize_mutalyzer_to_hgvs(hgvs, accession)

#     @classmethod
#     def _build_hgvs_input(cls, mutation: Mutation) -> str:
#         """Builds a variant annotation to be mutalyzed."""
#         prefix = f"{mutation.chromosome_ac}:g."
#         # g_pos1 = genomic_pos + 1
#         if mutation.mutation_type == "SNP":
#             return (
#                 f"{prefix}{mutation.genomic_pos}{mutation.ref}>{mutation.alt}"
#             )
#         elif mutation.mutation_type.startswith("DEL"):
#             return f"{prefix}{mutation.genomic_pos}_{mutation.genomic_pos + len(mutation.ref) - 1}del"  # noqa: E501
#         elif mutation.mutation_type.startswith("INS"):
#             return f"{prefix}{mutation.genomic_pos}_{mutation.genomic_pos + 1}ins{mutation.alt}"  # noqa: E501
#         elif mutation.mutation_type in ("nonsense", "missense", "DELINS"):
#             return f"{prefix}{mutation.genomic_pos}_{mutation.genomic_pos + len(mutation.ref) - 1}delins{mutation.alt}"  # noqa: E501
#         return f"{prefix}{mutation.genomic_start-1}_{mutation.genomic_end}delins{mutation.genomic}"  # noqa: E501

#     def annotate(
#         self,
#         mutation: Mutation,
#     ) -> Mutation:
#         # var_g = HGVSMutationAnnotator._build_g_variant(
#         #     mutation.chromosome_ac,
#         #     mutation.genomic_pos,
#         #     mutation.ref,
#         #     mutation.alt,
#         # )
#         # try:
#         mutalizer_input = self._build_hgvs_input(mutation)
#         # try:
#         #     normalized = self.normalize(mutalizer_input)
#         # except:
#         #     print(mutalizer_input)
#         #     print("MUTATION", mutation)
#         #     print(mutation.coding, mutation.genomic)
#         #     print("input for normalizer", mutalizer_input)
#         #     raise
#         return replace(
#             mutation,
#             hgvs_g=mutalizer_input,  # normalized["hgvs_g"],
#             # hgvs_c=normalized["hgvs_c"],
#             # hgvs_p=normalized["hgvs_p"],
#             # prot_pos=normalized.get("prot_pos"),
#         )
