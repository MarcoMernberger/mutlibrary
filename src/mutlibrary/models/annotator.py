from abc import ABC, abstractmethod
from dataclasses import dataclass, replace

import hgvs.assemblymapper  # type: ignore[import]
import hgvs.dataproviders.uta  # type: ignore[import]
import hgvs.edit  # type: ignore[import]
import hgvs.exceptions  # type: ignore[import]
import hgvs.location  # type: ignore[import]
import hgvs.parser  # type: ignore[import]
import hgvs.posedit  # type: ignore[import]
import hgvs.sequencevariant  # type: ignore[import]
from mavehgvs.variant import Variant as HgvsVariant  # type: ignore[import]


@dataclass
class Mutation:
    seq_id: str  # a name for the sequence
    genomic_id: str  # the genomic_id for normalized hgvs, e.g. NC_000017.11
    chromosome: str
    region_type: str  # intron or exon
    mutation_pos: int  #   the position in the input_sequence
    genomic_pos: int  #   the genomic position of the mutation (1-based)
    cds_pos: int | None  # the position in the CDS (0-based), None if non-coding
    ref: str  # the reference allele (nucleotide(s) in the input sequence), as in vcf
    alt: str  # the alternate allele (nucleotide(s) in the input sequence), as in vcf
    mutation_type: str  # e.g. "SNV", "Del", "Ins", "MNV", "Delins"
    seq: str  # the mutated sequence
    ref_codon: (
        str | None
    )  # the reference codon (3 nucleotides in the input sequence), None if non-coding
    alt_codon: (
        str | None
    )  # the alternate codon (3 nucleotides in the input sequence), None if non-coding
    ref_aa: str | None  # the reference amino acid, None if non-coding
    alt_aa: str | None  # the alternate amino acid, None if non-coding
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


class MutationAnnotator(ABC):

    def __init__(self):
        pass

    @abstractmethod
    def annotate(
        self,
        seq_id: str,
        genomic_id: str,
        chrom: str,
        mutation_pos: int,
        genomic_pos: int,
        region_type: str,
        ref: str,
        alt: str,
        mut_type: str,
        seq: str,
        cds_pos: int | None = None,
        ref_codon: str | None = None,
        alt_codon: str | None = None,
        ref_aa: str | None = None,
        alt_aa: str | None = None,
    ) -> Mutation:
        raise NotImplementedError()


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

    def __init__(self, assembly: str = "GRCh38"):
        # Connect to UTA — uses UTA_DB_URL env var if set, otherwise the
        # public UTA REST endpoint (requires network access).
        hdp = hgvs.dataproviders.uta.connect()
        self.am = hgvs.assemblymapper.AssemblyMapper(
            hdp,
            assembly_name=assembly,
            alt_aln_method="splign",
            replace_reference=True,
        )
        self.hp = hgvs.parser.Parser()
        self.CHROM_TO_NC = (
            self._CHROM_TO_NC_GRCH37
            if assembly in ("GRCh37", "hg19")
            else self._CHROM_TO_NC_GRCH38
        )

    @staticmethod
    def _build_hgvs_edit(ref: str, alt: str):
        """Baut das passende hgvs Edit-Objekt für SNV, Del, Ins, Delins."""
        ref_len = len(ref)
        alt_len = len(alt)

        if ref_len == 1 and alt_len == 1:
            # SNV
            return hgvs.edit.NARefAlt(ref=ref, alt=alt)
        elif ref_len == 0 and alt_len > 0:
            # pure Insertion — sollte nicht vorkommen wenn ref die Anchor-Base enthält
            return hgvs.edit.NARefAlt(ref=None, alt=alt)
        elif ref_len > 0 and alt_len == 0:
            # Deletion
            return hgvs.edit.NARefAlt(ref=ref, alt=None)
        else:
            # MNV / Delins
            return hgvs.edit.NARefAlt(ref=ref, alt=alt)

    @staticmethod
    def _build_g_variant(chrom_ac: str, genomic_pos: int, ref: str, alt: str):
        """
        Baut ein SequenceVariant-Objekt für die genomische g.-Darstellung.
        genomic_pos ist 0-based → wird intern zu 1-based konvertiert.
        """
        g_pos1 = genomic_pos + 1
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

        return hgvs.sequencevariant.SequenceVariant(
            ac=chrom_ac,
            type="g",
            posedit=posedit,
        )

    def _annotate_variant(
        self,
        info: Mutation,
    ) -> tuple[str, str, str, int | None]:

        # use_cds = cds_pos is not None
        chrom_ac = self.CHROM_TO_NC[info.chromosome]

        # ── 1. Genomisches Variant-Objekt bauen ───────────────────────────
        var_g = self._build_g_variant(
            chrom_ac, info.genomic_pos, info.ref, info.alt
        )

        # hgvs normalisiert und korrigiert die Referenzbase automatisch
        try:
            var_g = self.am.normalize(var_g)
        except Exception:
            pass  # bei synthetischen Sequenzen kann Normalisierung fehlschlagen

        hgvs_g = str(var_g)

        # ── 2. g. → c. ────────────────────────────────────────────────────
        hgvs_c = ""
        var_c = None
        if info.cds_pos and info.genomic_id:
            try:
                var_c = self.am.g_to_c(var_g, info.genomic_id)
                hgvs_c = str(var_c)
            except hgvs.exceptions.HGVSError as e:
                # z.B. intronic / UTR variants die kein c. haben
                hgvs_c = f"?({e})"

        # ── 3. c. → p. ────────────────────────────────────────────────────
        hgvs_p = ""
        if var_c is not None and info.region_type == "exon":
            try:
                var_p = self.am.c_to_p(var_c)
                hgvs_p = str(var_p)
            except hgvs.exceptions.HGVSError:
                hgvs_p = ""

        prot_pos = (info.cds_pos // 3 + 1) if info.cds_pos else None

        return hgvs_g, hgvs_c, hgvs_p, prot_pos

    def annotate(
        self,
        info: Mutation,
    ) -> Mutation:

        hgvs_g, hgvs_c, hgvs_p, prot_pos = self._annotate_variant(info)
        return replace(
            info, hgvs_g=hgvs_g, hgvs_c=hgvs_c, hgvs_p=hgvs_p, prot_pos=prot_pos
        )


class MaveMutationAnnotator(MutationAnnotator):

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

    def __init__(self, assembly: str = "GRCh38"):
        # No UTA database connection required — mavehgvs is purely string-based.
        self.CHROM_TO_NC = (
            self._CHROM_TO_NC_GRCH37
            if assembly in ("GRCh37", "hg19")
            else self._CHROM_TO_NC_GRCH38
        )

    def _build_g_hgvs(
        self, chrom_ac: str, genomic_pos: int, ref: str, alt: str
    ) -> str:
        """
        Build and validate a genomic g. HGVS string.
        genomic_pos is 0-based and is converted to 1-based internally.
        Returns e.g. 'NC_000017.11:g.7674220A>T'
        """
        g1 = genomic_pos + 1
        ref_len, alt_len = len(ref), len(alt)

        if ref_len == 1 and alt_len == 1:
            raw = f"g.{g1}{ref}>{alt}"
        elif ref_len == 0:
            # insertion between g1 and g1+1
            raw = f"g.{g1}_{g1 + 1}ins{alt}"
        elif alt_len == 0:
            # deletion
            raw = (
                f"g.{g1}del"
                if ref_len == 1
                else f"g.{g1}_{g1 + ref_len - 1}del"
            )
        else:
            # MNV / delins
            raw = (
                f"g.{g1}delins{alt}"
                if ref_len == 1
                else f"g.{g1}_{g1 + ref_len - 1}delins{alt}"
            )

        return f"{chrom_ac}:{HgvsVariant(raw)}"

    def _annotate_variant(
        self,
        info: Mutation,
    ) -> tuple[str, str, str, int | None]:
        # use_cds = cds_pos is not None
        chrom_ac = self.CHROM_TO_NC[str(info.chromosome)]
        c_pos1 = (info.cds_pos + 1) if info.cds_pos else None  # type: ignore[operator]

        # ── 1. g. ─────────────────────────────────────────────────────────
        hgvs_g = self._build_g_hgvs(
            chrom_ac, info.genomic_pos, info.ref, info.alt
        )

        # ── 2. c. ─────────────────────────────────────────────────────────
        hgvs_c = ""
        if info.cds_pos and info.genomic_id:
            if not info.ref and info.alt:
                # insertion: 1-based adjacent positions
                raw_c = f"c.{c_pos1}_{c_pos1 + 1}ins{info.alt}"  # type: ignore[operator]
            elif info.ref and not info.alt:
                if len(info.ref) == 1:
                    raw_c = f"c.{c_pos1}del"
                else:
                    end = c_pos1 + len(info.ref) - 1  # type: ignore[operator]
                    raw_c = f"c.{c_pos1}_{end}del"
            elif info.ref_codon and info.alt_codon:
                diffs = [
                    i
                    for i in range(3)
                    if info.ref_codon[i] != info.alt_codon[i]
                ]
                if len(diffs) == 1:
                    pos_nt = info.cds_pos + diffs[0] + 1  # type: ignore[operator]
                    raw_c = (
                        f"c.{pos_nt}"
                        f"{info.ref_codon[diffs[0]]}>{info.alt_codon[diffs[0]]}"
                    )
                else:
                    start_nt = info.cds_pos + diffs[0] + 1  # type: ignore[operator]
                    end_nt = info.cds_pos + diffs[-1] + 1  # type: ignore[operator]
                    alt_sub = info.alt_codon[diffs[0] : diffs[-1] + 1]
                    raw_c = f"c.{start_nt}_{end_nt}delins{alt_sub}"
            else:
                raw_c = f"c.{c_pos1}{info.ref}>{info.alt}"
            hgvs_c = f"{info.genomic_id}:{HgvsVariant(raw_c)}"

        # ── 3. p. ─────────────────────────────────────────────────────────
        hgvs_p = ""
        if (
            info.region_type == "exon"
            and info.ref_aa
            and info.alt_aa
            and info.cds_pos
        ):
            codon_index = info.cds_pos // 3 + 1  # type: ignore[operator]
            ref_aa_3 = self._AA_1TO3.get(info.ref_aa, info.ref_aa)
            alt_aa_3 = self._AA_1TO3.get(info.alt_aa, info.alt_aa)
            raw_p = f"p.{ref_aa_3}{codon_index}{alt_aa_3}"
            hgvs_p = f"{info.genomic_id}:{HgvsVariant(raw_p)}"

        prot_pos = (info.cds_pos // 3 + 1) if info.cds_pos else None  # type: ignore[operator]
        return hgvs_g, hgvs_c, hgvs_p, prot_pos

    def annotate(
        self,
        info: Mutation,
    ) -> Mutation:
        hgvs_g, hgvs_c, hgvs_p, prot_pos = self._annotate_variant(info)
        return replace(
            info, hgvs_g=hgvs_g, hgvs_c=hgvs_c, hgvs_p=hgvs_p, prot_pos=prot_pos
        )


class ManualMutationAnnotator(MutationAnnotator):

    def __generate_hgvs_g(self, info: Mutation) -> str:
        # ============================================================
        # 🧬 GENOMIC HGVS (g.)
        # ============================================================
        g_start = info.genomic_pos  # this should be 1-based
        if not info.ref and info.alt:
            # insertion
            hgvs_g = f"{info.genomic_id}:g.{g_start}_{g_start+1}ins{info.alt}"

        elif info.ref and not info.alt:
            # deletion
            if len(info.ref) == 1:
                hgvs_g = f"{info.genomic_id}:g.{g_start}del"
            else:
                end = g_start + len(info.ref) - 1
                hgvs_g = f"{info.genomic_id}:g.{g_start}_{end}del"

        elif info.ref and info.alt and len(info.ref) == len(info.alt) == 1:
            # SNP
            hgvs_g = f"{info.genomic_id}:g.{g_start}{info.ref}>{info.alt}"

        else:
            # complex replacement
            end = g_start + len(info.ref) - 1
            hgvs_g = f"{info.genomic_id}:g.{g_start}_{end}delins{info.alt}"
        return hgvs_g

    def __generate_hgvs_c(self, info: Mutation) -> str:
        if info.cds_pos:
            c_start = info.cds_pos + 1  # this may be 0-based? should be 1-based
            # ---- insertion ----
            if not info.ref and info.alt:
                hgvs_c = f"{info.seq_id}:c.{c_start}_{c_start+1}ins{info.alt}"

            # ---- deletion ----
            elif info.ref and not info.alt:
                if len(info.ref) == 1:
                    hgvs_c = f"{info.seq_id}:c.{c_start}del"
                else:
                    end = c_start + len(info.ref) - 1
                    hgvs_c = f"{info.seq_id}:c.{c_start}_{end}del"

            # ---- codon-aware substitution ----
            elif info.ref_codon and info.alt_codon:

                diffs = [
                    i
                    for i in range(3)
                    if info.ref_codon[i] != info.alt_codon[i]
                ]

                if len(diffs) == 1:
                    pos_nt = info.cds_pos + diffs[0] + 1
                    hgvs_c = (
                        f"{info.seq_id}:c.{pos_nt}"
                        f"{info.ref_codon[diffs[0]]}>{info.alt_codon[diffs[0]]}"
                    )

                else:
                    start_nt = info.cds_pos + diffs[0] + 1
                    end_nt = info.cds_pos + diffs[-1] + 1
                    alt_sub = info.alt_codon[diffs[0] : diffs[-1] + 1]

                    hgvs_c = (
                        f"{info.seq_id}:c.{start_nt}_{end_nt}"
                        f"delins{alt_sub}"
                    )

            # ---- fallback SNP ----
            elif info.ref and info.alt and len(info.ref) == len(info.alt) == 1:
                hgvs_c = f"{info.seq_id}:c.{c_start}{info.ref}>{info.alt}"

            else:
                end = c_start + len(info.ref) - 1
                hgvs_c = f"{info.seq_id}:c.{c_start}_{end}delins{info.alt}"

        else:
            hgvs_c = ""
        return hgvs_c

    def __generate_hgvs_p(self, info: Mutation) -> tuple[str, int | None]:
        if (
            info.region_type == "exon"
            and info.ref_aa
            and info.alt_aa
            and info.cds_pos
        ):

            prot_pos = info.cds_pos // 3 + 1

            if info.alt_aa == "*":
                hgvs_p = f"{info.seq_id}:p.{info.ref_aa}{prot_pos}Ter"

            elif info.ref_aa == info.alt_aa:
                hgvs_p = f"{info.seq_id}:p.{info.ref_aa}{prot_pos}="

            else:
                hgvs_p = f"{info.seq_id}:p.{info.ref_aa}{prot_pos}{info.alt_aa}"

        else:
            hgvs_p = ""
            prot_pos = None
        return hgvs_p, prot_pos

    def _annotate_variant(
        self,
        info: Mutation,
    ) -> tuple[str, str, str, int | None]:
        hgvs_g = self.__generate_hgvs_g(info)
        hgvs_c = self.__generate_hgvs_c(info)
        hgvs_p, prot_pos = self.__generate_hgvs_p(info)

        return hgvs_g, hgvs_c, hgvs_p, prot_pos

    def annotate(
        self,
        info: Mutation,
    ) -> Mutation:

        hgvs_g, hgvs_c, hgvs_p, prot_pos = self._annotate_variant(info)
        return replace(
            info, hgvs_g=hgvs_g, hgvs_c=hgvs_c, hgvs_p=hgvs_p, prot_pos=prot_pos
        )
