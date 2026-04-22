from __future__ import annotations

import re
from pathlib import Path
from typing import Literal

import hgvs.assemblymapper
import hgvs.edit
import hgvs.location
import hgvs.parser
import hgvs.posedit
import hgvs.sequencevariant
from Bio import SeqIO
from Bio.Seq import Seq
from biocommons.seqrepo import SeqRepo
from biocommons.seqrepo.dataproxy import SeqRepoRESTDataProxy
from cdot.hgvs.dataproviders import JSONDataProvider, RESTDataProvider
from hgvs.sequencevariant import SequenceVariant
from mmalignments.utils.utils import reverse_complement

from .annotator import Mutation


def _mut_type(ref: str | None, alt: str | None) -> str:
    if not ref:
        return "Ins"
    if not alt:
        return "Del"
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    return "Delins"


def _apply_edit(seq: str, start_0: int, ref: str, alt: str) -> str:
    """
    Appliziert ref→alt auf seq ab Position start_0 (0-based).
    ref darf leer sein (Insertion), alt darf leer sein (Deletion).
    """
    ref_len = len(ref)
    end_0 = start_0 + ref_len
    if ref:
        actual = seq[start_0:end_0].upper()
        if actual != ref.upper():
            raise ValueError(
                f"Ref-Mismatch an Position {start_0}: "
                f"erwartet {ref!r}, gefunden {actual!r}"
            )
    return seq[:start_0] + alt + seq[end_0:]


def _codon_info(
    cds_seq: str, cds_pos_0: int, alt: str | None, ref: str | None
) -> tuple[str | None, str | None, str | None, str | None, int | None]:
    """
    Berechnet ref_codon, alt_codon, ref_aa, alt_aa, prot_pos für SNVs im CDS.
    Gibt None-Tuple zurück wenn nicht anwendbar (Indels, nicht-kodierend).
    """
    if not ref or not alt or len(ref) != 1 or len(alt) != 1:
        return None, None, None, None, None
    codon_start = (cds_pos_0 // 3) * 3
    if codon_start + 3 > len(cds_seq):
        return None, None, None, None, None
    ref_codon = cds_seq[codon_start : codon_start + 3]
    mutated_cds = _apply_edit(cds_seq, cds_pos_0, ref, alt)
    alt_codon = mutated_cds[codon_start : codon_start + 3]
    ref_aa = str(Seq(ref_codon).translate())
    alt_aa = str(Seq(alt_codon).translate())
    prot_pos = codon_start // 3 + 1
    return ref_codon, alt_codon, ref_aa, alt_aa, prot_pos


class VariantSequenceMaterializer:
    """
    Converts hgvs SequenceVariant objects into actual DNA sequences.

    Requires only:
    - A SeqRepo instance (local or REST)
    - The SequenceVariant object (g. or c. type)
    - A flanking window if context is desired
    """

    def __init__(
        self,
        seqrepo_dir: str | None = "/usr/local/share/seqrepo/latest",
        seqrepo_rest_url: str | None = None,
        flanking: int = 0,
    ):
        """
        Parameters
        ----------
        seqrepo_dir
            Path to local SeqRepo snapshot. Takes precedence over REST.
        seqrepo_rest_url
            URL of seqrepo-rest-service, e.g. "http://localhost:5000".
            Used only if seqrepo_dir is None.
        flanking
            Bases of reference context to include left and right of the variant.
        """
        if seqrepo_dir:
            self._sr = SeqRepo(seqrepo_dir)
        elif seqrepo_rest_url:
            self._sr = SeqRepoRESTDataProxy(seqrepo_rest_url)
        else:
            raise ValueError("Provide either seqrepo_dir or seqrepo_rest_url")

        self.flanking = flanking

    def _apply_edit(
        self, ref_region: str, variant_pos: int, ref_len: int, edit
    ) -> str:
        """
        Apply any hgvs Edit subclass to a reference region string.

        Parameters
        ----------
        ref_region  : the fetched reference sequence (with flanking)
        variant_pos : 0-based start of the variant within ref_region
        ref_len     : length of the reference allele (end = variant_pos + ref_len)
        edit        : any hgvs.edit.Edit subclass instance
        """
        before = ref_region[:variant_pos]
        ref_seq = ref_region[variant_pos : variant_pos + ref_len]
        after = ref_region[variant_pos + ref_len :]

        match edit:

            case hgvs.edit.NARefAlt():
                # SNV, Del, Ins, Delins — standard ref→alt substitution
                alt = edit.alt or ""
                return before + alt + after

            case hgvs.edit.Dup():
                # g.100_102dup → duplicates the ref_seq in tandem
                return before + ref_seq + ref_seq + after

            case hgvs.edit.Inv():
                # g.100_102inv → reverse complement of ref_seq
                inv_seq = str(Seq(ref_seq).reverse_complement())
                return before + inv_seq + after

            case hgvs.edit.Repeat():
                # g.100_102[4] → repeat ref_seq N times
                count = edit.seq.ref_n  # the repeat count
                return before + (ref_seq * count) + after

            case _:
                raise NotImplementedError(
                    f"Edit type {type(edit).__name__} is not supported for "
                    f"sequence materialization"
                )

    def _ref_len_from_variant(self, var: SequenceVariant) -> int:
        """
        Determine the reference allele length from the variant position interval.
        For point variants (SNV, single-base del), start == end → length 1.
        For intervals (multi-base), end - start + 1.
        """
        pos = var.posedit.pos
        start = pos.start.base
        end = pos.end.base

        # Insertions: start and end are adjacent (e.g. g.5_6ins → ref_len=0)
        edit = var.posedit.edit
        if isinstance(edit, hgvs.edit.NARefAlt) and not edit.ref and edit.alt:
            return 0  # pure insertion, nothing consumed from reference

        return max(end - start + 1, 1)

    def _fetch(self, ac: str, start_0: int, end_0: int) -> str:
        """Fetch a sequence slice. Handles both local SeqRepo and REST proxy."""
        return self._sr[ac][start_0:end_0]

    def materialize(
        self, var: hgvs.sequencevariant.SequenceVariant
    ) -> tuple[str, int]:
        ac = var.ac
        pos = var.posedit.pos
        edit = var.posedit.edit

        start_1 = pos.start.base
        ref_len = self._ref_len_from_variant(var)

        # Fetch reference region (0-based for SeqRepo)
        fetch_start = max(0, start_1 - 1 - self.flanking)
        fetch_end = start_1 - 1 + max(ref_len, 1) + self.flanking
        ref_region = self._fetch(ac, fetch_start, fetch_end)

        variant_pos = (start_1 - 1) - fetch_start

        # Optional ref validation where possible
        ref_from_edit = getattr(edit, "ref", None)
        if ref_from_edit:
            actual = ref_region[variant_pos : variant_pos + ref_len]
            if actual.upper() != ref_from_edit.upper():
                print(
                    variant_pos,
                    ref_len,
                )
                raise ValueError(
                    f"Ref mismatch at {ac}:{start_1}: "
                    f"expected {ref_from_edit!r}, found {actual!r}"
                )

        mutated = self._apply_edit(ref_region, variant_pos, ref_len, edit)
        return mutated, variant_pos


class HGVSVariantGenerator:
    """
    Erzeugt Mutation-Instanzen aus HGVS-Annotationen (hgvs_c oder hgvs_g).

    Parameters
    ----------
    genomic_fasta : str | Path
        FASTA-Datei mit der genomischen Sequenz des Gens (Plus-Strand).
    cds_fasta : str | Path
        FASTA-Datei mit der CDS-Sequenz (ATG...Stop, ohne UTR).
    gene_start_0 : int
        0-based chromosomaler Start des genomischen Fragments
        (d.h. erstes Nukleotid in genomic_fasta entspricht dieser Chromosomenposition).
    transcript_ac : str
        RefSeq Transkript-Accession, z.B. "NM_000546.6".
    chrom_ac : str
        RefSeq Chromosom-Accession, z.B. "NC_000017.11".
    chromosome : str
        Chromosomname für die Mutation-Instanz, z.B. "17".
    cdot_json : str | Path
        Pfad zur cdot JSON-Datei für den Dataprovider.
    assembly : str
        Assemblyname, z.B. "GRCh38".
    strand : int
        +1 für Plus-Strand, -1 für Minus-Strand (TP53 → -1).
    cds_start_in_genomic_0 : int
        0-based Position des CDS-Starts (A von ATG) im genomic_fasta-Fragment.
        Wird gebraucht um cds_pos aus genomischer Position zu berechnen.
    flanking : int
        Anzahl Basen Flanking-Region die links/rechts an seq angehängt werden.
        0 = kein Flanking (Standard).
    seq_id_prefix : str
        Präfix für seq_id in der Mutation-Instanz.
    """

    def __init__(
        self,
        genomic: SeqIO.SeqRecord,
        cds: SeqIO.SeqRecord,
        transcript_ac: str,
        chrom_ac: str,
        cdot_json: str | Path | None = None,
        # strand: int = -1,
        cds_start_in_genomic_0: int = 0,
        flanking: int = 0,
        seq_id_prefix: str = "var",
        annotation: dict | None = None,
        genomic_flanks: tuple[int, int] | None = None,
    ):
        self.genomic = genomic
        desc = genomic.description
        splits = desc.split(":")
        strand = splits[-1]
        end = splits[-2]
        start = splits[-3]
        chromosome = splits[-4]
        assembly = splits[-5]
        self.strand = int(strand)
        self.gene_start_1 = int(start)
        self.gene_end_1 = int(end)
        if self.strand == -1:
            self.genomic.seq = self.genomic.seq.reverse_complement()
            self.strand = 1
        self.cds = cds
        self.genomic_flanks = genomic_flanks
        self.assembly = assembly
        self.transcript_ac = transcript_ac
        self.chrom_ac = chrom_ac
        self.chromosome = chromosome
        self.cds_start_in_genomic_0 = cds_start_in_genomic_0
        self.flanking = flanking
        self.seq_id_prefix = seq_id_prefix

        # hgvs setup
        self.hdp = (
            JSONDataProvider([str(cdot_json)])
            if cdot_json
            else RESTDataProvider()  # Uses API server at cdot.cc
            # else hgvs.dataproviders.uta.connect()
        )
        self.am = hgvs.assemblymapper.AssemblyMapper(
            self.hdp,
            assembly_name=assembly,
            alt_aln_method="splign",
            replace_reference=True,
        )
        self.hp = hgvs.parser.Parser()

    ############################################################################
    # Helper
    ############################################################################

    def _safe_g_to_c(self, var_g):
        try:
            return self.am.g_to_c(var_g, self.transcript_ac)
        except Exception:
            return None

    def _safe_c_to_p(self, var_c):
        if var_c is None:
            return None
        try:
            return self.am.c_to_p(var_c)
        except Exception:
            return None

    def _extract_edit(self, var) -> tuple[str, str]:
        """Extrahiert (ref, alt) aus einem SequenceVariant-Edit-Objekt."""
        print(var, type(var))
        edit = var.posedit.edit
        print("edit", edit, type(edit))

        # print(edit.ref, edit.alt)
        ref = ""
        if hasattr(edit, "ref") and edit.ref:
            ref = edit.ref if self.strand == 1 else reverse_complement(edit.ref)
        alt = ""
        if hasattr(edit, "alt") and edit.alt:
            alt = edit.alt if self.strand == 1 else reverse_complement(edit.alt)
        return ref, alt

    def _genomic_pos_to_local_0(self, g_pos_1based: int) -> int:
        """Konvertiert 1-based chromosomale Position → 0-based Index in genomic_seq."""
        print(
            "pos loc",
            g_pos_1based,
            self.strand,
            self.gene_start_1,
            self.cds_start_in_genomic_0,
            self.gene_end_1,
        )
        print("gene_Start", self.gene_start_1, "gene_end", self.gene_end_1)
        if self.strand == 1:
            return g_pos_1based - self.gene_start_1
        else:
            return self.gene_end_1 - g_pos_1based

    def _c_pos_to_cds_0(self, c_pos_1based: int) -> int:
        """Konvertiert 1-based c.-Position → 0-based Index in cds_seq."""
        return c_pos_1based - 1

    def _region_type(self, var_c) -> str:
        """Bestimmt ob die Mutation exonisch oder intronisch ist."""
        if var_c is None:
            return "intergenic"
        pos = var_c.posedit.pos.start
        # Intronische Positionen haben einen Offset != 0
        if hasattr(pos, "offset") and pos.offset != 0:
            return "intron"
        # UTR
        base = pos.base
        if base <= 0:
            return "5UTR"
        # Prüfe ob nach Stop (c.*N)
        if hasattr(pos, "datum") and str(pos.datum) == "SEQ_STOP":
            return "3UTR"
        return "exon"

    def mutalyzer_to_hgvs(
        self, mutalyzer_str: str, translate_to: Literal["c", "g"] = "c"
    ) -> tuple[str, str | None, str | None]:
        """
        Converts Mutalyzer-style strings into valid HGVS strings, hopefully.

        Examples:
            NC_000017.11(NM_000546.6):c.920-3_923del
            → NM_000546.6:c.920-3_923del

        Parameters
        ----------
        s : str
        prefer : "c" or "g"
            whether to prefer transcript or genomic reference

        Returns
        -------
        str
        """

        mutalyzer_str = mutalyzer_str.strip()
        # remove protein annotation
        # mutalyzer_str = re.sub(r"\s*\(p\.[^)]+\)", "", mutalyzer_str)

        # -----------------------------
        # Case 1: NC_...(NM_...):c....
        # -----------------------------
        match = re.match(r"(NC_[^()]+)\((NM_[^)]+)\):(c\..+)", mutalyzer_str)
        if match:
            genomic, transcript, change = match.groups()
            hgvs = (
                f"{transcript}:{change}"
                if translate_to == "c"
                else f"{genomic}:{change}"
            )
            return hgvs, transcript, genomic
        # -----------------------------
        # Case 2: NM_...(GENE):c....
        # -----------------------------
        match = re.match(r"(NM_[^(]+)\([^)]*\):(c\..+)", mutalyzer_str)
        if match:
            transcript, change = match.groups()
            return f"{transcript}:{change}", transcript, None

        # -----------------------------
        # Case 3: already valid
        # -----------------------------
        if re.match(r"[A-Z]{2}_[0-9]+\.[0-9]+:(c\.|g\.).+", mutalyzer_str):
            return mutalyzer_str, None, None

        raise ValueError(f"Unrecognized format: {mutalyzer_str}")

    ############################################################################
    # HGVS translation
    ############################################################################

    def from_hgvs_c(self, hgvs_or_mutalizer_c_str: str) -> Mutation:
        """Erzeugt eine Mutation-Instanz aus einem c.-HGVS-String."""
        hgvs_c_str, transcript, _ = self.mutalyzer_to_hgvs(
            hgvs_or_mutalizer_c_str, translate_to="c"
        )
        var_c = self.hp.parse_hgvs_variant(hgvs_c_str)
        var_g = self.am.c_to_g(var_c)  # , alt_ac=transcript)
        var_p = self._safe_c_to_p(var_c)
        return self._build_mutation(
            var_c=var_c,
            var_g=var_g,
            var_p=var_p,
            hgvs_c_str=hgvs_c_str,
        )

    def from_hgvs_g(self, hgvs_or_mutalizer_g_str: str) -> Mutation:
        """Erzeugt eine Mutation-Instanz aus einem g.-HGVS-String."""
        hgvs_g_str, _, _ = self.mutalyzer_to_hgvs(
            hgvs_or_mutalizer_g_str, translate_to="g"
        )
        var_g = self.hp.parse_hgvs_variant(hgvs_g_str)
        var_c = self._safe_g_to_c(var_g)
        var_p = self._safe_c_to_p(var_c) if var_c else None
        return self._build_mutation(
            var_c=var_c,
            var_g=var_g,
            var_p=var_p,
            hgvs_g_str=hgvs_g_str,
        )

    ############################################################################
    # Mutation Factory
    ############################################################################

    # def _build_seq_with_flanking(
    #     self,
    #     local_start_0: int,
    #     ref: str,
    #     alt: str,
    # ) -> tuple[str, int]:
    #     """
    #     Baut die mutierte Sequenz mit optionalen Flanking-Regionen.

    #     Returns
    #     -------
    #     seq : str
    #         Mutiertes Sequenzfragment (mit Flanking wenn self.flanking > 0).
    #     mutation_pos : int
    #         0-based Position der Mutation im zurückgegebenen seq-String.
    #     """
    #     if self.genomic_flanks is not None:
    #         flank_left, flank_right = self.genomic_flanks
    #         frag_start = self._genomic_pos_to_local_0(flank_left)
    #         frag_end = self._genomic_pos_to_local_0(flank_right)
    #         frag_start, frag_end = min(frag_start, frag_end), max(
    #             frag_start, frag_end
    #         )
    #         fragment = self.genomic.seq[frag_start:frag_end]
    #         local_in_frag = local_start_0 - frag_start
    #         print(
    #             "from local",
    #             self.genomic.seq[local_start_0 : local_start_0 + 20],
    #         )
    #         print(frag_start, frag_end, self.genomic_flanks, self.gene_start_1)
    #         print(ref, fragment, alt, frag_start, local_start_0, frag_end)
    #     else:
    #         flank = self.flanking
    #         frag_start = max(0, local_start_0 - flank)
    #         frag_end = min(
    #             len(self.genomic),
    #             local_start_0 + max(len(ref), 1) + flank,
    #         )
    #         print(
    #             "from local",
    #             self.genomic.seq[local_start_0 : local_start_0 + 20],
    #         )
    #         print(frag_start, frag_end)
    #         fragment = self.genomic.seq[frag_start:frag_end]
    #         print(ref, fragment, alt, frag_start, local_start_0, frag_end)
    #         # apply mutation in fragment
    #         local_in_frag = local_start_0 - frag_start
    #     mutated_fragment = _apply_edit(fragment, local_in_frag, ref, alt)
    #     print(",mutated_fragment", mutated_fragment)
    #     mutation_pos = local_in_frag
    #     return mutated_fragment, mutation_pos

    def _get_local_reference(
        self,
        genomic_pos_1: int,
        # local_start_0: int,
        # ref: str,
    ) -> tuple[str, int]:
        """
        Baut die mutierte Sequenz mit optionalen Flanking-Regionen.

        Returns
        -------
        seq : str
            Mutiertes Sequenzfragment (mit Flanking wenn self.flanking > 0).
        mutation_pos : int
            0-based Position der Mutation im zurückgegebenen seq-String.
        """
        local_pos_0 = self._genomic_pos_to_local_0(genomic_pos_1)
        if self.genomic_flanks is not None:
            flank_left, flank_right = self.genomic_flanks
            frag_start = self._genomic_pos_to_local_0(flank_left)
            frag_end = self._genomic_pos_to_local_0(flank_right)
            local_pos_0 = self._genomic_pos_to_local_0(genomic_pos_1)
            frag_start, frag_end = min(frag_start, frag_end), max(
                frag_start, frag_end
            )
            print(
                flank_left,
                flank_right,
                frag_start,
                frag_end,
                local_pos_0,
                genomic_pos_1,
            )
            fragment = self.genomic.seq[frag_start:frag_end]
            print("Fragment", fragment)
            local_in_frag = local_pos_0 - frag_start
        else:

            fragment = self.genomic.seq
            local_in_frag = local_pos_0
        return fragment, local_in_frag

    def _apply_edit(
        self, ref_region: str, variant_pos: int, ref_len: int, edit
    ) -> str:
        """
        Apply any hgvs Edit subclass to a reference region string.

        Parameters
        ----------
        ref_region  : the fetched reference sequence (with flanking)
        variant_pos : 0-based start of the variant within ref_region
        ref_len     : length of the reference allele (end = variant_pos + ref_len)
        edit        : any hgvs.edit.Edit subclass instance
        """
        before = ref_region[:variant_pos]
        ref_seq = ref_region[variant_pos : variant_pos + ref_len]
        after = ref_region[variant_pos + ref_len :]

        match edit:

            case hgvs.edit.NARefAlt():
                # SNV, Del, Ins, Delins — standard ref→alt substitution
                alt = edit.alt or ""
                return before + alt + after

            case hgvs.edit.Dup():
                # g.100_102dup → duplicates the ref_seq in tandem
                return before + ref_seq + ref_seq + after

            case hgvs.edit.Inv():
                # g.100_102inv → reverse complement of ref_seq
                inv_seq = str(Seq(ref_seq).reverse_complement())
                return before + inv_seq + after

            case hgvs.edit.Repeat():
                # g.100_102[4] → repeat ref_seq N times
                count = edit.seq.ref_n  # the repeat count
                return before + (ref_seq * count) + after

            case _:
                raise NotImplementedError(
                    f"Edit type {type(edit).__name__} is not supported for "
                    f"sequence materialization"
                )

    def _ref_len_from_variant(self, var: SequenceVariant) -> int:
        """
        Determine the reference allele length from the variant position interval.
        For point variants (SNV, single-base del), start == end → length 1.
        For intervals (multi-base), end - start + 1.
        """
        pos = var.posedit.pos
        start = pos.start.base
        end = pos.end.base

        # Insertions: start and end are adjacent (e.g. g.5_6ins → ref_len=0)
        edit = var.posedit.edit
        if isinstance(edit, hgvs.edit.NARefAlt) and not edit.ref and edit.alt:
            return 0  # pure insertion, nothing consumed from reference

        return max(end - start + 1, 1)

    def materialize(
        self, var: hgvs.sequencevariant.SequenceVariant
    ) -> tuple[str, int]:
        ac = var.ac
        pos = var.posedit.pos
        edit = var.posedit.edit

        start_1 = pos.start.base
        g_start_1 = (
            pos.start.base if self.strand == 1 else pos.end.base
        )  # 1-based chromosomal

        ref_len = self._ref_len_from_variant(var)

        # Fetch reference region (0-based for SeqRepo)
        # fetch_start = max(0, start_1 - 1 - self.flanking)
        # fetch_end = start_1 - 1 + max(ref_len, 1) + self.flanking
        # ref_region = self.genomic.seq  # self._fetch(ac, fetch_start, fetch_end)
        ref_region, local_variant_pos_0 = self._get_local_reference(g_start_1)
        # variant_pos = (start_1 - 1) - fetch_start

        # Optional ref validation where possible
        print(edit)
        ref_from_edit = getattr(edit, "ref", None)
        if ref_from_edit:
            # actual = ref_region[variant_pos : variant_pos + ref_len]
            actual = ref_region[
                local_variant_pos_0 : local_variant_pos_0 + ref_len
            ]
            if actual.upper() != ref_from_edit.upper():
                print(ref_len, ref_region, local_variant_pos_0, g_start_1)
                raise ValueError(
                    f"Ref mismatch at {ac}:{g_start_1}: "
                    f"expected {ref_from_edit!r}, found {actual!r}"
                )

        mutated = self._apply_edit(
            ref_region, local_variant_pos_0, ref_len, edit
        )
        return mutated, local_variant_pos_0
        # # Sequenzen laden
        # self.genomic_seq: str = str(
        #     next(SeqIO.parse(genomic_fasta, "fasta")).seq
        # ).upper()
        # self.cds_seq: str = str(
        #     next(SeqIO.parse(cds_fasta, "fasta")).seq
        # ).upper()

    def _build_mutation(
        self,
        var_g,
        var_c=None,
        var_p=None,
        hgvs_c_str: str | None = None,
        hgvs_g_str: str | None = None,
    ) -> Mutation:

        # HGVS Strings
        hgvs_g = hgvs_g_str or str(var_g)
        hgvs_c = hgvs_c_str or (str(var_c) if var_c else None)
        hgvs_p = str(var_p) if var_p else None
        print("var_g", var_g)
        print("var_c", var_c)
        print("var_p", var_p)
        print("hgvs_g", hgvs_g)

        # extract edit from g
        g_ref, g_alt = self._extract_edit(var_g)
        g_start_1 = (
            var_g.posedit.pos.start.base
            if self.strand == 1
            else var_g.posedit.pos.end.base
        )  # 1-based chromosomal
        local_start_0 = self._genomic_pos_to_local_0(g_start_1)
        print(
            "Where it wants to place it",
            self.genomic.seq[local_start_0 : local_start_0 + 20],
            self.genomic.seq[7673605 - 20 : local_start_0],
        )
        print(
            "this is the local start",
            local_start_0,
            g_start_1,
            self.genomic.seq[3435 - len(g_ref) : 3435],
        )
        print("g_ref", g_ref, "g_alt", g_alt)
        # extract edit from c (if available) and determine cds_pos
        print("strand", self.strand)
        cds_pos_0: int | None = None
        c_ref, c_alt = (
            (g_ref, g_alt)
            if self.strand == 1
            else (
                reverse_complement(g_ref),
                reverse_complement(g_alt),
            )
        )  # Fallback - this is dump
        print("c_ref", c_ref, "c_alt", c_alt)
        # print(
        #     var_c.posedit,
        #     dir(var_c.posedit),
        #     var_c.posedit.pos,
        #     dir(var_c.posedit.pos),
        #     var_c.posedit.edit,
        #     dir(var_c.posedit.edit),
        #     "ref",
        #     var_c.posedit.edit.ref,
        #     type(var_c.posedit.edit.ref),
        #     "alt",
        #     var_c.posedit.edit.alt,
        #     type(var_c.posedit.edit.alt),
        #     "base",
        #     var_c.posedit.pos.start,
        #     type(var_c.posedit.pos.start),
        # )
        if var_c is not None:
            c_ref, c_alt = self._extract_edit(var_c)

            c_start = var_c.posedit.pos.start
            # Nur für exonische Positionen ohne Intron-Offset
            if not (hasattr(c_start, "offset") and c_start.offset != 0):
                base = c_start.base
                if 1 <= base <= len(self.cds):
                    cds_pos_0 = self._c_pos_to_cds_0(base)

        print(
            "aaaaa",
            g_ref,
            c_ref,
            "alt",
            g_alt,
            c_alt,
            reverse_complement(g_ref),
            reverse_complement(c_ref),
        )
        print("c_ref", c_ref, "c_alt", c_alt)

        # Mutation type
        mutation_type = _mut_type(g_ref, g_alt)

        # Codon-Info (only for exonic SNVs)
        ref_codon = alt_codon = ref_aa = alt_aa = None
        prot_pos: int | None = None

        if cds_pos_0 is not None and mutation_type == "SNV":
            ref_codon, alt_codon, ref_aa, alt_aa, prot_pos = _codon_info(
                self.cds, cds_pos_0, c_alt, c_ref
            )

        # ── Region-Typ ────────────────────────────────────────────────────
        region_type = self._region_type(var_c)

        # ── Sequenz mit Flanking ──────────────────────────────────────────
        seq, mutation_pos = self.materialize(var_g)
        # seq, mutation_pos = self._build_seq_with_flanking(
        #     local_start_0, g_ref, g_alt
        # )
        print(seq)
        # ── seq_id ────────────────────────────────────────────────────────
        seq_id = f"{self.seq_id_prefix}_{hgvs_c or hgvs_g}"
        return Mutation(
            seq_id=seq_id,
            genomic_id=self.chrom_ac,
            chromosome=self.chromosome,
            region_type=region_type,
            mutation_pos=mutation_pos,
            genomic_pos=g_start_1,  # 1-based
            cds_pos=cds_pos_0,  # 0-based
            ref=g_ref,
            alt=g_alt,
            mutation_type=mutation_type,
            seq=seq,
            ref_codon=ref_codon,
            alt_codon=alt_codon,
            ref_aa=ref_aa,
            alt_aa=alt_aa,
            prot_pos=prot_pos,
            hgvs_g=hgvs_g,
            hgvs_c=hgvs_c,
            hgvs_p=hgvs_p,
        )


############################################################################
# ClinVar
############################################################################


class ClinVarMutator:

    def __init__(self, annotator, cds_to_genomic_map):
        """
        cds_to_genomic_map: dict[int, int]
            maps CDS index (0-based) → genomic index (0-based)
        """
        self.annotator = annotator
        self.map = cds_to_genomic_map

    # --------------------------------------------------
    # 🧬 Parse HGVS c. string
    # --------------------------------------------------
    def parse_hgvs_c(self, hgvs: str):
        """
        Returns structured dict
        """
        hgvs = hgvs.split(":")[-1]  # remove transcript prefix

        # remove protein part
        hgvs = hgvs.split(" ")[0]

        # patterns
        patterns = {
            "del": r"c\.(.+)del$",
            "dup": r"c\.(.+)dup$",
            "delins": r"c\.(.+)delins([ACGT]+)",
            "sub": r"c\.(\d+)([ACGT])>([ACGT])",
        }

        for k, p in patterns.items():
            m = re.match(p, hgvs)
            if m:
                return k, m.groups()

        raise ValueError(f"Unsupported HGVS: {hgvs}")

    # --------------------------------------------------
    # 🧬 Convert c. position → genomic
    # --------------------------------------------------
    def cds_to_genomic(self, cds_pos):
        """
        cds_pos can be:
        - "920"
        - "993+4"
        - "920-2"
        """
        m = re.match(r"(\d+)([+-]\d+)?", cds_pos)

        base = int(m.group(1)) - 1
        offset = int(m.group(2)) if m.group(2) else 0

        genomic = self.map[base] + offset
        return genomic

    # --------------------------------------------------
    # 🧬 Handle ranges
    # --------------------------------------------------
    def parse_range(self, pos_str):
        if "_" in pos_str:
            start, end = pos_str.split("_")
            return self.cds_to_genomic(start), self.cds_to_genomic(end)
        else:
            g = self.cds_to_genomic(pos_str)
            return g, g

    # --------------------------------------------------
    # 🧬 Apply variant
    # --------------------------------------------------
    def apply_variant(self, seq, start, end, op, alt=None):

        if op == "del":
            return seq[:start] + seq[end + 1 :]

        elif op == "dup":
            dup = seq[start : end + 1]
            return seq[: end + 1] + dup + seq[end + 1 :]

        elif op == "delins":
            return seq[:start] + alt + seq[end + 1 :]

        elif op == "sub":
            return seq[:start] + alt + seq[start + 1 :]

        else:
            raise ValueError(op)

    # --------------------------------------------------
    # 🚀 MAIN
    # --------------------------------------------------
    def generate_from_hgvs(
        self,
        hgvs_string,
        record,
        seq_id,
        genomic_id,
    ):

        seq = str(record.seq)

        op, groups = self.parse_hgvs_c(hgvs_string)

        if op == "sub":
            pos, ref, alt = groups
            start = self.cds_to_genomic(pos)
            end = start

        elif op == "del":
            pos = groups[0]
            start, end = self.parse_range(pos)
            alt = ""

        elif op == "dup":
            pos = groups[0]
            start, end = self.parse_range(pos)
            alt = None

        elif op == "delins":
            pos, alt = groups
            start, end = self.parse_range(pos)

        else:
            raise ValueError(op)

        mutated_seq = self.apply_variant(seq, start, end, op, alt)

        ref = seq[start : end + 1] if op != "sub" else seq[start]

        return self.annotator.annotate(
            seq_id=seq_id,
            genomic_id=genomic_id,
            mutation_pos=start,
            genomic_pos=start,
            region_type="exon",  # TODO: map properly
            ref=ref,
            alt=alt if alt else "",
            mut_type=op.upper(),
            seq=mutated_seq,
        )
