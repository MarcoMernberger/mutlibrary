from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from functools import cached_property
from pathlib import Path

import hgvs.assemblymapper  # type: ignore[import]
import hgvs.edit  # type: ignore[import]
import hgvs.location  # type: ignore[import]
import hgvs.parser  # type: ignore[import]
import hgvs.posedit  # type: ignore[import]
import hgvs.sequencevariant  # type: ignore[import]
import requests  # type: ignore[import]
from Bio import SeqIO  # type: ignore[import]
from Bio.Seq import Seq  # type: ignore[import]
from hgvs.sequencevariant import SequenceVariant  # type: ignore[import]
from mmalignments.services.io import from_json  # type: ignore[import]
from mmalignments.utils.utils import reverse_complement  # type: ignore[import]

from .annotator import (
    HGVSMutationAnnotator,
    # Mutalyzer,
    Mutation,
    get_info_from_description,
)


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


# def fix_hgvs_c_annotation(hgvs_c: str) -> str:
#     """A helper to correct the names in the table. remove this."""
#     return hgvs_c.replace("(TP53)", "").split()[0]


class VariantG:

    def __init__(self, var: SequenceVariant):
        if var is None:
            raise ValueError("Nontype variant provided")
        self.var = var
        if var.type != "g":
            raise ValueError("Variant is not genomic (g.)")
        self.pe = var.posedit
        self.pos = self.pe.pos
        self.edit = self.pe.edit
        self.start = self.pos.start.base
        self.end = self.pos.end.base if self.pos.end else self.start
        self.edit_type = (
            self.edit.type if hasattr(self.edit, "type") else "has none"
        )
        self.ref = getattr(self.edit, "ref", "")
        self.alt = getattr(self.edit, "alt", "")

    @cached_property
    def genomic_start_1(self) -> int:
        g_start_1 = self.start
        if self.edit_type == "ins":
            g_start_1 = self.end if self.end else g_start_1 + 1
        return g_start_1

    @cached_property
    def ref_len(self):
        edit = self.var.posedit.edit
        pos = self.var.posedit.pos

        if self.edit_type == "delins":
            return self.end - self.start + 1
        # 1. Deletion (del)
        if isinstance(edit, hgvs.edit.NARefAlt) and edit.ref and not edit.alt:
            return len(edit.ref)

        # 2. Insertion (ins) - WICHTIG: ref_len = 0!
        if isinstance(edit, hgvs.edit.NARefAlt) and not edit.ref and edit.alt:
            return 0

        # 3. Substitution / delins
        if isinstance(edit, hgvs.edit.NARefAlt) and edit.ref:
            return len(edit.ref)

        # 4. Duplication
        if isinstance(edit, hgvs.edit.Dup):
            start = pos.start.base
            end = pos.end.base if pos.end and pos.end.base else start
            return end - start + 1

        # 5. Fallback
        start = pos.start.base
        end = pos.end.base if pos.end and pos.end.base else start
        return end - start + 1

    def validate_ref(self, ref_region, local_variant_pos_0):
        if self.ref:
            # actual = ref_region[variant_pos : variant_pos + ref_len]
            actual = ref_region[
                local_variant_pos_0 : local_variant_pos_0 + self.ref_len
            ]
            if actual.upper() != self.ref.upper():
                print(
                    self.ref_len,
                    ref_region,
                    local_variant_pos_0,
                    self.genomic_start_1,
                )
                raise ValueError(
                    f"Ref mismatch at {self.var.ac}:{self.genomic_start_1} for {self.var} "  # noqa: E501
                    f"expected {self.ref!r}, found {actual!r}"
                )

    def print_hgvs_details(self):
        print("=== HGVS Variant Details ===")

        print(f"Accession (ac): {self.var.ac}")
        print(f"Type: {self.var.type}")

        pe = self.var.posedit
        if pe is None:
            print("No posedit information")
            return

        # --- Position ---
        pos = pe.pos
        print("\n--- Position ---")
        print(f"Raw pos object: {pos}")

        if hasattr(pos, "start") and hasattr(pos, "end"):
            print(f"Start: {pos.start}")
            print(f"End: {pos.end}")

            if hasattr(pos.start, "base"):
                print(f"Start base: {pos.start.base}")
            if hasattr(pos.end, "base"):
                print(f"End base: {pos.end.base}")

        # --- Edit ---
        edit = pe.edit
        print("\n--- Edit ---")
        print(f"Edit type: {type(edit).__name__}")
        print(f"Raw edit: {edit}")
        print(f"Ref_len: {self.ref_len}")

        if hasattr(edit, "ref"):
            print(f"Reference: {edit.ref}")
        if hasattr(edit, "alt"):
            print(f"Alternate: {edit.alt}")

        if hasattr(edit, "type"):
            print(f"Edit type string: {edit.type}")

        print("\n=== END ===")


class HGVSVariantGenerator:
    """
    A class to generate Mutation instances from HGVS annotations. It can
    handle both hgvs_c and hgvs_g as input, and will annotate the resulting
    Mutation with hgvs_g, hgvs_c, hgvs_p, and hgvs_r if possible.

    Parameters
    ----------
    genomic : SeqIO.SeqRecord
        Genomic sequence record containing the gene of interest. The
        description should contain information about the strand, start, end,
        chromosome, and assembly (e.g., "chr17:7661779-7687550:-1").
    cds : SeqIO.SeqRecord
        CDS sequence record containing the coding sequence of the gene.
    cdot_json : str | Path | None, optional
        Path to the cdot JSON file for the data provider, by default None.
    flanking : int, optional
        Number of bases in the flanking region to be added to the left/right
        of the sequence, by default 0 (no flanking).
    genomic_flanks : tuple[int, int] | None, optional
        Tuple specifying the number of bases to add as flanks to the genomic
        sequence (left_flank, right_flank), by default None. If provided,
        this will override the `flanking` parameter for the genomic
        sequence.
    annotation_json : Path | str | None, optional
        Path to the annotation JSON file, by default None. This file should
        contain a mapping from CDS IDs to transcript and chromosome
        accessions, as well as protein accessions. If not provided, these
        fields will be left empty in the Mutation instances.
    """

    def __init__(
        self,
        genomic: SeqIO.SeqRecord,
        cds: SeqIO.SeqRecord,
        cdot_json: str | Path | None = None,
        flanking: int = 0,
        genomic_flanks: tuple[int, int] | None = None,
        annotation_json: Path | str | None = None,
    ):
        self.genomic = genomic
        info = get_info_from_description(genomic)
        self.strand = int(info["strand"])
        self.gene_start_1 = int(info["start"])
        self.gene_end_1 = int(info["end"])
        self.chromosome = info["chromosome"]
        self.assembly = info["assembly"]
        if self.strand == -1:
            self.genomic.seq = self.genomic.seq.reverse_complement()
        self.cds = cds
        self.genomic_flanks = genomic_flanks
        # self.transcript_ac = transcript_ac
        # self.chrom_ac = chrom_ac
        self.flanking = flanking
        self.cds_start_in_genomic_0 = self._find_cds_start()
        # hgvs setup
        # self.hdp = (
        #     JSONDataProvider([str(cdot_json)])
        #     if cdot_json
        #     else RESTDataProvider()  # Uses API server at cdot.cc
        #     # else hgvs.dataproviders.uta.connect()
        # )
        # self.am = hgvs.assemblymapper.AssemblyMapper(
        #     self.hdp,
        #     assembly_name=self.assembly,
        #     alt_aln_method="splign",
        #     replace_reference=True,
        # )
        # self.hp = hgvs.parser.Parser()
        # self.hn = hgvs.normalizer.Normalizer(self.hdp)
        annotation = from_json(Path(annotation_json)) if annotation_json else {}
        self.transcript_id = annotation[self.cds.id].get("transcript_id", "")
        self.transcript_ac = annotation[self.cds.id].get("transcript_ac", "")
        self.chrom_ac = annotation[self.cds.id].get("chromosome_ac", "")
        self.protein_ac = annotation[self.cds.id].get("protein_ac", "")
        self.annotator = HGVSMutationAnnotator(
            assembly=self.assembly, cdot_json=cdot_json
        )

    ############################################################################
    # Helper
    ############################################################################

    def _find_cds_start(self) -> int:  # TODO infer from coordinates
        cds_start_in_genomic_0 = self.genomic.seq.find(self.cds.seq[:15])
        return cds_start_in_genomic_0

    # def _safe_g_to_c(self, var_g):
    #     try:
    #         return self.am.g_to_c(var_g, self.transcript_ac)
    #     except Exception:
    #         return None

    # def _safe_c_to_p(self, var_c):
    #     if var_c is None:
    #         return None
    #     try:
    #         return self.am.c_to_p(var_c)
    #     except Exception:
    #         return None

    def _extract_edit(self, var) -> tuple[str, str]:
        """Extrahiert (ref, alt) aus einem SequenceVariant-Edit-Objekt."""
        # print(var, type(var))
        edit = var.posedit.edit
        # print("edit", edit, type(edit))

        # print(edit.ref, edit.alt)
        ref = ""
        if hasattr(edit, "ref") and edit.ref:
            ref = edit.ref if self.strand == 1 else reverse_complement(edit.ref)
        alt = ""
        if hasattr(edit, "alt") and edit.alt:
            alt = edit.alt if self.strand == 1 else reverse_complement(edit.alt)
        return ref, alt

    def _genomic_pos_to_local_0(self, g_pos_1based: int) -> int:
        """convert 1-based chromosomal position → 0-based index in
        genomic_seq."""
        return g_pos_1based - self.gene_start_1

    def _c_pos_to_cds_0(self, c_pos_1based: int) -> int:
        """convert 1-based c.-position → 0-based index in cds_seq."""
        return c_pos_1based - 1

    def _region_type(self, var_c) -> str:
        """determine if the mutation is exonic or intronic."""
        if var_c is None:
            return "intergenic"
        pos = var_c.posedit.pos.start
        # Intronic positions have an offset != 0
        if hasattr(pos, "offset") and pos.offset != 0:
            return "intron"
        # UTR
        base = pos.base
        if base <= 0:
            return "5UTR"
        # Check if after stop (c.*N)
        if hasattr(pos, "datum") and str(pos.datum) == "SEQ_STOP":
            return "3UTR"
        return "exon"

    # def mutalyzer_to_hgvs(
    #     self, mutalyzer: str, translate_to: Literal["c", "g"] = "c"
    # ) -> tuple[str, str | None, str | None]:
    #     """
    #     Converts Mutalyzer-style strings into valid HGVS strings, hopefully.

    #     Examples:
    #         NC_000017.11(NM_000546.6):c.920-3_923del
    #         → NM_000546.6:c.920-3_923del

    #     Parameters
    #     ----------
    #     s : str
    #     prefer : "c" or "g"
    #         whether to prefer transcript or genomic reference

    #     Returns
    #     -------
    #     str
    #     """

    #     mutalyzer_str = mutalyzer_str.strip()
    #     # remove protein annotation
    #     # mutalyzer_str = re.sub(r"\s*\(p\.[^)]+\)", "", mutalyzer_str)

    #     # -----------------------------
    #     # Case 1: NC_...(NM_...):c....
    #     # -----------------------------
    #     match = re.match(r"(NC_[^()]+)\((NM_[^)]+)\):(c\..+)", mutalyzer_str)
    #     if match:
    #         genomic, transcript, change = match.groups()
    #         hgvs = (
    #             f"{transcript}:{change}"
    #             if translate_to == "c"
    #             else f"{genomic}:{change}"
    #         )
    #         return hgvs, transcript, genomic
    #     # -----------------------------
    #     # Case 2: NM_...(GENE):c....
    #     # -----------------------------
    #     match = re.match(r"(NM_[^(]+)\([^)]*\):(c\..+)", mutalyzer_str)
    #     if match:
    #         transcript, change = match.groups()
    #         return f"{transcript}:{change}", transcript, None

    #     # -----------------------------
    #     # Case 3: already valid
    #     # -----------------------------
    #     if re.match(r"[A-Z]{2}_[0-9]+\.[0-9]+:(c\.|g\.).+", mutalyzer_str):
    #         return mutalyzer_str, None, None

    #     raise ValueError(f"Unrecognized format: {mutalyzer_str}")

    ############################################################################
    # Fetch ClinVar
    ############################################################################

    def vcv_to_variation_id(self, vcv_id: str) -> str:
        return vcv_id.replace("VCV", "").lstrip("0")

    def vcv_to_id(self, variation_id: str) -> str:
        """VCV000428872 -> 428872"""

        return re.sub(r"^VCV0*", "", variation_id)

    # def get_spdi_variation_api(self, variation_id):
    #     url = f"https://api.ncbi.nlm.nih.gov/variation/v0/vcv/{variation_id}"
    #     params = {"assembly": "GRCh38"}
    #     r = requests.get(url, params=params)
    #     print(r.status_code)
    #     print(r.text[:3000])
    #     data = r.json()

    #     # SPDI is directly available
    #     spdi_list = data.get("data", {}).get("spdi", [])
    #     for spdi in spdi_list:
    #         print(spdi)
    #     return spdi_list

    def get_spdi_from_clinvar(
        self, variation_id: str, assembly: str = "GRCh38"
    ) -> list[str]:
        vid = self.vcv_to_id(variation_id)
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/?redir=vcv&variant={variation_id}"
        params = {
            # "db": "clinvar",
            # "id": vid,
            # "rettype": "clinvarset",  # "vcv",
        }
        r = requests.get(url, params=params)
        r.raise_for_status()
        print("url", r.url)
        print(r.text)
        root = ET.fromstring(r.text)
        spdis = []
        for loc in root.iter("SequenceLocation"):
            if loc.get("Assembly") != assembly:
                continue
            accession = loc.get("Accession")
            pos_vcf = loc.get("positionVCF")  # 1-based
            ref = loc.get("referenceAlleleVCF")
            alt = loc.get("alternateAlleleVCF")

            if accession and pos_vcf and ref and alt:
                # SPDI is 0-based, VCF is 1-based
                spdi = f"{accession}:{int(pos_vcf) - 1}:{ref}:{alt}"
                spdis.append(spdi)

        return spdis

    ############################################################################
    # HGVS translation
    ############################################################################

    def from_vcv_hgvs(
        self,
        vcv_accession: str,
        hgvs: str,
        genomic_range: tuple[int, int] | None = None,
    ) -> Mutation:
        """Erzeugt eine Mutation-Instanz aus einem c.-HGVS-String."""
        spdi = self.get_spdi_from_clinvar(vcv_accession)
        if not spdi:
            raise ValueError(f"No SPDI found for {vcv_accession}")
        return self.from_hgvs(
            spdi[0], genomic_range=genomic_range
        )  # Assuming the SPDI can be directly parsed as HGVS

    def from_vcv(
        self,
        vcv_accession: str,
        genomic_range: tuple[int, int] | None = None,
        hgvs: str | None = None,
    ) -> Mutation:
        """Erzeugt eine Mutation-Instanz aus einem c.-HGVS-String."""
        spdi = self.get_spdi_from_clinvar(vcv_accession)
        if not spdi:
            raise ValueError(f"No SPDI found for {vcv_accession}")
        return self.from_spdi(
            spdi[0],
            genomic_range=genomic_range,
            var_id=vcv_accession,
            hgvs=hgvs,
        )

    def from_spdi(
        self,
        spdi_str: str,
        genomic_range: tuple[int, int] | None = None,
        var_id: str | None = None,
        hgvs: str | None = None,
    ) -> Mutation:
        """Erzeugt eine Mutation-Instanz aus einem c.-HGVS-String."""
        var_g = self.annotator.spdi_to_var_g(spdi_str)
        if var_g is None:
            raise ValueError(f"Could not parse genomic HGVS for {spdi_str}")
        print("var_g", var_g, dir(var_g))
        # print("var_g", var_g, dir(var_g))
        # if (
        #     hasattr(var_g, "posedit")
        #     and var_g.posedit
        #     and hasattr(var_g.posedit, "pos")
        # ):
        #     print("post", var_g.posedit.pos, dir(var_g.posedit.pos))
        var_c = (
            self.annotator.g_to_c(var_g, self.transcript_ac)
            or self.annotator.hp.parse_hgvs_variant(hgvs)
            if hgvs
            else None
        )
        var_p = self.annotator.c_to_p(var_c, self.protein_ac)
        # try:
        #     # print("var_c", var_c, dir(var_c))
        # except hgvs.exceptions.HGVSParseError as e:
        #     print(f"Error parsing HGVS c variant: {e}")
        #     raise
        # if var_c:
        # try:
        #     var_r = self.annotator.c_to_r(var_c)  # , self.protein_ac)
        # except Exception as e:
        #     print(f"Error parsing HGVS r variant: {e}")
        #     var_r = None
        # try:
        # except Exception as e:
        #     print(f"Error parsing HGVS p variant: {e}")
        #     var_p = None
        return self._build_mutation(
            var_g=var_g,
            var_c=var_c,
            var_p=var_p,
            genomic_range=genomic_range,
            var_id=var_id,
            # normalized=None,
            # var_r=var_r,
        )

    # def from_hgvs_g(self, hgvs_or_mutalyzer_g_str: str) -> Mutation:
    #     """Erzeugt eine Mutation-Instanz aus einem g.-HGVS-String."""
    #     hgvs_g_str, _, _ = self.mutalyzer_to_hgvs(
    #         hgvs_or_mutalyzer_g_str, translate_to="g"
    #     )
    #     var_g = self.hp.parse_hgvs_variant(hgvs_g_str)
    #     var_c = self._safe_g_to_c(var_g)
    #     var_p = self._safe_c_to_p(var_c) if var_c else None
    #     return self._build_mutation(
    #         var_c=var_c,
    #         var_g=var_g,
    #         var_p=var_p,
    #         hgvs_g_str=hgvs_g_str,
    #     )

    def from_hgvs(
        self,
        hgvs_or_mutalyzer_c_str: str,
        genomic_range: tuple[int, int] | None = None,
        var_id: str | None = None,
    ) -> Mutation:
        """Erzeugt eine Mutation-Instanz aus einem c.-HGVS-String."""
        try:
            var_c = self.annotator.hp.parse_hgvs_variant(
                hgvs_or_mutalyzer_c_str
            )
        except hgvs.exceptions.HGVSParseError as e:
            print(
                f"Error parsing HGVS g variant: {e}, input: {hgvs_or_mutalyzer_c_str}"
            )
            raise
        var_g = self.annotator.c_to_g(var_c)
        if var_g is None:
            raise ValueError(
                f"Could not parse genomic HGVS for {hgvs_or_mutalyzer_c_str}"
            )
        print(hgvs_or_mutalyzer_c_str, var_g, var_c)
        var_p = self.annotator.c_to_p(
            var_c, self.protein_ac
        )  # , self.protein_ac)
        # hgvs_c_str, transcript, _ = self.mutalyzer_to_hgvs(
        #     hgvs_or_mutalyzer_c_str, translate_to="c"
        # )
        # normalized = Mutalyzer.normalize(hgvs_or_mutalyzer_c_str)
        # print("mutalyzer", normalized)
        # hgvs_compliant = Mutalyzer.mutalyzer_to_hgvs(normalized)
        # print("compliant", hgvs_compliant)
        # var_c = None
        # var_g = None
        # try:
        #     # print("trying normalized", normalized["hgvs_g"])
        #     var_g = self.annotator.hp.parse_hgvs_variant(normalized["hgvs_g"])
        # except Exception:
        #     print(
        #         "Input HGVS/Mutalyzer string:",
        #         hgvs_or_mutalyzer_c_str,
        #         normalized["hgvs_g"],
        #     )
        # print(f"Error parsing HGVS g variant: {e}")

        #     try:
        #         # print("trying compliant", hgvs_compliant["hgvs_g"])
        #         var_c = self.annotator.hp.parse_hgvs_variant(
        #             hgvs_or_mutalyzer_c_str
        #         )
        #         var_g = self.annotator.am.c_to_g(var_c)
        #     except Exception as e:
        #         print(
        #             f"Error parsing HGVS g variant: {e}, input: {hgvs_compliant['hgvs_g']}"
        #         )
        #         # raise
        # if var_g is None:
        #     raise ValueError(
        #         f"Could not parse genomic HGVS for {hgvs_or_mutalyzer_c_str}"
        #     )
        # print("var_g", var_g, dir(var_g))
        # if (
        #     hasattr(var_g, "posedit")
        #     and var_g.posedit
        #     and hasattr(var_g.posedit, "pos")
        # ):
        #     print("post", var_g.posedit.pos, dir(var_g.posedit.pos))
        # try:
        #     var_c = self.annotator.hp.parse_hgvs_variant(normalized["hgvs_c"])
        #     # print("var_c", var_c, dir(var_c))
        # except hgvs.exceptions.HGVSParseError as e:
        #     print(f"Error parsing HGVS c variant: {e}")
        # try:
        #     var_c = self.annotator.hp.parse_hgvs_variant(
        #         hgvs_compliant["hgvs_c"]
        #     )
        #     # print("var_c compliant", var_c, dir(var_c))
        # except Exception as e:
        #     print(f"Error parsing HGVS c variant: {e}")
        #     raise
        # var_p = (
        #     self.hp.parse_hgvs_variant(hgvs_compliant["hgvs_p"])
        #     if hgvs_compliant["hgvs_p"]
        #     else None
        # )
        # var_r = (
        #     self.hp.parse_hgvs_variant(hgvs_compliant["hgvs_r"])
        #     if hgvs_compliant["hgvs_r"]
        #     else None
        # )
        # var_g = self.am.c_to_g(var_c)  # , alt_ac=transcript)
        # var_p = self._safe_c_to_p(var_c)
        # var_r = None
        return self._build_mutation(
            var_g=var_g,
            var_c=var_c,
            var_p=var_p,
            genomic_range=genomic_range,
            var_id=var_id,
            # var_r=var_r,
            # normalized=normalized,
            # hgvs_compliant=hgvs_compliant,
        )

    # def from_hgvs_g(self, hgvs_or_mutalyzer_g_str: str) -> Mutation:
    #     """Erzeugt eine Mutation-Instanz aus einem g.-HGVS-String."""
    #     hgvs_g_str, _, _ = self.mutalyzer_to_hgvs(
    #         hgvs_or_mutalyzer_g_str, translate_to="g"
    #     )
    #     var_g = self.hp.parse_hgvs_variant(hgvs_g_str)
    #     var_c = self._safe_g_to_c(var_g)
    #     var_p = self._safe_c_to_p(var_c) if var_c else None
    #     return self._build_mutation(
    #         var_c=var_c,
    #         var_g=var_g,
    #         var_p=var_p,
    #         hgvs_g_str=hgvs_g_str,
    #     )

    ############################################################################
    # Mutation Factory
    ############################################################################

    def _get_local_reference(
        self,
        genomic_pos_1: int,
        genomic_range: tuple[int, int] | None = None,
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
        flanks_around_position = (
            (
                genomic_pos_1 - self.flanking,
                genomic_pos_1 + self.flanking,
            )
            if self.flanking
            else None
        )
        flanks = genomic_range or self.genomic_flanks or flanks_around_position
        # priority: genomic_range > default flanks > local_cut > no cut
        # print(
        #     "flanks:",
        #     flanks,
        #     "genomic_range:",
        #     genomic_range,
        #     "self.genomic_flanks:",
        #     self.genomic_flanks,
        #     "flanks_around_position:",
        #     flanks_around_position,
        # )
        if flanks:
            flank_left = min(flanks)
            flank_right = max(flanks)
            frag_start = self._genomic_pos_to_local_0(
                flank_left
            )  # do not assume this is the smaller index
            frag_end = self._genomic_pos_to_local_0(flank_right)
            # print(
            #     "flank coordinates",
            #     flank_left,
            #     flank_right,
            #     "frag coordinates",
            #     frag_start,
            #     frag_end,
            # )
            local_pos_0 = self._genomic_pos_to_local_0(genomic_pos_1)
            frag_start, frag_end = min(frag_start, frag_end), max(
                frag_start, frag_end
            )
            # print(
            #     flank_left,
            #     flank_right,
            #     frag_start,
            #     frag_end,
            #     local_pos_0,
            #     genomic_pos_1,
            # )
            fragment = self.genomic.seq[frag_start : frag_end + 1]
            # print("Fragment", fragment, str(Seq(fragment).reverse_complement()))
            local_in_frag = local_pos_0 - frag_start
        else:
            fragment = self.genomic.seq
            local_in_frag = local_pos_0
        return fragment, local_in_frag

    def _apply_edit(
        self,
        ref_region: str,
        variant_pos: int,
        ref_len: int,
        edit,
    ) -> str:
        """
        Apply any hgvs Edit subclass to a reference region string.

        Parameters
        ----------
        ref_region  : the fetched reference sequence (with flanking)
        variant_pos : 0-based start of the variant within ref_region
        ref_len     : length of the ref allele (end = variant_pos + ref_len)
        edit        : any hgvs.edit.Edit subclass instance

        Supported edit types:
        ---------------------
        - NARefAlt: SNV (ref=1, alt=1), Deletion (ref>0, alt=""),
                    Insertion (ref="", alt>0),
                    DelIns (ref>0, alt>0, different lengths)
        - Dup: Tandem duplication of reference sequence
        - Inv: Inversion (reverse complement) of reference sequence
        - Repeat: Repeat reference sequence N times
        """
        before = ref_region[:variant_pos]
        ref_seq = ref_region[variant_pos : variant_pos + ref_len]
        after = ref_region[variant_pos + ref_len :]
        # print(
        #     variant_pos,
        #     ref_len,
        #     ref_seq,
        #     before,
        #     after,
        #     Seq(before).reverse_complement(),
        #     Seq(after).reverse_complement(),
        # )
        match edit:
            case hgvs.edit.NARefAlt():
                # SNV, Del, Ins, Delins — standard ref→alt substitution
                # DelIns: both ref and alt present, replaces ref_len with alt
                alt = edit.alt or ""
                # print("ALT", alt)
                return before + alt + after

            case hgvs.edit.Dup():
                # g.100_102dup → duplicates the ref_seq in tandem
                # print("DUP of", ref_seq)
                return before + ref_seq + ref_seq + after

            case hgvs.edit.Inv():
                # g.100_102inv → reverse complement of ref_seq
                inv_seq = str(Seq(ref_seq).reverse_complement())
                # print("INV of", ref_seq, "->", inv_seq)
                return before + inv_seq + after

            case hgvs.edit.Repeat():
                # g.100_102[4] → repeat ref_seq N times
                count = edit.seq.ref_n  # the repeat count
                # print("REPEAT of", ref_seq, "->", ref_seq * count)
                return before + (ref_seq * count) + after

            case _:
                raise NotImplementedError(
                    f"Edit type {type(edit).__name__} is not supported for "
                    f"sequence materialization"
                )

    def _ref_len_from_variant(self, var):
        edit = var.posedit.edit
        pos = var.posedit.pos

        # 1. Deletion (del)
        if isinstance(edit, hgvs.edit.NARefAlt) and edit.ref and not edit.alt:
            return len(edit.ref)

        # 2. Insertion (ins) ref_len = 0!
        if isinstance(edit, hgvs.edit.NARefAlt) and not edit.ref and edit.alt:
            return 0

        # 3. Substitution / delins
        if isinstance(edit, hgvs.edit.NARefAlt) and edit.ref:
            return len(edit.ref)

        # 4. Duplication
        if isinstance(edit, hgvs.edit.Dup):
            start = pos.start.base
            end = pos.end.base if pos.end and pos.end.base else start
            return end - start + 1

        # 5. Fallback
        start = pos.start.base
        end = pos.end.base if pos.end and pos.end.base else start
        return end - start + 1

    def _normalize_hgvs(
        self, var: SequenceVariant, hgvs: str | None
    ) -> str | None:
        try:
            normalized_hgvs = self.annotator.hn.normalize(var) if var else hgvs
        except Exception:
            normalized_hgvs = hgvs
        return normalized_hgvs

    def materialize(
        self,
        gvar: VariantG,
        genomic_range: tuple[int, int] | None = None,
    ) -> tuple[str, str, int]:
        # gvar.print_hgvs_details()
        ref_len = gvar.ref_len  # self._ref_len_from_variant(var)
        ref_region, local_variant_pos_0 = self._get_local_reference(
            gvar.genomic_start_1,
            genomic_range=genomic_range,
        )

        # ac = var.ac
        # pos = var.posedit.pos
        # g_start_1 = pos.start.base

        # if isinstance(edit, hgvs.edit.NARefAlt) and not edit.ref and edit.alt:
        #     # Insertion: Verwende die Position nach start
        #     g_start_1 = pos.end.base if pos.end else pos.start.base + 1
        # else:
        #     # Alle anderen: start ist die erste betroffene Base
        #     g_start_1 = pos.start.base

        # print(
        #     "We use",
        #     gvar.genomic_start_1,
        #     "with",
        #     gvar.edit_type,
        #     "to infer the local position",
        #     local_variant_pos_0,
        # )
        gvar.validate_ref(ref_region, local_variant_pos_0)
        # Optional ref validation where possible

        # self.__validate_ref(
        #     var.pe,
        #     ref_region,
        #     local_variant_pos_0,
        #     ref_len,
        #     var,
        #     g_start_1,
        # )

        mutated = self._apply_edit(
            ref_region, local_variant_pos_0, ref_len, gvar.edit
        )
        return mutated, ref_region, local_variant_pos_0

    def __infer_cds_pos_0_from_var_c(
        self, var_c: SequenceVariant
    ) -> int | None:
        cds_pos_0: int | None = None
        if var_c:
            c_start = var_c.posedit.pos.start
            # Nur für exonische Positionen ohne Intron-Offset
            if not (hasattr(c_start, "offset") and c_start.offset != 0):
                base = c_start.base
                if 1 <= base <= len(self.cds):
                    cds_pos_0 = self._c_pos_to_cds_0(base)
        return cds_pos_0

    def _build_mutation(
        self,
        var_g: SequenceVariant,
        var_c: SequenceVariant | None,
        var_p: SequenceVariant | None,
        genomic_range: tuple[int, int] | None = None,
        var_id: str | None = None,
        # var_r: SequenceVariant | None,
        # normalized: dict[str, str | None],
        # hgvs_compliant: dict[str, str | None],
    ) -> Mutation:

        # HGVS Strings
        # hgvs_g = hgvs_g_str or str(var_g)
        # hgvs_c = hgvs_c_str or (str(var_c) if var_c else None)
        # hgvs_p = str(var_p) if var_p else None

        # extract edit from g
        gvar = VariantG(var_g)
        # variant_info = self.get_g_variant_info(var_g)
        # g_ref, g_alt = gvar.ref, gvar.alt #self._extract_edit(var_g)
        edit_type = gvar.edit_type
        if var_c is not None:
            c_ref, c_alt = self._extract_edit(var_c)
        else:
            c_ref, c_alt = (
                (gvar.ref, gvar.alt)
                if self.strand == 1
                else (
                    reverse_complement(gvar.ref),
                    reverse_complement(gvar.alt),
                )
            )
        # mutation_type = _mut_type(g_ref, g_alt)
        # print(g_ref, g_alt, mutation_type)
        # print("var_g_in", hgvs_compliant["hgvs_g"])
        # positions = self.get_g_variant_info(var_g) if var_g else None
        # g_start_1 = positions["start"] if positions else None
        # g_start_1 = var_g.posedit.pos.start.base
        # local_start_0 = self._genomic_pos_to_local_0(g_start_1)
        # print(
        #     "Where it wants to place it",
        #     self.genomic.seq[local_start_0 : local_start_0 + 20],
        #     self.genomic.seq[7673605 - 20 : local_start_0],
        # )
        # print(
        #     "this is the local start",
        #     local_start_0,
        #     g_start_1,
        #     self.genomic.seq[3435 - len(g_ref) : 3435],
        # )
        # print("g_ref", g_ref, "g_alt", g_alt)
        # extract edit from c (if available) and determine cds_pos
        # print("strand", self.strand)
        g_start_1 = gvar.genomic_start_1
        cds_pos_0 = self.__infer_cds_pos_0_from_var_c(var_c) if var_c else None
        # print(
        #     "aaaaa",
        #     g_ref,
        #     c_ref,
        #     "alt",
        #     g_alt,
        #     c_alt,
        #     reverse_complement(g_ref),
        #     reverse_complement(c_ref),
        # )
        # print("c_ref", c_ref, "c_alt", c_alt)

        # Mutation type

        # Codon-Info (only for exonic SNVs)
        ref_codon = alt_codon = ref_aa = alt_aa = None
        prot_pos: int | None = None
        if cds_pos_0 is not None and edit_type == "snv":
            ref_codon, alt_codon, ref_aa, alt_aa, prot_pos = _codon_info(
                self.cds, cds_pos_0, c_alt, c_ref
            )

        # ── Region-Typ ────────────────────────────────────────────────────
        region_type = self._region_type(var_c)

        # ── Sequenz mit Flanking ──────────────────────────────────────────
        seq, ref_region, local_variant_pos_0 = self.materialize(
            gvar, genomic_range
        )
        # seq, mutation_pos = self._build_seq_with_flanking(
        #     local_start_0, g_ref, g_alt
        # )
        # hgvs_g = normalized.get("hgvs_g")
        # hgvs_c = normalized.get("hgvs_c")
        # hgvs_p = normalized.get("hgvs_p")
        # hgvs_r = normalized.get("hgvs_r")

        # hgvs_g_normalized = self._normalize_hgvs(
        #     var_g, hgvs_compliant["hgvs_g"]
        # )
        # hgvs_c_normalized = self._normalize_hgvs(
        #     var_c, hgvs_compliant["hgvs_c"]
        # )
        # hgvs_p_normalized = self._normalize_hgvs(
        #     var_p, hgvs_compliant["hgvs_p"]
        # )
        # hgvs_r_normalized = self._normalize_hgvs(
        #     var_r, hgvs_compliant["hgvs_r"]
        # )
        hgvs_g = self.annotator.to_hgvs(var_g) if var_g else None
        hgvs_c = self.annotator.to_hgvs(var_c) if var_c else None
        hgvs_p = self.annotator.to_hgvs(var_p) if var_p else None
        hgvs_r = self.annotator.hgvs_c_to_r(hgvs_c) if hgvs_c else None
        seq_id = var_id or f"{hgvs_c or hgvs_g}"
        # print("before mutation", self.strand, type(seq))
        if self.strand == -1:
            mutated = str(Seq(seq).reverse_complement())
            local_variant_pos_0 = (
                len(mutated) - local_variant_pos_0 - gvar.ref_len
            )
        else:
            mutated = seq

        return Mutation(
            seq_id=seq_id,
            genomic=seq,
            coding=mutated,
            genomic_id=self.genomic.id,
            chromosome=self.chromosome,
            strand=self.strand,
            genomic_start=self.gene_start_1,
            genomic_end=self.gene_end_1,
            region_type=region_type,
            mutation_pos=local_variant_pos_0,
            genomic_pos=g_start_1,  # 1-based
            cds_pos=cds_pos_0,  # 0-based
            ref=gvar.ref,
            alt=gvar.alt,
            mutation_type=gvar.edit_type,
            ref_codon=ref_codon,
            alt_codon=alt_codon,
            ref_aa=ref_aa,
            alt_aa=alt_aa,
            chromosome_ac=self.chrom_ac or "",
            transcript_ac=self.transcript_ac or "",
            protein_ac=self.protein_ac or "",
            prot_pos=prot_pos,
            hgvs_g=hgvs_g,
            hgvs_c=hgvs_c,
            hgvs_p=hgvs_p,
            hgvs_r=hgvs_r,
            refseq=ref_region,
        )
