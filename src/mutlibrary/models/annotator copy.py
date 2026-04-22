from dataclasses import dataclass

from mavehgvs.variant import Variant as HgvsVariant  # type: ignore[import]

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


@dataclass
class Mutation:
    seq_id: str
    genomic_id: str
    region_type: str
    mutation_pos: int
    genomic_pos: int
    cds_pos: int | None
    prot_pos: int | None
    ref: str
    alt: str
    mutation_type: str
    hgvs_g: str
    hgvs_c: str
    hgvs_p: str
    seq: str
    ref_aa: str | None
    alt_aa: str | None


def _build_g_hgvs(chrom_ac: str, genomic_pos: int, ref: str, alt: str) -> str:
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
        raw = f"g.{g1}del" if ref_len == 1 else f"g.{g1}_{g1 + ref_len - 1}del"
    else:
        # MNV / delins
        raw = (
            f"g.{g1}delins{alt}"
            if ref_len == 1
            else f"g.{g1}_{g1 + ref_len - 1}delins{alt}"
        )

    return f"{chrom_ac}:{HgvsVariant(raw)}"


class MutationAnnotator:

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
        # No UTA database connection required — mavehgvs is purely string-based.
        self.CHROM_TO_NC = (
            self._CHROM_TO_NC_GRCH37
            if assembly in ("GRCh37", "hg19")
            else self._CHROM_TO_NC_GRCH38
        )

    def _annotate_variant(
        self,
        genomic_id: str,
        chrom: str | int,
        genomic_pos: int,  # 0-based
        region_type: str,
        ref: str,
        alt: str,
        cds_pos: int | None = None,
        ref_codon: str | None = None,
        alt_codon: str | None = None,
        ref_aa: str | None = None,
        alt_aa: str | None = None,
    ) -> tuple[str, str, str, int | None]:
        use_cds = cds_pos is not None
        chrom_ac = self.CHROM_TO_NC[str(chrom)]
        c_pos1 = (cds_pos + 1) if use_cds else None  # type: ignore[operator]

        # ── 1. g. ─────────────────────────────────────────────────────────
        hgvs_g = _build_g_hgvs(chrom_ac, genomic_pos, ref, alt)

        # ── 2. c. ─────────────────────────────────────────────────────────
        hgvs_c = ""
        if use_cds and genomic_id:
            if not ref and alt:
                # insertion: 1-based adjacent positions
                raw_c = f"c.{c_pos1}_{c_pos1 + 1}ins{alt}"  # type: ignore[operator]
            elif ref and not alt:
                if len(ref) == 1:
                    raw_c = f"c.{c_pos1}del"
                else:
                    end = c_pos1 + len(ref) - 1  # type: ignore[operator]
                    raw_c = f"c.{c_pos1}_{end}del"
            elif ref_codon and alt_codon:
                diffs = [i for i in range(3) if ref_codon[i] != alt_codon[i]]
                if len(diffs) == 1:
                    pos_nt = cds_pos + diffs[0] + 1  # type: ignore[operator]
                    raw_c = (
                        f"c.{pos_nt}"
                        f"{ref_codon[diffs[0]]}>{alt_codon[diffs[0]]}"
                    )
                else:
                    start_nt = cds_pos + diffs[0] + 1  # type: ignore[operator]
                    end_nt = cds_pos + diffs[-1] + 1  # type: ignore[operator]
                    alt_sub = alt_codon[diffs[0] : diffs[-1] + 1]
                    raw_c = f"c.{start_nt}_{end_nt}delins{alt_sub}"
            else:
                raw_c = f"c.{c_pos1}{ref}>{alt}"
            hgvs_c = f"{genomic_id}:{HgvsVariant(raw_c)}"

        # ── 3. p. ─────────────────────────────────────────────────────────
        hgvs_p = ""
        if region_type == "exon" and ref_aa and alt_aa and use_cds:
            codon_index = cds_pos // 3 + 1  # type: ignore[operator]
            ref_aa_3 = _AA_1TO3.get(ref_aa, ref_aa)
            alt_aa_3 = _AA_1TO3.get(alt_aa, alt_aa)
            raw_p = f"p.{ref_aa_3}{codon_index}{alt_aa_3}"
            hgvs_p = f"{genomic_id}:{HgvsVariant(raw_p)}"

        prot_pos = (cds_pos // 3 + 1) if use_cds else None  # type: ignore[operator]
        return hgvs_g, hgvs_c, hgvs_p, prot_pos

    def annotate(
        self,
        seq_id: str,
        genomic_id: str,
        chrom: str | int,
        mutation_pos: int,
        genomic_pos: int,  # 0-based
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
        hgvs_g, hgvs_c, hgvs_p, prot_pos = self._annotate_variant(
            genomic_id=genomic_id,
            chrom=chrom,
            genomic_pos=genomic_pos,
            region_type=region_type,
            ref=ref,
            alt=alt,
            cds_pos=cds_pos,
            ref_codon=ref_codon,
            alt_codon=alt_codon,
            ref_aa=ref_aa,
            alt_aa=alt_aa,
        )
        return Mutation(
            seq_id=seq_id,
            genomic_id=genomic_id,
            region_type=region_type,
            mutation_pos=mutation_pos,
            genomic_pos=genomic_pos,
            cds_pos=cds_pos,
            prot_pos=prot_pos,
            ref=ref,
            alt=alt,
            mutation_type=mut_type,
            hgvs_g=hgvs_g,
            hgvs_c=hgvs_c,
            hgvs_p=hgvs_p,
            seq=seq,
            ref_aa=ref_aa,
            alt_aa=alt_aa,
        )
