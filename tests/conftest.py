import sys
from pathlib import Path

import pandas as pd
import pytest

root = Path(__file__).parent.parent
sys.path.append(str(root / "src"))

from mutlibrary.models.mutator import MutatorHGVS


@pytest.fixture(scope="session")
def test_fasta():
    return Path(__file__).parent / "data" / "test.fasta"


@pytest.fixture(scope="session")
def tp53_201_genomic():
    infile = (
        Path(__file__).parent / "data" / "tp53_ENST00000269305_genomic.fasta"
    )
    records = MutatorHGVS().read_records(infile)
    return records["17"]


@pytest.fixture(scope="session")
def tp53_201_cds():
    infile = Path(__file__).parent / "data" / "tp53_ENST00000269305_cds.fasta"
    records = MutatorHGVS().read_records(infile)
    return records["TP53-201_cds_protein_coding"]


@pytest.fixture(scope="session")
def test_region():
    return Path(__file__).parent / "data" / "regions.json"


@pytest.fixture(scope="session")
def test_annotation():
    return Path(__file__).parent / "data" / "annotations.json"


@pytest.fixture(scope="session")
def clinvar_cases():
    df = pd.read_csv(
        Path(__file__).parent / "data" / "test_cases_clinvar.tsv", sep="\t"
    )
    df = df.where(df.notna(), None)
    return df
