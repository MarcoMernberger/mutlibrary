import sys
from pathlib import Path

import pandas as pd  # type: ignore[import]
import pytest  # type: ignore[import]

from mutlibrary.models.mutator import MutatorHGVS  # type: ignore[import]

root = Path(__file__).parent.parent
sys.path.append(str(root / "src"))


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
def real_fasta():
    return Path(__file__).parent / "data" / "test_real.fasta"


@pytest.fixture(scope="session")
def real_region():
    return Path(__file__).parent / "data" / "real_regions.json"


@pytest.fixture(scope="session")
def test_region():
    return Path(__file__).parent / "data" / "regions.json"


@pytest.fixture(scope="session")
def annotation():
    return Path(__file__).parent / "data" / "annotations.json"


@pytest.fixture(scope="session")
def clinvar_cases():
    df = pd.read_csv(
        Path(__file__).parent / "data" / "test_cases_clinvar.tsv", sep="\t"
    )
    df = df.where(df.notna(), None)
    return df


@pytest.fixture
def real_records(real_fasta):
    major = MutatorHGVS()
    return major.read_records(real_fasta)


@pytest.fixture
def real_sample(real_records):
    return real_records["TP53-201_Ex9_realtest"]


@pytest.fixture
def real_cds(real_records):
    print(real_records.keys())
    return real_records["TP53-201_cds_protein_coding"]


@pytest.fixture
def real_genomic(real_records):
    return real_records["17"]


@pytest.fixture
def real_regions(real_region, real_records):
    major = MutatorHGVS()
    return major.region_definitions(real_region, real_records)
