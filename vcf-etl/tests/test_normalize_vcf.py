import pandas as pd
from src.normalize_vcf import (
    normalize_vcf,
    read_vcf,
    screen_filter,
    screen_alt_allele,
    get_AF,
    trim_info,
    trim_format,
    trim_format_rowfunc,
    get_target_field_indexes,
    extract_by_index,
    write_vcf,
    get_outfile,
    normalize_vcf
)
import pytest

filepath = "data/NA12878_cod_ex.genome.vcf"
filepath_multicol = "data/Multi_sample_cod_ex.genome.vcf"

def test_read_vcf():
    df, metadata = read_vcf(filepath)
    assert isinstance(df, pd.DataFrame)
    assert isinstance(metadata, str)


@pytest.fixture
def vcf_dataframe():
    df, _ = read_vcf(filepath)
    return df


# Fixture for testing functions that handle multi-sample VCF files
#@pytest.fixture
#def vcf_dataframe_multicol():
#    df, _ = read_vcf(filepath_multicol)
#    return df


def test_screen_filter(vcf_dataframe):
    filtered_df = screen_filter(vcf_dataframe)
    n_post_filter = len(filtered_df.index)
    assert (
        n_post_filter == 30
    ), "Incorrect number of variants removed from {} with 'FILTER == PASS or LowGQX'.".format(
        filepath
    )


def test_screen_alt_allele(vcf_dataframe):
    altered_df = screen_alt_allele(vcf_dataframe)
    n_post_altered = len(altered_df.index)
    assert (
        n_post_altered > 0
    ), "ALL variants removed from {} based on not having an alternate allele.".format(
        filepath
    )


# Mock DataFrame for testing trim_info function
@pytest.fixture
def mock_dataframe():
    data = {
        "INFO": ["AF1000G=0.25;SomeOtherInfo", "AF1000G=0.75;SomeMoreInfo", "SomeInfo"]
    }
    return pd.DataFrame(data)


def test_get_AF():
    # Test the get_AF function
    assert get_AF("AF1000G=0.25;SomeOtherInfo") == "AF1000G=0.25"
    assert get_AF("AF1000G=0.75;SomeMoreInfo") == "AF1000G=0.75"
    assert get_AF("SomeInfo") == "."

def test_trim_info(vcf_dataframe):
    trim_df = trim_info(vcf_dataframe)
    test_value = ''.join(trim_df['INFO'][1:5].tolist())
    expected_value  = '...AF1000G=0.494209'
    assert test_value == expected_value

# Mock DataFrame for testing trim_info function
@pytest.fixture
def mock_dataframe():
    data = {
        'INFO': ['SomeAnnotations;AF1000G=0.25', 'AF1000G=0.75;SomeMoreInfo', 'SomeInfo']
    }
    return pd.DataFrame(data)

def test_get_AF():
    # Test the get_AF function
    assert get_AF('AF1000G=0.25;SomeOtherInfo') == 'AF1000G=0.25'
    assert get_AF('AF1000G=0.75;SomeMoreInfo') == 'AF1000G=0.75'
    assert get_AF('SomeInfo') == '.'

def test_get_target_field_indexes():
    assert get_target_field_indexes('FOO:GT:BAR:AD:BAZ:DP') == [1,3,5]

def test_extract_by_index():
    assert extract_by_index('FOO:GT:BAR:AD:BAZ:DP', [1,3,5]) == 'GT:AD:DP'

