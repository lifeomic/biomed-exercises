from src.normalize_vcf import normalize_vcf, read_vcf, filter_Q1, filter_Q2, filter_Q3, filter_Q4

import tempfile
import pytest
import re

@pytest.fixture
def dataframe():
    return read_vcf("data/NA12878_cod_ex.genome.vcf")[0]


def test_vcf_reader(dataframe):
    """
    Check that .vcf was read in properly.
    """
    assert dataframe.shape == (37,10)



def test_filter_Q1(dataframe):
    """
    Check that FILTER column only contains PASS or LowGQX values.
    """
    df_Q1 = filter_Q1(dataframe)
    assert set(filter_df['FILTER']) == set(['PASS', 'LowGQX'])


def test_filter_Q2(dataframe):
    """
    Check that INFO column only contains AF values or a period value
    """
    df_Q2 = filter_Q2(dataframe)
    assert df_Q2['INFO'].str.contains(f'^.$|{AF_match}').all() == True


def test_filter_Q3(dataframe):
    """
    test the 3rd filter function.
    
    """
    df_Q3 = format_Q3(dataframe)

    # Check that FORMAT column only contains the desired values (GT, AD, DP)
    for sub_list in _df_Q3['FORMAT']:
        assert all(x in ['GT','AD','DP'] for x in sub_list.split(':'))

    # Check that FORMAT column and NA12878 column contain the same number of values
    for i, j in zip(df_Q3['FORMAT'], df_Q3['NA12878']):
        assert (len(list(i.split(':'))) == len(list(j.split(':'))))


def test_filter_Q4(dataframe):
    """
    Check that ALT column is not empty and not equal to the REF column
    """
    df_Q4 = filter_Q4(dataframe)
    for i, j in zip(df_Q4['REF'], df_Q4['ALT']):
        assert i != j
        assert j != '.'


def test_normalize_vcf():
    actual = normalize_vcf("file-in", "file-out")
    assert actual != "Hello World"
