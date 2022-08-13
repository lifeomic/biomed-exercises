from src.normalize_vcf import vcf_reader, filter_filter, alt_filter, info_filter, format_filter, vcf_writer, normalize_vcf
import subprocess
import re
import tempfile
import pytest


@pytest.fixture
def dataframe():
    return vcf_reader("data/NA12878_cod_ex.genome.vcf")


def test_vcf_reader(dataframe):
    # Check that .vcf was read in properly
    assert dataframe.shape == (37,10)


def test_filter_filter(dataframe):
    # Check that FILTER column only contains PASS or LowGQX values
    filter_df = filter_filter(dataframe)
    for i in filter_df['FILTER']:
        assert i == 'PASS' or i == 'LowGQX'


def test_alt_filter(dataframe):
    # Check that ALT column only contains alternative alleles from the REF column
    alt_df = alt_filter(dataframe)
    for i, j in zip(alt_df['REF'], alt_df['ALT']):
        assert i != j
        assert j != '.'


def test_info_filter(dataframe):
    # Check that INFO column only contains AF values or a period value
    info_df = info_filter(dataframe)
    regex = r'AF1000G=0\.\d{6}'
    for i in info_df['INFO']:
        res = False
        if re.match(regex, i) or i == '.':
            res = True
        assert res == True


def test_format_filter(dataframe):
    # Check that FORMAT column only contains desired values (GT, AD, and DP)
    format_df = format_filter(dataframe)
    for i in format_df['FORMAT']:
        for j in i.split(':'):
            assert j == 'GT' or j == 'AD' or j == 'DP'

    # Check that FORMAT column and VARIANT column have equal amount of items
    for i in zip(format_df['FORMAT'], format_df.iloc[:, -1]):
        assert (len(list(i[0].split(':'))) == len(list(i[1].split(':'))))


def test_normalize_vcf():
    # Check that .vcf was properly written out
    with tempfile.TemporaryDirectory() as temp:
        normalize_vcf('data/NA12878_cod_ex.genome.vcf', temp + 'test.vcf')
        vcf_out_line_count = subprocess.check_output('wc -l '+temp+'test.vcf', shell=True)
        assert vcf_out_line_count.strip().decode('utf-8')[0:3] == '105'
