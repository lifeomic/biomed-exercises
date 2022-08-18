import pytest
import pandas as pd
from src.normalize_vcf import normalize_vcf
from src.normalize_vcf import read_vcf
from src.normalize_vcf import remove_filter_not_pass_lowgqx
from src.normalize_vcf import reformat_info_af_only
from src.normalize_vcf import reformat_format_gt_ad_dp
from src.normalize_vcf import remove_no_alternate_alleles
from src.normalize_vcf import write_vcf


@pytest.fixture
def input_minimal_vcf():
    """create a minimal test vcf dataset that is the same as data/minimal.vcf"""
    # set metadata lines
    meta_lines = "##fileformat=VCFv4.1"
    # construct DataFrame of vcf data
    column_names = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NA12878"]
    vcf_data = pd.DataFrame(columns = column_names)
    vcf_data.loc[0] = ["chr1", 10496, ".", "C", ".", ".", "PASS", "END=10514;BLOCKAVG_min30p3a", "GT:GQX:DP:DPF", "0/0:15:6:0"]

    return meta_lines, vcf_data


def test_read_vcf(input_minimal_vcf):
    """read_vcf"""
    # set expected values for reading in data/minimal.vcf
    expected_meta_lines, expected_vcf_data = input_minimal_vcf
    # run read_vcf for minimal test set and capture output
    actual_meta_lines, actual_vcf_data = read_vcf("data/minimal.vcf")

    # check that the outputs of read_vcf match the expected
    assert expected_meta_lines == actual_meta_lines
    assert pd.testing.assert_frame_equal(expected_vcf_data, actual_vcf_data)


def test_remove_filter_not_pass_lowgqx():
    """remove_filter_not_pass_lowgqx"""
    pass


def test_reformat_info_af_only():
    """reformat_info_af_only"""
    pass


def test_reformat_format_gt_ad_dp():
    """reformat_format_gt_ad_dp"""
    pass


def test_remove_no_alternate_alleles():
    """remove_no_alternate_alleles"""
    pass


def test_write_vcf():
    """write_vcf"""
    pass


def test_normalize_vcf():
    actual = normalize_vcf("file-in", "file-out")
    assert actual == "Hello World"
