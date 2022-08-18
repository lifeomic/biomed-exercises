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


@pytest.fixture
def input_vcf_df():
    """creates a fake vcf DataFrame"""
    # construct DataFrame of vcf data
    column_names = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NA12878"]
    vcf_data = pd.DataFrame(columns = column_names)
    # add variants that need no normalization
    vcf_data.loc[0] = ["chr1", 1, ".", "A", "C", ".", "PASS", "AF1000G=0.494209", "GT:AD:DP", "0/0:15:6:0"]
    vcf_data.loc[1] = ["chr2", 1, ".", "G", "A", ".", "LowGQX", ".", "GT:DP", "0/0:124:44:0"]
    # add variants that require normalization with each required change
    vcf_data.loc[2] = ["chr3", 1, ".", "T", ".", ".", "HighDPFRatio", "SNVHPOL=4;MQ=60;GMAF=G|0.1154;AF1000G=0.884585;phyloP=0.241", "GT:GQX:DP:DPF", "0/0:15:6:0"]
    vcf_data.loc[3] = ["chr4", 1, ".", "C", ".", ".", "LowGQX;HighDPFRatio", "END=29094578;BLOCKAVG_min30p3a", "GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL", "0/0:124:44:0"]

    return vcf_data


def test_read_vcf(input_minimal_vcf):
    """read_vcf"""
    # set expected values for reading in data/minimal.vcf
    expected_meta_lines, expected_vcf_data = input_minimal_vcf
    # run read_vcf for minimal test set and capture output
    actual_meta_lines, actual_vcf_data = read_vcf("data/minimal.vcf")

    # check that the outputs of read_vcf match the expected
    assert expected_meta_lines == actual_meta_lines
    assert pd.testing.assert_frame_equal(expected_vcf_data, actual_vcf_data)


def test_remove_filter_not_pass_lowgqx(input_vcf_df):
    """remove_filter_not_pass_lowgqx"""
    # run remove_filter_not_pass_lowgqx on the test set
    actual_vcf_df = remove_filter_not_pass_lowgqx(input_vcf_df)
    # expected and actual series of #CHROM values
    expected_chrom = pd.Series(["chr1", "chr2"], name="#CHROM")
    actual_chrom = actual_vcf_df["#CHROM"]

    # check that rows 3 and 4 are dropped
    assert len(actual_vcf_df) == 2
    pd.testing.assert_series_equal(expected_chrom, actual_chrom)


def test_reformat_info_af_only(input_vcf_df):
    """reformat_info_af_only"""
    # run reformat_info_af_only on the test set
    actual_vcf_df = reformat_info_af_only(input_vcf_df)
    # expected and actual series of INFO values
    expected_info = pd.Series(["AF1000G=0.494209", ".", "AF1000G=0.884585", "."], name="INFO")
    actual_info = actual_vcf_df["INFO"]

    # check that only specific allele frequencines are in INFO
    assert len(actual_vcf_df) == 4
    pd.testing.assert_series_equal(expected_info, actual_info)


def test_reformat_format_gt_ad_dp(input_vcf_df):
    """reformat_format_gt_ad_dp"""
    # run reformat_format_gt_ad_dp on test set
    actual_vcf_df = reformat_format_gt_ad_dp(input_vcf_df)
    # expected and actual series of FORMAT values
    expected_format = pd.Series(["GT:AD:DP", "GT:DP", "GT:DP", "GT:DP:AD"], name="FORMAT")
    actual_format = actual_vcf_df["FORMAT"]

    # check that only GT:AD:DP in FORMAT column
    assert len(actual_vcf_df) == 4
    pd.testing.assert_series_equal(expected_format, actual_format)


def test_remove_no_alternate_alleles():
    """remove_no_alternate_alleles"""
    pass


def test_write_vcf():
    """write_vcf"""
    pass


def test_normalize_vcf():
    actual = normalize_vcf("file-in", "file-out")
    assert actual == "Hello World"
