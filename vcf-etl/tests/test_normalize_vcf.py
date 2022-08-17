from src.normalize_vcf import normalize_vcf
from src.normalize_vcf import read_vcf
from src.normalize_vcf import remove_filter_not_pass_lowgqx
from src.normalize_vcf import reformat_info_af_only
from src.normalize_vcf import reformat_format_gt_ad_dp
from src.normalize_vcf import remove_no_alternate_alleles
from src.normalize_vcf import write_vcf


def test_read_vcf():
    """read_vcf"""
    pass


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
