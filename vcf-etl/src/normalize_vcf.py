import logging


def read_vcf(vcf_in: str):
    """Read VCF file as pandas DataFrame"""
    pass


def remove_filter_not_pass_lowgqx(vcf_data):
    """Remove any records that don't have a FILTER of 'PASS' or 'LowGQX'"""
    pass


def reformat_info_af_only(vcf_data):
    """The INFO column should only contain this variant's specific allele frequency, i.e. 'AF={actual value}'"""
    pass


def reformat_format_gt_ad_dp(vcf_data):
    """The FORMAT column should only include the following values 'GT:AD:DP'"""
    pass


def remove_no_alternate_alleles(vcf_data):
    """The final file should only have variant rows with an alternate allele."""
    pass


def write_vcf(vcf_meta_info, vcf_data, vcf_out: str):
    """Write vcf to output file"""
    pass


def normalize_vcf(vcf_in: str, vcf_out: str):
    print("Hello World")
    return "Hello World"
