import logging
import pandas as pd


def read_vcf(vcf_in: str):
    """Read VCF file as pandas DataFrame"""
    # read meta info lines and join by new lines and header line
    with open(vcf_in) as file:
        meta_lines = []
        header_line = []
        for line in file:
            if line.startswith("##"):
                meta_lines.append(line.rstrip())
            elif line.startswith("#CHROM"):
                header_line.append(line.rstrip())
        meta_lines = "\n".join(meta_lines)
        header_line = header_line[0].split("\t")
    # read vcf data into a DataFrame
    vcf_data = pd.read_csv(vcf_in, names=header_line, sep="\t", comment="#")
    
    return meta_lines, vcf_data


def remove_filter_not_pass_lowgqx(vcf_data):
    """Remove any records that don't have a FILTER of 'PASS' or 'LowGQX'"""
    # set values that are allowed for FILTER
    allowed_filter_values = ["PASS", "LowGQX"]

    # remove rows not in allowed FILTER values
    vcf_data_filtered = vcf_data.loc[vcf_data["FILTER"].isin(allowed_filter_values)]

    return vcf_data_filtered


def reformat_info_af_only(vcf_data):
    """The INFO column should only contain this variant's specific allele frequency, i.e. 'AF={actual value}'"""
    pass


def reformat_format_gt_ad_dp(vcf_data):
    """The FORMAT column should only include the following values 'GT:AD:DP'"""
    pass


def remove_no_alternate_alleles(vcf_data):
    """The final file should only have variant rows with an alternate allele."""
    # remove rows with no alternate alleles
    vcf_data_filtered = vcf_data.loc[vcf_data["ALT"] != "."]

    return vcf_data_filtered


def write_vcf(vcf_meta_info, vcf_data, vcf_out: str):
    """Write vcf to output file"""
    # open output file
    with open(vcf_out, "w") as file:
        # write meta info lines
        file.write(vcf_meta_info + "\n")
    # append header and vcf data to file
    vcf_data.to_csv(vcf_out, mode="a", index=False, sep="\t")


def normalize_vcf(vcf_in: str, vcf_out: str):
    pass
