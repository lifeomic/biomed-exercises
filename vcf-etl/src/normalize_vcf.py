import logging
import pandas as pd
import re


def read_vcf(vcf_in: str):
    """Read VCF file as pandas DataFrame"""
    # read meta info lines and header line
    with open(vcf_in) as file:
        meta_lines = []
        header_line = []
        for line in file:
            # store meta info lines in a list (these start with "##")
            if line.startswith("##"):
                meta_lines.append(line.rstrip())
            # store header line in a list (this starts with "#CHROM")
            elif line.startswith("#CHROM"):
                header_line.append(line.rstrip())
        # store meta lines as a string
        meta_lines = "\n".join(meta_lines)
        # generate list from header
        header_line = header_line[0].split("\t")
    # read vcf data into a DataFrame and skip meta info and header lines
    vcf_data = pd.read_csv(vcf_in, names=header_line, sep="\t", comment="#", on_bad_lines="warn")
    
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
    def match_af(info):
        """use a regular expression to match specific allele frquency i.e. AF={actual value}"""
        # split info into a list
        info_list = info.split(";")
        # setup the regular expression pattern
        string_pattern = "^AF"
        regex_pattern = re.compile(string_pattern)
        # match info list with regex pattern
        af_list = list(filter(regex_pattern.match, info_list))
        # convert allele frequency to a string or "." if there was no match for allele frequency
        allele_frequency = "".join(af_list) if af_list else "."

        return allele_frequency


    vcf_data["INFO"] = vcf_data["INFO"].apply(match_af)

    return vcf_data


def reformat_format_gt_ad_dp(vcf_data):
    """The FORMAT column should only include the following values 'GT:AD:DP'"""
    # set values that are allowed for FORMAT
    allowed_format_values = ["GT", "AD", "DP"]

    vcf_data["FORMAT"] = vcf_data["FORMAT"].apply(lambda x: ":".join([form for form in x.split(":") if form in allowed_format_values]))

    return vcf_data


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
    """normalize vcf"""
    # read vcf
    vcf_meta, vcf_data = read_vcf(vcf_in)

    # remove records where FILTER is not PASS or LowGQX
    vcf_data = remove_filter_not_pass_lowgqx(vcf_data=vcf_data)
    # remove variant rows with no alternate allele
    vcf_data = remove_no_alternate_alleles(vcf_data=vcf_data)
    # reformat INFO to only contain specific allele frequency
    vcf_data = reformat_info_af_only(vcf_data=vcf_data)
    # reformat FORMAT to only include GT:AD:AP
    vcf_data = reformat_format_gt_ad_dp(vcf_data=vcf_data)

    # write vcf
    write_vcf(vcf_meta_info=vcf_meta, vcf_data=vcf_data, vcf_out=vcf_out)
