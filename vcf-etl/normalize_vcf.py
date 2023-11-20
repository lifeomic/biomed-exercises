import logging
import pandas as pd
import tempfile



def read_vcf(vcf_in: str):
    vcf_file_meta_info = ""

    lines = list()

    with open(vcf_in) as f:

        for line in f:
            # preserve meta information lines
            if line.startswith('##'):
                vcf_file_meta_info += line

            else:
                lines.append(line)

    headers = lines[0].split('\t')

    headers[-1] = headers[-1].strip()


    df = pd.DataFrame(columns=headers)

    for i, line in enumerate(lines[1:]):
        df.loc[i] = line.split('\t')[:10]

    return df, vcf_file_meta_info


def write_vcf(vcf_file_meta_info: str, df: pd.DataFrame, vcf_out: str):
    """
    Write VCF data file.
    """
    with tempfile.TemporaryDirectory() as temp:
        df.to_csv(temp + 'data.vcf', sep='\t', index=False)

    # append the header from the original VCF file, to the transformed VCF file, and write to the specified vcf_out
    with open(temp + 'data.vcf') as f, open(vcf_out, 'a+') as file_out:
        file_out.write(vcf_file_meta_info)

        for line in f:
            file_out.write(line)



def filter_vcf(df: pd.DataFrame):
    """
    filter vcf
    """
    filter_df = filter_status_filter(df)

    info_filter_df = info_filter(filter_df)

    format_filter_df = format_filter(info_filter_df)

    alt_filter_df = alt_allele_filter(info_filter_df)

    # print(alt_filter_df)

    return alt_filter_df


def filter_status_filter(df):
    """
    Remove any records that don't have a FILTER of "PASS" or "LowGQX"
    """
    df = df[df['FILTER'].isin(['PASS', 'LowGQX'])]

    return df


def info_filter(df: pd.DataFrame):
    """
    Tranfrom INFO to only contain the variant's specific allele frequency. 
    Assign INFO as '.' to records missing an allele frequency.
    """
    # split INFO by semi colon
    df.loc[:, 'INFO'] = df['INFO'].apply(lambda x: x.split(';'))
    
    filter_pattern = ["AF1000G"]

    df.loc[:, 'INFO'] = filter_helper(df['INFO'], filter_pattern)

    # print(df['INFO'])

    return df


def filter_helper(values: list, filter_pattern: list):
    """
    Filters a column based on the specified filter_pattern
    """
    def filter_elements(element):
        return any(pattern in element for pattern in filter_pattern)


    # Assign INFO as '.' for records with missing allele frequency
    filtered_column = [
        ";".join([element for element in value if filter_elements(element)]) or "."
        for value in values
    ]

    return filtered_column


def format_filter(df: pd.DataFrame):
    """
    Transform FORMAT to only include the following values "GT:AD:DP"
    """
    desired_tags = ["GT", "DP", "AD"]

    # Apply the filtering logic to the DataFrame
    filtered_tags = []
    filtered_values = []

    for index, row in df.iterrows():
        tags_list = row["FORMAT"].split(":")
        values_list = row["NA12878"].split(":")

        filtered_tags_row = [tag for tag in desired_tags if tag in tags_list]
        filtered_values_row = [values_list[tags_list.index(tag)] for tag in filtered_tags_row]

        filtered_tags.append(":".join(filtered_tags_row))
        filtered_values.append(":".join(filtered_values_row))

    # Update the original DataFrame with filtered tags and values
    df.loc[:, "FORMAT"] = filtered_tags
    df.loc[:, "NA12878"] = filtered_values

    return df


def alt_allele_filter(df: pd.DataFrame):
    """
    Remove variant rows with that don't have an alternate allele
    """
    return df[df['ALT'] != '.']


def normalize_vcf(vcf_in: str, vcf_out: str):

    df, vcf_file_meta_info = read_vcf(vcf_in)

    # filtered_df = filter_vcf(df)
    # print(filtered_df)

    filter_df = filter_status_filter(df)

    info_filter_df = info_filter(filter_df)

    format_filter_df = format_filter(info_filter_df)

    alt_filter_df = alt_allele_filter(info_filter_df)

    write_vcf(vcf_file_meta_info, alt_filter_df, vcf_out)


# def main():
#     normalize_vcf('../data/NA12878_cod_ex.genome.vcf', 'test.vcf')


# if __name__ == "__main__":
#     main()
