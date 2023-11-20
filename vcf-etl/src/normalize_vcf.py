import logging
import pandas as pd
import tempfile
import re

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('normalize_vcf')


def read_vcf(vcf_in: str) -> (pd.DataFrame, str):
    """
    Read VCF into pandas dataframe. Preserve VCF header.
    """
    vcf_header = ""
    lines = list()

    with open(vcf_in) as f:
        for line in f:
            if line.startswith('##'):
                vcf_header += line
            else:
                lines.append(line)

    headers = lines[0].split('\t')
    headers[-1] = headers[-1].strip()

    df = pd.DataFrame(columns=headers)

    for i, line in enumerate(lines[1:]):
        df.loc[i] = line.split('\t')[:10]

    return df, vcf_header


def write_vcf(vcf_file_meta_info: str, df: pd.DataFrame, vcf_out: str) -> None:
    """
    Write VCF to file (vcf_out).
    """
    with tempfile.TemporaryDirectory() as temp:
        df.to_csv(temp + 'data.vcf', sep='\t', index=False)

    # append the header from the original VCF file, to the transformed VCF file, and write to the specified vcf_out
    with open(temp + 'data.vcf') as f, open(vcf_out, 'a+') as file_out:
        file_out.write(vcf_file_meta_info)

        for line in f:
            file_out.write(line)


def filter_Q1(df: pd.DataFrame) -> pd.DataFrame:
    return df[df['FILTER'].isin(['PASS', 'LowGQX'])]


def filter_Q2(df: pd.DataFrame) -> pd.DataFrame:
    df.loc[:, 'INFO'] = df['INFO'].apply(lambda x: x.split(';'))
    
    filter_pattern = ["AF1000G"]

    df.loc[:, 'INFO'] = filter_helper(df['INFO'], filter_pattern)

    return df


def filter_helper(values: list, filter_pattern: list) -> list:
    def filter_elements(element):
        return any(pattern in element for pattern in filter_pattern)


    # Assign INFO as '.' for records with missing allele frequency
    filtered_column = [
        ";".join([element for element in value if filter_elements(element)]) or "."
        for value in values
    ]

    return filtered_column


def filter_Q3(df: pd.DataFrame) -> pd.DataFrame:
    desired_tags = ["GT", "AD", "DP"]

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


def filter_Q4(df: pd.DataFrame) -> pd.DataFrame:
    return df[df['ALT'] != '.']


def normalize_vcf(vcf_in: str, vcf_out: str) -> None:
    """
    Normalize an input VCF file (vcf_in) according to the following.

    Q1. Remove any records that don't have a FILTER of "PASS" or "LowGQX".
    Q2. The INFO column should only contain this variant's specific allele frequency, i.e. "AF={actual value}"
        Assign INFO as '.' to records missing an allele frequency.
    Q3. The FORMAT column should only include the following values "GT:AD:DP"
    Q4. The final file should only have variant rows with an alternate allele
    """

    logger.info(f' Reading input VCF file: {vcf_in}')
    _df, vcf_header = read_vcf(vcf_in)
    assert _df.shape == (37,10)

    logger.info(' Retaining records with PASS or LowGQX.')
    _df_Q1 = filter_Q1(_df)
    assert set(_df_Q1['FILTER']) == set(['PASS', 'LowGQX'])

    logger.info(' Filtering the INFO column for variant\'s allele frequency.')
    _df_Q2 = filter_Q2(_df_Q1)

    AF_match = r'AF1000G=0\.\d{6}'
    assert _df_Q2['INFO'].str.contains(f'^.$|{AF_match}').all() == True


    logger.info(' Filtering the FORMAT column for the following values: "GT:AD:DP".')
    _df_Q3 = filter_Q3(_df_Q2)
    
    for sub_list in _df_Q3['FORMAT']:
        assert all(x in ['GT','AD','DP'] for x in sub_list.split(':'))

    # check that FORMAT column and NA12878 column contain the same number of values
    for i, j in zip(_df_Q3['FORMAT'], _df_Q3['NA12878']):
        assert (len(list(i.split(':'))) == len(list(j.split(':'))))
        

    logger.info(' Retaining records that have an alternative allele.')
    _df_Q4 = filter_Q4(_df_Q3)
    # print(_df_Q4['ALT'])
    for i, j in zip(_df_Q4['REF'], _df_Q4['ALT']):
        assert i != j
        assert j != '.'

    logger.info(f' Writing the normalized VCF to file: {vcf_out}.')
    write_vcf(vcf_header, _df_Q4, vcf_out)


# def main():
#     normalize_vcf('../data/NA12878_cod_ex.genome.vcf', 'test.vcf')


# if __name__ == "__main__":
#     main()
