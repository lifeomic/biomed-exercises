import logging
import pandas as pd
import tempfile
import re


def vcf_reader(vcf_in):
    with open(vcf_in) as f:
        lines = [line for line in f if not line.startswith('##')]

    headers = lines[0].split('\t')
    headers[-1] = headers[-1].strip()

    df = pd.DataFrame(columns=headers)

    for line in lines[1:]:
        # Ensure proper number of columns per entry / removes trailing empty columns
        df.loc[len(df.index)] = line.split('\t')[:10]

    return df


def filter_filter(dataframe):
    new_dataframe_1 = dataframe[dataframe.FILTER.isin(['PASS', 'LowGQX'])]
    return new_dataframe_1


def alt_filter(dataframe):
    new_dataframe_2 = dataframe.drop(dataframe[dataframe['ALT'] == '.'].index)
    return new_dataframe_2


def info_filter(dataframe):
    # Replacing non-AF containing entries with a period value
    dataframe.INFO = dataframe.INFO.apply(lambda x: '.' if 'AF1000G' not in x else x)

    # Removing all text before AF value
    front_end = re.compile(r'^.*AF1000G')
    dataframe.INFO = dataframe.INFO.str.replace(front_end, 'AF1000G')

    # Removing all text after AF value
    dataframe.INFO = dataframe.INFO.str[:16]

    new_dataframe_3 = dataframe
    return new_dataframe_3


def format_filter(dataframe):
    # Isolating FORMAT column and corresponding VARIANT column values
    f_list = list(dataframe['FORMAT'])
    na_list = list(dataframe.iloc[:, -1])
    zipped = zip(f_list, na_list)

    new_f_list = []
    new_na_list = []
    for i in zipped:
        # Pruning FORMAT to only include GT, AD, and DP values.
        # Pruning VARIANT column to only include the corresponding values.
        zip_dict = dict(zip(i[0].split(':'), i[1].split(':')))
        new_dict = {k: v for k, v in zip_dict.items() if k in ['GT', 'AD', 'DP']}

        key_string = ':'.join(list(new_dict.keys()))
        new_f_list.append(key_string)

        value_string = ':'.join(list(new_dict.values()))
        new_na_list.append(value_string)

    # Replacing FORMAT and VARIANT column with new pruned lists
    dataframe['FORMAT'] = new_f_list
    dataframe.iloc[:, -1] = new_na_list

    new_dataframe_4 = dataframe
    return new_dataframe_4


def vcf_writer(dataframe, vcf_in, vcf_out):
    with tempfile.TemporaryDirectory() as temp:
        dataframe.to_csv(temp + 'data.vcf', sep='\t', index=False)

    # Grabbing header from vcf_in, data from dataframe, and writing to specified vcf_out
    with open(vcf_in) as file_one, open(temp + 'data.vcf') as file_two, open(vcf_out, 'a+') as file_out:
        for line in file_one:
            if line.startswith('##'):
                file_out.write(line)
        for line in file_two:
            file_out.write(line)


def normalize_vcf(vcf_in: str, vcf_out: str):
    df = vcf_reader(vcf_in)

    filter_df = filter_filter(df)

    alt_df = alt_filter(filter_df)

    info_df = info_filter(alt_df)

    format_df = format_filter(info_df)

    vcf_writer(format_df, vcf_in, vcf_out)
