import logging
import pandas as pd
import numpy as np
import os

def read_vcf(file_path: str) -> (str, pd.DataFrame):

    if not os.path.isfile(file_path):
        raise FileNotFoundError
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
        header = []
        head_count = 0
        for line in lines:
            if line.startswith("##"):
                head_count +=1
                header.append(line)
            
    vcf_head = ''.join(header) 
    vcf_b = lines[head_count:]    

    vcf_body = pd.DataFrame([x.strip().split('\t') for x in vcf_b])
    vcf_body.columns = vcf_body.iloc[0]
    vcf_body = vcf_body.iloc[1:]

    return (vcf_head, vcf_body)

def write_vcf(file_path: str, vcf_header: str, df_body: pd.DataFrame) -> None:

    with open(file_path, 'w') as f:
        f.write(vcf_header)

    df_body.to_csv(file_path, sep="\t", mode='a', index=False)

    return

def filter_Q1(df_in: pd.DataFrame) -> pd.DataFrame:
    df_out = df_in[df_in['FILTER'].str.contains('PASS|LowGQX')].reset_index(drop=True)

    return df_out

def filter_Q2(df_in: pd.DataFrame) -> pd.DataFrame:
    df_out = df_in[df_in['INFO'].str.contains("AF")].reset_index(drop=True)
    return df_out

def filter_Q3(df_in: pd.DataFrame, words: list) -> pd.DataFrame:
    
    df_out = df_in[df_in['FORMAT'].apply(lambda x: x.split(':')).
                 apply(lambda x: all([item in x for item in words]))]
    
    return df_out

def filter_Q4(df_in: pd.DataFrame) -> pd.DataFrame:
    df_nodot = df_in[df_in['ALT'] != "."] ## dot
    df_nodot_nona = df_nodot[df_nodot['ALT'].notna()] ## NA

    return df_nodot_nona

def normalize_vcf(vcf_in: str, vcf_out: str):

    print(f"Reading input vcf file: {vcf_in}")
    # read vcf file into pandas dataframe
    (vcf_head, df_vcf_body) = read_vcf(file_path=vcf_in)

    ## Q1: Remove any records that don't have a FILTER of "PASS" or "LowGQX".
    print('Retaining records with PASS or LowGQX...')
    df_pass_ = filter_Q1(df_vcf_body)

    ## Q2: The INFO column should only contain this variant's specific allele frequency, i.e. "AF={actual value}"
    print('Retaining records with AF=...')
    df_af_ = filter_Q2(df_pass_)

    ## Q3: The FORMAT column should only include the following values "GT:AD:DP"
    print(f'Retaining records with GT, AD, DP...')
    df_format_ = filter_Q3(df_in=df_af_, words=['GT', 'AD', 'DP'])

    ## Q4: The final file should only have variant rows with an alternate allele
    print(f'Final check...')
    df_norm_ = filter_Q4(df_format_)

    # write results to output file 
    print(f"Saving to normalized output vcf file: {vcf_out}")
    write_vcf(file_path=vcf_out, vcf_header=vcf_head, df_body=df_norm_)
    
    return "Hello World"
