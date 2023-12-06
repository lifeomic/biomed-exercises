import pandas as pd
import io

def replace_with_AF(info):
    allInfo = info.split(';')

    af_value = next((item for item in allInfo if 'AF1000G' in item), None)
    if af_value is not None:
        return af_value.replace('AF1000G','AF')
    else:
        return af_value
   
def replace_format(format):
     formats = format.split(':')
     gt_ad_dp = [element for element in formats if element == 'GT' or element == 'AD' or element == 'DP']
     return ':'.join(filter(None, gt_ad_dp))

def normalize_vcf(vcf_in: str, vcf_out: str):

    vcf_df = read_vcf(vcf_in)

    filtered_df = vcf_df.query("(FILTER == 'PASS' or FILTER == 'LowGQX') and ALT != '.'")

    filtered_df['INFO'] = filtered_df['INFO'].apply(lambda x:replace_with_AF(x))
    filtered_df['FORMAT'] = filtered_df['FORMAT'].apply(lambda x:replace_format(x))

    filtered_df.to_csv(vcf_out, sep='\t', index = False)

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        sep='\t'
    )
