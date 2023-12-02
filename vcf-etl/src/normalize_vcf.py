import logging
import re
import pandas as pd

# Note:  https://samtools.github.io/hts-specs/VCFv4.1.pdf

# READ VCF RETURN DATAFRAME AND METADATA
def read_vcf(vcf_in: str):

    # Extract metadata and colnames from vcf file
    with open(vcf_in, 'rt') as file:
        metadata = ''
        colnames = []

        for line in file:
            if line.startswith('##'):
                metadata = metadata + line

            elif line.startswith('#CHROM'):
                colnames = line.strip('#\n').split()
                break
        
        # Extract data as dataframe
        vcf = pd.read_csv(vcf_in, comment='#',header=None, names = colnames, delim_whitespace = True)
    
    return(vcf, metadata)

# SCREEN BY FILTER VALUES
def screen_filter(vcf: pd.DataFrame):
    vcf_new = vcf[vcf["FILTER"].isin(["PASS","LowGQX"])]
    return vcf_new 

# SCREEN BY ALTERNATE ALLELE
def screen_alt_allele(vcf: pd.DataFrame):
    vcf_new = vcf.drop(vcf[vcf.ALT == '.'].index)
    return vcf_new

# TRIM INFO COLUMN
def trim_info(vcf: pd.DataFrame):
    vcf['INFO'] = vcf['INFO'].map(get_AF)
    return vcf 

# Helper function to extract allele frequency
def get_AF(INFO: str):
    AF = ''
    af_pattern = r'AF1000G=([\d.]+)(;|$)'
    af_match = re.search(af_pattern, INFO)
    if af_match:
        AF = 'AF1000G=' + af_match.group(1)
    else:
        AF = '.'

    return AF


# TRIM FORMAT COLUMN
def trim_format(vcf: pd.DataFrame):
    return vcf.apply(trim_format_rowfunc, axis=1)

# Helper: apply per row
def trim_format_rowfunc(row: pd.Series):
    field_indexes = get_target_field_indexes(row['FORMAT'])
    row[8:] = row[8:].map(lambda x: extract_by_index(x, field_indexes))
    return row

def get_target_field_indexes(FORMAT: str):
    targets = ['GT','AD','DP']
    format_fields = FORMAT.split(':')
    indexes = list()
    
    # The target fields may not be all presents
    for t in targets:
        if t in format_fields:
            indexes.append(format_fields.index(t))

    return indexes

# Helper function to extract specific field data from entries
def extract_by_index(col: str, indexes: list):
    fields = col.split(':')
    return ":".join([fields[i] for i in indexes])


def write_vcf(vcf: pd.DataFrame, metadata: str, vcf_out: str):
    with open(vcf_out, 'wt') as file:
        file.write(metadata + '#')
    vcf.to_csv(vcf_out, sep = '\t', index = False, mode='a')

    logging.debug(f"Wrote '{vcf_out}'")
    return()

def get_outfile(vcf_in: str, vcf_out: str):
    if not vcf_out:
        if vcf_in.endswith('.vcf'):
            vcf_out = vcf_in[:-4] + '.normalized.vcf'
        else:
            vcf_out = vcf_in + '.normalized.vcf'

    return vcf_out

def normalize_vcf(vcf_in: str, vcf_out: str):

    vcf_out = get_outfile(vcf_in, vcf_out)

    vcf, metadata = read_vcf(vcf_in)
    
    vcf = screen_filter(vcf)
    
    vcf = screen_alt_allele(vcf)
    
    vcf = trim_info(vcf)
    
    vcf = trim_format(vcf)
    
    write_vcf(vcf, metadata, vcf_out)
    return()

# Run script solo with debug output
if __name__ == "__main__":
    import argparse
    import pprint
    pp = pprint.PrettyPrinter(indent=4)
    FORMAT = "%(filename)s:%(funcName)s:%(lineno)s - %(message)s"
    logging.basicConfig(format=FORMAT, level=logging.DEBUG)

    parser = argparse.ArgumentParser(description ='Normalize a VCF file.')

    parser.add_argument('-i','--infile',  type = str, required=True, help='Path to VCF file to be normalized.')
    parser.add_argument('-o', '--outfile', type = str, help = 'Path to where normalized VCF will be saved.')
    args, positional = parser.parse_known_args()

    normalize_vcf(args.infile, args.outfile)
