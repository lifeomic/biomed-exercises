import logging
import pandas as pd
import argparse

# Specify input file when using argparse
#parser = argparse.ArgumentParser()
#parser.add_argument('-i',"--inputFile",help="File name of input VCF", type=str, default="../data/NA12878_cod_ex.genome.vcf")
#parser.add_argument('-o',"--outputFile",help="File name of normalized VCF", type=str, default="../data/normalized.vcf")
#args = parser.parse_args()
#vcf_in = args.inputFile
#vcf_out = args.outputFile

def read_vcf(vcf_in: str):
    with open(vcf_in,'rt') as inFile:
        # Variable to hold metadata lines for vcf writing
        metadataLines = []
        for line in inFile:
            # Pull metadata lines into list
            if line.startswith("##"):
                metadataLines.append(line)
            # Pull header names for dataframe
            elif line.startswith("#CHROM"):
                headNames = line.strip("\n#").split("\t")
                break
    # Construct dataframe using header and ignore # lines
    vcf = pd.read_csv(vcf_in, comment='#', delim_whitespace=True, header=None, names=headNames)
    return(vcf, metadataLines)

def write_vcf(vcf, metadata, vcf_out: str):
    # Open write file and output metadata lines
    with open(vcf_out,'wt') as outFile:
        outFile.write("".join(metadata) + '#')
    # Append write file with new vcf dataframe
    vcf.to_csv(vcf_out, sep="\t", mode = 'a', index=False)
    return()

def filter_filter(vcf):
    # Define dataframe to select PASS or LowGQX
    vcf_filter = vcf[vcf["FILTER"].isin(["PASS","LowGQX"])]
    return (vcf_filter)

def filter_info(vcf):
    # Split the elements of INFO into iterable elements
    vcf["INFO"] = vcf["INFO"].transform(lambda x: x.split(';'))
    # Apply list comprehension to filter for AF element
    vcf["INFO"] = vcf["INFO"].transform(lambda x:[y for y in x if 'AF1000G' in y])
    # Pull out of list
    vcf["INFO"] = vcf["INFO"].transform(lambda x: ';'.join(x))
    return (vcf)

def filter_format(vcf):
    # Split the elements of FORMAT into iterable elements
    vcf["FORMAT"] = vcf["FORMAT"].transform(lambda x: x.split(':'))
    # Apply intersection to filter for GT:AD:DP elements
    vcf["FORMAT"] = vcf["FORMAT"].transform(lambda x: set(x).intersection({'GT','AD','DP'}))
    # Pull out of set
    vcf["FORMAT"] = vcf["FORMAT"].transform(lambda x: ':'.join(x))
    return (vcf)

def filter_allele(vcf):
    # Define dataframe to select for ALT that is not a '.'
    vcf_allele = vcf[vcf['ALT'] != '.']
    return (vcf_allele)

def normalize_vcf(vcf_in: str, vcf_out: str):
    #Read VCF
    vcf, metadata = read_vcf(vcf_in)
    
    #Remove any records that don't have a FILTER of "PASS" or "LowGQX"
    vcf = filter_filter(vcf)
    
    #The INFO column should only contain this variant's specific allele frequency, i.e. "AF={actual value}"
    vcf = filter_info(vcf)
    
    #The FORMAT column should only include the following values "GT:AD:DP"
    vcf = filter_format(vcf)
    
    #The final file should only have variant rows with an alternate allele.
    vcf = filter_allele(vcf)
    
    #Write out Normalized VCF
    write_vcf(vcf, metadata, vcf_out)

#testrun with argparse
#normalize_vcf(vcf_in, vcf_out)