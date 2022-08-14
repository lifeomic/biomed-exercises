import os
import logging
import pandas as pd

def read_input_file(vcf_in: str):
    """
    Read in a vcf file, and separate it out to the header, and the content

    Parameters
    ----------
    vcf_in : str
        Path to the input file to be read.

    Returns
    -------
    header: list
        The header rows, stored as a list
    df: dataframe
        The data, stored as a pandas dataframe

    """
    
    # check that the file exists
    text = f'The file {vcf_in} does not exist.'\
            ' Please specify an existing file'
    assert os.path.isfile(vcf_in), text
    
    # read the file via list comprehension
    with open(vcf_in) as read_file:
        data = [line.strip('\n') for line in read_file]
    
    # separate into header and non-header
    header = [line for line in data if line[0] == '#']
    data = [line.split('\t') for line in data if line[0] != '#']
    
    # convert data into a pandas dataframe for easier manipulation
    # we need column labels, and need them all capitalized for consistency
    columns = [label.upper() for label in header[-1][1:].split('\t')]
    data = [data_row[:len(columns)] for data_row in data]
    df = pd.DataFrame(data, columns=columns)
    
    return header, df

def get_allele_frequency(information: str):
    """
    Take the information data on a variant row, and only return allele
    frequency

    Parameters
    ----------
    information : str
        The full information data on a variant row

    Returns
    -------
    information : str
        The information data on a variant row, stripped to allele frequency

    """
    
    # the following can be merged into a single line command, but should be
    # kept separate for improved readability
    
    # split information into groups
    information_list = information.split(';')
    
    # use list comprehension to keep targets of interest
    information_list = [info for info in information_list if "AF=" in info]
    
    # merge and return
    information = ';'.join(information_list)
    
    return information

def clear_format(format_data: str):
    """
    Take the format data from a variant row, and remove invalid formats

    Parameters
    ----------
    format_data : str
        The format data for a variant row

    Returns
    -------
    format_data : str
        The fixed format data for a variant row

    """
    
    # the following can be merged into a single line command, but should be
    # kept separate for improved readability
    
    # split format inot groups
    format_list = format_data.split(':')
    
    # use list comprehension to keep valid formats only
    format_list = [form for form in format_list if form in ('GT', 'AD', 'DP')]
    
    # merge and return
    format_data = ':'.join(format_list)
    
    return format_data

def fix_vcf_file(df):
    """
    Get the data from a read in vcf file, stored as a pandas dataframe.
    Convert it into a standardized format
    
    Parameters
    ----------
    df : pandas dataframe
        A pandas dataframe of data to be reformatted

    Returns
    -------
    None.

    """
    
    # we remove all rows with a wrong FILTER value
    to_keep = df.FILTER.isin(['PASS', "LowGQX"])
    to_drop = df.index[~to_keep]
    df.drop(to_drop, axis=0, inplace=True)
    
    # remove all but allele frequency from INFO column
    df.INFO = [get_allele_frequency(info) for info in df.INFO]
    
    # Only keep GT, AD, and DP for formats
    df.FORMAT = [clear_format(format_data) for format_data in df.FORMAT]
    
    # Drop rows without an alternate allele
    to_drop = df.index[df.ALT == '.']
    df.drop(to_drop, axis=0, inplace=True)
    to_drop = df.index[df.ALT == df.REF]
    df.drop(to_drop, axis=0, inplace=True)
        
    return

def save_output_file(vcf_out: str, header, df):
    """
    Compile a header and variant row data into a vcf file and save it

    Parameters
    ----------
    vcf_out : str
        The path to save the output file to
    header : list
        A list of the header data of the vcf file
    df : Pandas dataframe
        a pandas dataframe of the normalized vcf data

    Returns
    -------
    None.

    """
    
    # make sure that the directory we are saving to exists
    directory = os.path.dirname(vcf_out)
    text = (f'The directory {directory} to save {vcf_out} into does not'\
            'exist. Please create the directory.'
            )
    
    if len(directory) > 0:
        assert os.path.isdir(directory), text
    
    # save first header, then the data
    with open(vcf_out, 'w') as write_file:
        write_file.write('\n'.join(header))
        write_file.write('\n')
        
        df.to_csv(write_file, sep='\t', index=False, header=False)
    
    return

def normalize_vcf(vcf_in: str, vcf_out: str):
    """
    Read in a vcf file, and normalize it to a standardized format
    
    Parameters
    ----------
    vcf_in : str
        Path to the input file to be read.
    vcf_out : str
        Path to the output file to be read

    Returns
    -------
    output : str
        The string "Hello World"
    """
    
    # read in the file
    header, df = read_input_file(vcf_in)
    
    # fix the file
    fix_vcf_file(df)
    
    # save the results
    save_output_file(vcf_out, header, df)
    
    return "Hello World"
