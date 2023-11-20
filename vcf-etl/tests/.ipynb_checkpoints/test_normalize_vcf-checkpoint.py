from src.normalize_vcf import normalize_vcf, read_vcf, write_vcf, filter_Q1, filter_Q2, filter_Q3, filter_Q4 

IN_FILE = "/full/path/to/biomed-exercises/vcf-etl/data/NA12878_cod_ex.genome.vcf"
(vcf_head, df_vcf_body) = read_vcf(file_path=IN_FILE)

def test_read_vcf():

    # verify all lines of header start with ##
    assert all([line[:2] == '##' for line in vcf_head.split('\n') if len(line) != 0]), 'VCF Header is incorrectly formatted!'
    
    # verify that dataframe consists of required columns for further processing
    assert all([c in df_vcf_body.columns for c in ['ALT', 'FILTER', 'INFO', 'FORMAT']]), 'VCF Body does not have requisite columns!'

    return

def test_filters():

    # verify that FILTER column has PASS or LowGQX in each row (Q1)
    df_pass_ = filter_Q1(df_vcf_body)
    assert all(['PASS' in c or 'LowGQX' in c for c in df_pass_['FILTER'].to_list()]), 'FILTER column does not meet requirements!'

    df_af_ = filter_Q2(df_pass_)
    # verify that INFO column has AF in each row (Q2)
    assert all(['AF' in c for c in df_af_['INFO'].to_list()]), 'INFO column does not meet requirements!'

    df_format_ = filter_Q3(df_in=df_af_, words=['GT', 'AD', 'DP'])
    # verify that FORMAT column has either 'GT', 'AD' or 'DP' in each row (Q3)
    assert all(['GT' in c or 'AD' in c or 'DP' in c for c in df_format_['FORMAT'].to_list()]), 'FORMAT column does not meet requirements!'

    df_norm_ = filter_Q4(df_format_)
    # verify that ALT column has no NAs or .
    assert all([c != '' and c != '.' for c in df_norm_['ALT'].to_list()]), 'ALT column does not meet requirements!'
    
    return

def test_normalize_vcf():
    actual = normalize_vcf("/full/path/to/biomed-exercises/vcf-etl/data/NA12878_cod_ex.genome.vcf", 
                           "/full/path/to/biomed-exercises/vcf-etl/data/NA12878_normalized.vcf")

    assert actual == "Hello World"
