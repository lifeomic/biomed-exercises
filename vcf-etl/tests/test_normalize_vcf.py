import os
from unittest.mock import patch
from src.normalize_vcf import normalize_vcf, replace_with_AF, replace_format
import pandas as pd

def test_replace_with_AF_exists():
    # Test when 'AF1000G' exists in the input string
    info = 'DP=100;AF1000G=0.05;GQ=30'
    result = replace_with_AF(info)
    assert result == 'AF=0.05'

def test_replace_with_AF_not_exists():
    # Test when 'AF1000G' does not exist in the input string
    info = 'DP=100;GQ=30'
    result = replace_with_AF(info)
    assert result == None

def test_replace_with_AF_multiple_values():
    # Test when 'AF1000G' exists multiple times in the input string
    info = 'DP=100;AF1000G=0.05;GQ=30;AF1000G=0.1'
    result = replace_with_AF(info)
    assert result =='AF=0.05'  # Only the first occurrence should be replaced

def test_replace_with_AF_empty_input():
    # Test with an empty input string
    info = ''
    result = replace_with_AF(info)
    assert result == None

def test_replace_format_all_present():
    input_format = 'GT:GQ:AD:DP:PL'
    result = replace_format(input_format)
    assert result == 'GT:AD:DP'

def test_replace_format_some_missing():
    input_format = 'GT:GQ:PL'
    result = replace_format(input_format)
    assert result == 'GT'

def test_replace_format_none_present():
    input_format = 'GQ:PL'
    result = replace_format(input_format)
    assert result == ''

def test_replace_format_empty_input():
    input_format = ''
    result = replace_format(input_format)
    assert result == ''

#@patch('pandas.read_csv')
@patch('src.normalize_vcf.read_vcf')
@patch('src.normalize_vcf.replace_with_AF')
@patch('src.normalize_vcf.replace_format')
def test_normalize_vcf(mock_replace_format, mock_replace_with_AF, mock_readvcf):
    test_data = {
    'FILTER': ['PASS', 'LowGQX', 'other', 'PASS'],
    'ALT': ['A', 'T', 'G','.'],
    'INFO': ['AF1000G=0.05', 'AF1000G=0.1', 'DP=100;AF1000G=0.2','DP=100'],
    'FORMAT': ['GT:GQ:AD:DP:PL', 'GT:GQ:AD:', 'GT:GQ:AD:DP:PL','GT:GQ:AD:DP:PL']
    }
    mock_readvcf.return_value = pd.DataFrame(test_data)
    mock_replace_with_AF.return_value = "AF=0.05"
    mock_replace_format.return_value = "GT:AD:DP"
    normalize_vcf('input.vcf', 'output.vcf')
    print('completed')
    if os.path.exists('output.vcf'):
        print('found................................................')
    # Assertions
    #mock_read_csv.assert_called_once_with('input.vcf', sep='\t', comment='#')
    #mock_replace_with_AF.assert_called()
   # mock_replace_format.assert_called()
    #mock_to_csv.assert_called_once_with('output.vcf', sep='\t', index=False)
    #with patch('pandas.read_csv', wraps=pd.read_csv):
    output_df = pd.read_csv('output.vcf', sep='\t', comment='#')
    assert len(output_df) == 2
    assert (output_df['INFO'] == 'AF=0.05').all()
    assert (output_df['FORMAT'] == 'GT:AD:DP').all()
    assert (output_df['FILTER'].isin(['PASS','LowGQX'])).all()
    assert (output_df['ALT'] != '.').all()
    if os.path.exists('output.vcf'):
        print('found................................................')
        os.remove('output.vcf')


