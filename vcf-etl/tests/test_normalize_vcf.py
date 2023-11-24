from src.normalize_vcf import normalize_vcf, write_vcf, read_vcf, filter_filter, filter_info, filter_format, filter_allele

import pytest
import os

@pytest.fixture
def testVCF():
    return read_vcf("../data/NA12878_cod_ex.genome.vcf")[0]

@pytest.fixture
def testMetaData():
    return read_vcf("../data/NA12878_cod_ex.genome.vcf")[1]

def test_write_vcf(testVCF,testMetaData):
    outputVCF = 'unitTest.vcf'
    write_vcf(testVCF, testMetaData, outputVCF)
    # Check if file was created
    assert os.path.exists("unitTest.vcf")
    
def test_read_vcf(testVCF):
    actual = testVCF
    # Check shape to determine if read worked
    assert actual.shape == (37,10)
    
def test_filter_filter(testVCF):
    actual = filter_filter(testVCF)
    # Check shape and FILTER contents
    assert actual.shape == (30,10)
    assert set(actual['FILTER']) == {'PASS','LowGQX'}
    
def test_filter_info(testVCF):
    actual = filter_info(testVCF)
    assert actual.shape == (37,10)
    # Check for AF in INFO
    af_Vals = [x for x in actual['INFO'] if "AF" in x]
    assert len(af_Vals) > 0
    # Check for multiple inputs in INFO
    other_Vals = [x for x in actual['INFO'] if len(x.split(';')) > 1]
    assert len(other_Vals) == 0

def test_filter_format(testVCF):
    actual = filter_format(testVCF)
    assert actual.shape == (37,10)
    # Check for inputs outside of the ref_Vals
    ref_Vals = {'GT','AD','DP'}
    format_Vals = [x for x in actual['FORMAT'] if len(set(x.split(':')).union(ref_Vals)) > 3]
    assert len(format_Vals) == 0

def test_filter_allele(testVCF):
    actual = filter_allele(testVCF)
    assert actual.shape == (18,10)
    # Check for elements that contain a "." there should be none
    alt_Vals = [x for x in actual['ALT'] if x == '.']
    assert len(alt_Vals) == 0

