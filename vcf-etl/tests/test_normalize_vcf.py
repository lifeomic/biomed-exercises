from src.normalize_vcf import normalize_vcf

def test_normalize_vcf():
    actual = normalize_vcf("tests/testData/header.vcf", "file-out")
    predicted = "##fileformat=VCFv4.1"
    assert actual == predicted
    print("Header OK")
    
    actual = normalize_vcf("tests/testData/colHeaders.vcf", "file-out")
    predicted = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878"
    assert actual == predicted
    print("Column Header OK")
    
    actual = normalize_vcf("tests/testData/normalData.vcf", "file-out")
    predicted = "chr1\t917495\trs13303369\tC\tT\t655\tPASS\tAF=0.494209\tGT:DP:AD\t1/1:51:0,51"
    assert actual == predicted
    print("Normal Data OK")
    
    actual = normalize_vcf("tests/testData/normalDataLowQual.vcf", "file-out")
    predicted = ""
    assert actual == predicted
    print("Low Quality Filter OK")
    
    actual = normalize_vcf("tests/testData/noVariant.vcf", "file-out")
    predicted = ""
    assert actual == predicted
    print("Vairant Filter OK")
    
    actual = normalize_vcf("tests/testData/nullData.vcf", "file-out")
    predicted = "" #skipped since doesn't pass variant filter or low quality filter
    assert actual == predicted
    print("Null Data OK")
    
    actual = normalize_vcf("tests/testData/emptyLine.vcf", "file-out")
    predicted = ""
    assert actual == predicted
    print("No Data OK")

    print("All tests passed!")




test_normalize_vcf()


    # print(repr(actual))
    # print(repr(predicted))

# colHeaders.vcf        emptyLine.vcf       header.vcf      noVariant.vcf       normalData.vcf      normalDataLowQual.vcf   nullData.vcf