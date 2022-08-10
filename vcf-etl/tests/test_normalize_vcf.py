from src.normalize_vcf import normalize_vcf


def test_normalize_vcf():
    actual = normalize_vcf("file-in", "file-out")
    assert actual == "Hello World"
