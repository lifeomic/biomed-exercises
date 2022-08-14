import os
from src.normalize_vcf import normalize_vcf
import src.normalize_vcf as normalize_module

def test_get_allele_frequency():
    """
    Tests the get_allele_frequency function

    Returns
    -------
    None.

    """
    
    # define a list of potential inputs
    inputs = ["END=10514;BLOCKAVG_min30p3a",
              "SNVHPOL=4;MQ=60;GMAF=G|0.1154;AF1000G=0.884585;phyloP=0.241",
              "SNVHPOL=2;MQ=48;AF=A|0.4778;AF1000G=0.522165;phyloP=-1.14;",
              "SNVHPOL=15;MQ=60;AA=T;GMF=G|0.07285;AF1000G=0.072848;phyloP=0.",
              "SNVHPOL=4;MQ=60;AA=G;GMAF=T|0.143;AF1000G=0.143046;"
        ]
    
    # define the required outputs
    outputs = ["",
               "GMAF=G|0.1154",
               "AF=A|0.4778",
               "",
               "GMAF=T|0.143"
        ]
    
    # generate predicted outputs
    predictions = [normalize_module.get_allele_frequency(input_value)
                   for input_value in inputs
        ]
    
    # compare results
    assert predictions[0] == outputs[0]
    assert predictions[1] == outputs[1]
    assert predictions[2] == outputs[2]
    assert predictions[3] == outputs[3]
    assert predictions[4] == outputs[4]
        
    return

def test_clear_format():
    """
    Test the clear_format function

    Returns
    -------
    None.

    """
    
    # define a list of potential inputs
    inputs = ["GT:AD:GQX:DP:DPF:AD:AD:ADR:SB:FT:PL",
              "GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL",
              "GT:GQX:DP:DPF",
              "GT:GQX:DP:DPF:GT:GQX:DP:DPF",
              "0/0:42:15:12"
              ]
    
    # define the correct outputs
    outputs = ["GT:AD:DP:AD:AD",
               "GT:DP:AD",
               "GT:DP",
               "GT:DP:GT:DP",
               ""
               ]
    
    # generate predicted outputs
    predictions = [normalize_module.clear_format(input_value)
                   for input_value in inputs
        ]
    
    # compare results
    assert predictions[0] == outputs[0]
    assert predictions[1] == outputs[1]
    assert predictions[2] == outputs[2]
    assert predictions[3] == outputs[3]
    assert predictions[4] == outputs[4]
        
    return

def test_normalize_vcf():
    """
    Tests the normalize_vcf function. Makes sure that it runs without throwing
    an error, and that it produces the target file.

    Returns
    -------
    None.

    """
    file_in = os.path.join('data', 'NA12878_cod_ex.genome.vcf')
    file_out = os.path.join('data', 'NA12878_cod_ex.genome_norm.vcf')
    actual = normalize_vcf(file_in, file_out)
    assert actual == "Hello World"
    assert os.path.isfile(file_out)
    
    return
