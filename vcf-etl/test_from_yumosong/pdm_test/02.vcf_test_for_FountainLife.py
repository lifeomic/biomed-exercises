import vcf
import argparse

parser = argparse.ArgumentParser(description='sequencing saturation')
parser.add_argument('-i','--infile', type=str, help='''path for input file''')
parser.add_argument('-o','--outfile',type=str, help='''path for output file''')
args = parser.parse_args()

def filter_vcf(input_file, output_file):

    vcf_reader = vcf.Reader(open(input_file, 'r'))
    vcf_writer = vcf.Writer(open(output_file, 'w'), vcf_reader)

    for record in vcf_reader:
        # rule 1: keep only variant rows with an alternate allele
        if record.num_hom_ref == len(record.samples):
            continue
        if all(sample['GT'] != '.' for sample in record.samples):

            # rule 2: remove records that don't have a FILTER of "PASS" or "LowGQX"
            if record.FILTER and "PASS" not in record.FILTER and "LowGQX" not in record.FILTER:
                continue

            # rule 3: the INFO column only contain specific allele frequency
            af_1000g_value = record.INFO.get('AF1000G', '.')
            gmaf_value = record.INFO.get('GMAF', '.')
            record.INFO = {'AF1000G': af_1000g_value, 'GMAF': gmaf_value}

            # rule 4: the first sample output information according to the format "GT:AD:DP"
            if record.samples:
                sample_data = record.samples[0]
                # extract values for "GT", "AD", and "DP" 
                gt_value = sample_data['GT'] if 'GT' in sample_data.data._fields else '.'
                ad_value = sample_data['AD'] if 'AD' in sample_data.data._fields else '.'
                dp_value = sample_data['DP'] if 'DP' in sample_data.data._fields else '.'

                # create a new Call object with the modified values
                new_call_data = vcf.model.make_calldata_tuple(['GT', 'AD', 'DP'])
                new_call = vcf.model._Call(record, sample_data.sample, new_call_data(*[gt_value, ad_value, dp_value]))

                # replace the Call object in record.samples
                record.samples[0] = new_call

            # rule 5: the FORMAT column only include "GT", "AD", and "DP" 
            format_columns_to_keep = ["GT", "AD", "DP"]
            record.FORMAT = ":".join(format_columns_to_keep)

            vcf_writer.write_record(record)

    vcf_writer.close()


if __name__ == "__main__":
    input_file_path = args.infile
    output_file_path = args.outfile
    filter_vcf(input_file_path, output_file_path)
