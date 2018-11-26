#!/home/users/cjyoon/anaconda3/bin/python
import sys
import os
import subprocess
import shlex
import re
import argparse
import cyvcf2

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    parser.add_argument('-t', '--vaf_threshold', required=True, type=float,help='Remove variants where VAF in the panel of normal exceeds the threshold')

    args = vars(parser.parse_args())

    return args['input'], args['vaf_threshold'], args['output_dir']

def main():
    input, vaf_threshold, output_dir = argument_parser()
    outputfile = os.path.join(output_dir, re.sub('.vcf$', '.filtered.vcf', os.path.basename(input)))

    vcf_handle = cyvcf2.VCF(input)
    print(vcf_handle)
    writer = cyvcf2.Writer(outputfile, vcf_handle)

    for variant in cyvcf2.VCF(input):
        if variant.INFO['PON_VAF'] < vaf_threshold:
            writer.write_record(variant)

    vcf_handle.close()
    writer.close()
    
if __name__=='__main__':
    main()


