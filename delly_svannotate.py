
from svAnnotate import *
import cyvcf2
import re
import argparse
import os

def argument_parser():
    parser = argparse.ArgumentParser(description='Writes a TSV file with break point annotations and orientations')
    parser.add_argument('-i', '--input_vcf', help='Input delly vcf to be annotated. Must have CHR2= and END= in the INFO column')
    parser.add_argument('--header', help='Do not write header if 0', type=int, default=1)
    parser.add_argument('-o', '--output_dir',help='Output directory', default=os.getcwd())

    args = vars(parser.parse_args())
    return args['input_vcf'], bool(args['header']), args['output_dir']

def main():
    vcfpath, header, output_dir = argument_parser()
    outputbasename = os.path.basename(vcfpath) + '.svanno.txt' 
    outputfile = os.path.join(output_dir, outputbasename)
    
    with open(outputfile ,'w') as f:
        for variant in cyvcf2.VCF(vcfpath):
            bp1 = f'{variant.CHROM}:{variant.POS}-{variant.POS + 1}'
            bp2 = f"{variant.INFO.get('CHR2')}:{variant.INFO.get('END')}-{variant.INFO.get('END')+1}"
            orientation = variant.INFO.get('CT')
            svtype = variant.INFO.get('SVTYPE')
            bp1_annotation, bp2_annotation = sv_annotation(bp1, bp2, '/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta')
            f.write(f"{bp1}\t{bp2}\t{svtype}\t{orientation}\t{bp1_annotation}\t{bp2_annotation}\n")

if __name__=='__main__':
    main()

