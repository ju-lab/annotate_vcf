#!/home/users/cjyoon/anaconda3/bin/python
'''Given a vcf and a list of normal bams, will return VAF of that position in panel of normals an will annotate the VCF with normal VAF
For each position will run mpileup at a single base. 

2018.07.09 CJY
'''
import sys
import os
import subprocess
import shlex
import re
import argparse
import cyvcf2
def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='Input vcf to be annotated with panel of normal')
    parser.add_argument('-r', '--reference',  help='Reference fasta', default='/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta', nargs='+')

    parser.add_argument('-n', '--normal_bams', help='Normal Bams to calculate panel of normal VAFS', nargs='+')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')

    args = vars(parser.parse_args())

    return args['input'], args['normal_bams'], args['output_dir'], args['reference']

def calculate_vaf(normal_bam_list, query_position, reference):
	'''calculates the VAF in the normal Bams in the list of a given query_position'''
	bamlistString = '\t'.join(normal_bam_list)
	cmd = f'samtools mpileup -f {reference} -r {query_position} {bamlistString}'
	mpileup = subprocess.check_output(shlex.split(cmd), stderr=subprocess.DEVNULL)
	split_mpileup = mpileup.decode("utf-8").split('\t')


	totalDepth = 0 
	mismatches = 0
	totalCharacters = 0
	for i in range(0, int(len(split_mpileup)/3) -1 ):
		bases = 3*i + 1
		depths = 3*i
		totalDepth += int(split_mpileup[depths])
		baseAStringCount = split_mpileup[bases].upper().count('A')
		baseTStringCount = split_mpileup[bases].upper().count('T')
		baseGStringCount = split_mpileup[bases].upper().count('G')
		baseCStringCount = split_mpileup[bases].upper().count('C')
		mismatchCount = baseAStringCount + baseTStringCount + baseGStringCount + baseCStringCount
		mismatches += mismatchCount

	return round(float(mismatches/totalDepth), 3)

def main():
    vcf, normal_bams, output_dir, reference = argument_parser()
    
    vcf_handle = cyvcf2.VCF(vcf)
    vcf_handle.add_info_to_header({'ID': 'PON_VAF', 'Description': 'VAF in Panel of Normals',
    'Type':'Float', 'Number': '1'})

    output_vcf = os.path.join(output_dir, re.sub(r'.vcf$', '.pon.vcf', os.path.basename(vcf)))

    output_handle = cyvcf2.Writer(output_vcf, vcf_handle)
    for variant in vcf_handle:
    	variant_position = f'{variant.CHROM}:{variant.POS}-{variant.POS}'
    	pon_vafs =  calculate_vaf(normal_bams, variant_position, reference)
    	variant.INFO['PON_VAF'] = str(pon_vafs)
    	output_handle.write_record(variant)


if __name__=='__main__':
    main()


