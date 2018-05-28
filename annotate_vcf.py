import subprocess
import shlex
import argparse
import re
import os 
def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_vcf', required=True, help='Input vcf or vcf.gz file to be annotated with VEP')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    args = vars(parser.parse_args())
    
    return args['input_vcf'], args['output_dir']

def outputfile(vcffile, output_dir):
    if vcffile.endswith('vcf.gz'):
        outputfile = os.path.join(output_dir, os.path.basename(re.sub(string=vcffile, pattern=r'vcf.gz$', repl='vep.vcf')))
    elif vcffile.endswith('vcf'):
        outputfile = os.path.join(output_dir, os.path.basename(re.sub(string=vcffile, pattern=r'vcf$', repl='vep.vcf')))
    else:
        print('Input must be a file that has a suffix of vcf or vcf.gz')
        raise ValueError

    return outputfile 

def vep_annotate(input_vcf, annotated_vcf):
    ''' if you want single annotation then use --per_gene option'''

    cmd = f'/home/users/cjyoon/ensembl-vep/vep --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length  --offline --no_progress --no_stats  --polyphen b  --regulatory --af --af_1kg --af_gnomad --af_esp --max_af -i {input_vcf} -o {annotated_vcf} --force_overwrite --nearest symbol'
    print(cmd)
    vep_cmd = subprocess.Popen(shlex.split(cmd))
    vep_cmd.wait()
    return 0

def bgzip_tabix(annotated_vcf):
    bgzipcmd=subprocess.Popen(shlex.split(f'bgzip {annotated_vcf}'))
    bgzipcmd.wait()
    tabixcmd = subprocess.Popen(shlex.split(f'tabix -p vcf {annotated_vcf}.gz'))
    tabixcmd.wait()

    return 0

def main():
    input_vcf, output_dir = argument_parser()
    
    # prepare annotated output file path 
    annotated_vcf = outputfile(input_vcf, output_dir)
    # run VEP
    vep_annotate(input_vcf, annotated_vcf)
    # bgzip and tabix output vcf
    bgzip_tabix(annotated_vcf)
    

if __name__=='__main__':
    main() 

