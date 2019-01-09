#!/home/users/cjyoon/anaconda3/bin/python
# 2018.09.19 cjyoon edit with isolated conda environment
# 2018.12.15 cjyoon edit to output only the canonical annotation in the CSQ INFO field. 
# 2018.12.15 Also supports GRCh37/GRCh38 flexibility with -g argument


import subprocess
import shlex
import argparse
import re
import os 
import cyvcf2

# Local GRChX installation and Cache version for assemblies
cache_version_mapper = dict({
    "GRCh37": 91, 
    "GRCh38": 93
    })


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_vcf', required=True, help='Input vcf or vcf.gz file to be annotated with VEP')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    parser.add_argument('-g', '--genome_assembly', required=False, default='GRCh37', help='Genome Assembly version to use. Default=GRCh37')
    args = vars(parser.parse_args())
    
    return args['input_vcf'], args['output_dir'], args['genome_assembly']

def outputfile(vcffile, output_dir):
    if vcffile.endswith('vcf.gz'):
        outputfile = os.path.join(output_dir, os.path.basename(re.sub(string=vcffile, pattern=r'vcf.gz$', repl='vep.vcf')))
    elif vcffile.endswith('vcf'):
        outputfile = os.path.join(output_dir, os.path.basename(re.sub(string=vcffile, pattern=r'vcf$', repl='vep.vcf')))
    else:
        print('Input must be a file that has a suffix of vcf or vcf.gz')
        raise ValueError

    return outputfile 

def vep_annotate(input_vcf, temp_annotated_vcf, assembly_version, cache_version, display=True):
    ''' if you want single annotation then use --per_gene option'''
    # DEFINE PERL5LIB PATH for other people's use
    os.environ['PERL5LIB'] = '/home/users/cjyoon/anaconda3/envs/vep/lib/perl5/site_perl/5.22.0/x86_64-linux-thread-multi' 
    os.environ['PATH'] = '/home/usrs/cjyoon/anaconda3/bin:' + os.environ['PATH']
    os.environ['PATH'] = '/home/users/cjyoon/anaconda3/envs/vep/bin:' + os.environ['PATH']
    print(os.environ['PATH'])    
# DEFINE VEP CACHE directory to use
#    dir = '/home/users/cjyoon/vep_dir'
    dir = '/home/users/cjyoon/.vep'
    activate_environment = os.system('source activate /home/users/cjyoon/anaconda3/envs/vep')

#    VEP_PATH = '/home/users/cjyoon/anaconda3/envs/vep/bin/variant_effect_predictor.pl'
    # 2018.11.13 cjyoon
    VEP_PATH = '/home/users/cjyoon/anaconda3/envs/vep/bin/vep'
    CACHE_VER = f' --cache_version {cache_version}'
    ASSEMBLY_VER = assembly_version
    cmd = f'{VEP_PATH} --dir {dir} --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length  --offline --no_progress --no_stats  --polyphen b  --regulatory --af --af_1kg --af_gnomad --af_esp --max_af -i {input_vcf} -o {temp_annotated_vcf} --force_overwrite --nearest symbol {CACHE_VER} --assembly {ASSEMBLY_VER}'

    if display==True:
        print(cmd)

    vep_cmd = subprocess.Popen(shlex.split(cmd))
    vep_cmd.wait()
    return 0

def find_canonical_annotation(vep_annotation_string):
    """VEP annotates with many alternative transcripts as well as canonical transcript
    this function finds the canonical transcript within vep_annotation_string. 
    If there is no canonical transcript, which is usually the case fore intergenic, 
    will just report the first annotation. 
    """
    annotations = vep_annotation_string.split(',')
    return_status = 0
    for annotation in annotations:
        CANONICAL = annotation.split('|')[26] # CANONICAL
        if CANONICAL == 'YES':
            return_status = 1
            return annotation
    
    if return_status == 0:
        return vep_annotation_string.split(',')[0]

def get_canonical_annotation(temp_annotated_vcf, annotated_vcf):
    """get only the canonical annotation into the CSQ INFO field"""
    vcf_handle = cyvcf2.VCF(temp_annotated_vcf)
    vcf_writer = cyvcf2.Writer(annotated_vcf, vcf_handle)
    print('writing ' + annotated_vcf)
    for variant in vcf_handle:
        canonical_annotation = find_canonical_annotation(variant.INFO['CSQ'])
        variant.INFO['CSQ'] = canonical_annotation
        vcf_writer.write_record(variant)

    vcf_writer.close() 

    return 0

def bgzip_tabix(annotated_vcf):
    bgzipcmd=subprocess.Popen(shlex.split(f'bgzip -f {annotated_vcf}'))
    bgzipcmd.wait()
    tabixcmd = subprocess.Popen(shlex.split(f'tabix -p vcf {annotated_vcf}.gz'))
    tabixcmd.wait()

    return 0

def main():
    input_vcf, output_dir, genome_assembly = argument_parser()
    
    cache_version = cache_version_mapper[genome_assembly]

    # prepare annotated output file path 
    annotated_vcf = outputfile(input_vcf, output_dir)
    temp_annotated_vcf = annotated_vcf + '.temp.vcf'

    # run VEP
    vep_annotate(input_vcf, temp_annotated_vcf, genome_assembly, cache_version)

    # re-write with only canonical variants
    get_canonical_annotation(temp_annotated_vcf, annotated_vcf)
    
    # bgzip and tabix output vcf
    bgzip_tabix(annotated_vcf)
    

if __name__=='__main__':
    main() 

