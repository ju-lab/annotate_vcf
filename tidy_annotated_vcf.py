"""This script reads in VEP annotated VCF and parses into tidy format for analysis
April 8 2018 Jongsoo Yoon (cjyoon@kaist.ac.kr)
"""

import cyvcf2
import re
import argparse
import os, sys
import subprocess
import shlex

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_vcf', required=True, help='Input VCF with VEP annotation that will be tidied')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')

    args = vars(parser.parse_args())
    return args['input_vcf'], args['output_dir']


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

def prepare_outputfile(input_vcf, output_dir):
    """Prepares outputfile path. If an outputfile already exists, overwrite by removing original file"""
    output_file = os.path.join(output_dir, os.path.basename(input_vcf) + '.tidy.txt')
    if os.path.isfile(output_file):
        subprocess.call(shlex.split('rm -rf ' + output_file))

    return output_file

def tidy_annotation(input_vcf, output_dir):
    vcfHandle = cyvcf2.VCF(input_vcf)
    outputfile = prepare_outputfile(input_vcf, output_dir)
    with open(outputfile, 'w') as f:
        # write header
        f.write('chromosome\tposition\tref\talt\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tALLELE_NUM\tDISTANCE\tSTRAND\tFLAGS\tPICK\tVARIANT_CLASS\tHGNC_ID\tCANONICAL\tTSL\tCCDS\tENSP\tSWISSPROT\tTREMBL\tUNIPARC\tRefSeq\tGENE_PHENO\tSIFT\tPolyPhen\tDOMAINS\tHGVS_OFFSET\tAF\tAFR_AF\tAMR_AF\tEAS_AF\tEUR_AF\tSAS_AF\tAA_AF\tEA_AF\tgnomAD_AF\tgnomAD_AFR_AF\tgnomAD_AMR_AF\tgnomAD_ASJ_AF\tgnomAD_EAS_AF\tgnomAD_FIN_AF\tgnomAD_NFE_AF\tgnomAD_OTH_AF\tgnomAD_SAS_AF\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tMOTIF_NAME\tMOTIF_POS\tHIGH_INF_POS\tMOTIF_SCORE_CHANGE')

        for variant in vcfHandle:
            is_canonical_present = False
            annotation = find_canonical_annotation(variant.INFO['CSQ'])
            Allele, Consequence, IMPACT, SYMBOL, Gene, Feature_type, Feature, BIOTYPE, \
            EXON, INTRON, HGVSc, HGVSp, cDNA_position, CDS_position, Protein_position, \
            Amino_acids, Codons, Existing_variation, ALLELE_NUM, DISTANCE, STRAND, FLAGS, \
            PICK, VARIANT_CLASS, SYMBOL_SOURCE, HGNC_ID, CANONICAL, TSL, CCDS, ENSP, \
            SWISSPROT, TREMBL, UNIPARC, RefSeq, GENE_PHENO, SIFT, PolyPhen, DOMAINS, \
            HGVS_OFFSET, AF, AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF, AA_AF, EA_AF, \
            gnomAD_AF, gnomAD_AFR_AF, gnomAD_AMR_AF, gnomAD_ASJ_AF, gnomAD_EAS_AF, \
            gnomAD_FIN_AF, gnomAD_NFE_AF, gnomAD_OTH_AF, gnomAD_SAS_AF, MAX_AF, \
            MAX_AF_POPS, CLIN_SIG, SOMATIC, PHENO, PUBMED, MOTIF_NAME, MOTIF_POS, \
            HIGH_INF_POS, MOTIF_SCORE_CHANGE = annotation.split('|')
            # write either the canonical variant, or if not available, first annotated info
            f.write(f'\n{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t{Consequence}\t{IMPACT}\t{SYMBOL}\t{Gene}\t{Feature_type}:{Feature}\t{BIOTYPE}\t{EXON}\t{INTRON}\t{HGVSc}\t{HGVSp}\t{cDNA_position}\t{CDS_position}\t{Protein_position}\t{Amino_acids}\t{Codons}\t{Existing_variation}\t{ALLELE_NUM}\t{DISTANCE}\t{STRAND}\t{FLAGS}\t{PICK}\t{VARIANT_CLASS}\t{SYMBOL_SOURCE}:{HGNC_ID}\t{CANONICAL}\t{TSL}\t{CCDS}\t{ENSP}\t{SWISSPROT}\t{TREMBL}\t{UNIPARC}\t{RefSeq}\t{GENE_PHENO}\t{SIFT}\t{PolyPhen}\t{DOMAINS}\t{HGVS_OFFSET}\t{AF}\t{AFR_AF}\t{AMR_AF}\t{EAS_AF}\t{EUR_AF}\t{SAS_AF}\t{AA_AF}\t{EA_AF}\t{gnomAD_AF}\t{gnomAD_AFR_AF}\t{gnomAD_AMR_AF}\t{gnomAD_ASJ_AF}\t{gnomAD_EAS_AF}\t{gnomAD_FIN_AF}\t{gnomAD_NFE_AF}\t{gnomAD_OTH_AF}\t{gnomAD_SAS_AF}\t{MAX_AF}\t{MAX_AF_POPS}\t{CLIN_SIG}\t{SOMATIC}\t{PHENO}\t{PUBMED}\t{MOTIF_NAME}\t{MOTIF_POS}\t{HIGH_INF_POS}\t{MOTIF_SCORE_CHANGE}')




    print(f'Tidy Variant Annotation File Ready for Analysis: {outputfile}')
    return 0
 

def main():
    input_vcf, output_dir = argument_parser()
    tidy_annotation(input_vcf, output_dir)

    return 0

if __name__ == '__main__':
    main()

