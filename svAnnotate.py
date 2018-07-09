import gzip 
import re
import pysam 
import os, sys
import random
import subprocess, shlex
import argparse
# to allow vep to annotate ALT column cannot be in SV annotation format in vcf. Change it to SNV annotation so that I can still utilize VEP's annotation


class Position():
    ''' python class for handling genomic positions
    0-based
    '''
    def __init__(self, chromosome, start, end, is_bp=None, clipped_reads=None):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.is_bp = is_bp
        self.clipped_reads = None

    def __repr__(self):
        return str(self.chromosome + ":" + str(self.start) + '-' + str(self.end))

    def __str__(self):
        return str(self.chromosome + ":" + str(self.start) + '-' + str(self.end))

    def __len__(self):
        return int(self.end - self.start)

    def __iter__(self):
        for i in range(self.start, self.end):
            yield Position(self.chromosome, i, i+1)

    def __next__(self):
        self.start += 1
        self.end += 1
        return self

    def __hash__(self):
        return hash((self.chromosome, self.start, self.end))

    def __eq__(self, other):
        if isinstance(other, Position):
            return self.chromosome == other.chromosome and self.start == other.start and self.end == other.end
        else:
            print("Not of the same class, cannot compare equality")
            return None

    @classmethod
    def fromstring(cls, position_string):
        if isinstance(position_string, str):
            chromosome = position_string.split(':')[0]
            start = int(position_string.split(':')[1].split('-')[0])
            end = int(position_string.split(':')[1].split('-')[1])
            return Position(chromosome, start, end)
        elif isinstance(position_string, Position):
            return position_string
        else:
            print('position_string has to be either a string class or Position class')
            raise TypeError
            
    @staticmethod
    def check_format(genomicPositionString):
        if re.search(r'[A-Za-z1-9]+:[0-9]+-[0-9]+', genomicPositionString):
            return True
        else:
            return False

    @staticmethod
    def overlap(position1, position2):
        '''true if position1 and position2 has more than 1 overlapping base'''
        try:
            if isinstance(position1, Position) and isinstance(position2, Position):
                if position1.chromosome == position2.chromosome:
                    if min(position1.end, position2.end) > max (position1.start, position2.start):
                        return True
                    else:
                        return False
                else:
                    return False  # cannot compare if two positions are in different chromosome
            else:
                return None # has to be Posiiton class.
        except:
            Exception

    def extend(self, direction, basepairs):
        """extends objects in by specified base pairs, either upstream, downstream, or both"""
        if direction=="up":
            return Position(self.chromosome, max(0, self.start-basepairs), end)
        elif direction=="down":
            return Position(self.chromosome, self.start, self.end + basepairs)
        elif direction=="both":
            return Position(self.chromosome, max(0, self.start - basepairs), self.end + basepairs)
        else:
            print('direction has to be either up, down, or both')
            raise ValueError

def get_canonical_annotation(long_annotation_string):
    '''vep output gives multiple annotations for a given loci, one is often only interested in the canonical sites
    thus looks for |YES| string to find canonical annotation
    if canonical annotation is not found, just return the very first annotation'''
    return_status = 0
    for anno in long_annotation_string.split(','):
        if re.search(r'|YES|', anno):
            return_status =1 
            return anno

    if return_status == 0:
        return annotation.split(',')[0]




def get_vep_annotation(genomic_position, reference_fasta):
    '''temporarily create a minimal SNV vcf for a given genomic position
    and then use this vcf to get variant information by running VEP
    then parse the VEP output vcf to get the necessary variant annotation info
    '''
    
    temp_vcf = write_temp_vcf(genomic_position, reference_fasta)
    # now run VEP annotation command
    cmd = f'python /home/users/cjyoon/scripts/annotate_vcf/annotate_vcf.py -i {temp_vcf}'
    subprocess.call(shlex.split(cmd))
    
    vep_vcf = re.sub(r'.vcf$', '.vep.vcf.gz', temp_vcf)
    with gzip.open(vep_vcf, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                annotation = line.strip().split()[7].split('CSQ=')[1]
                # multiple annotations per variant, now need to find which one is the canonical
                canonical_annotation = get_canonical_annotation(annotation)
    cleanup([vep_vcf, vep_vcf + '.tbi', temp_vcf])
    return canonical_annotation

def cleanup(fileList):
    '''list of files to clean up after done getting the annotation'''
    for afile in fileList:
        subprocess.call(shlex.split('rm -rf ' + afile))

    return 0

def make_temp_variant(genomic_position, reference_fasta):
    '''make a temporary VCF line that will be used as an input for VEP command'''
    gPosition = Position.fromstring(genomic_position)
    reference = pysam.FastaFile(reference_fasta)
    refbase = reference.fetch(gPosition.chromosome, gPosition.start, gPosition.start + 1)
    bases = ['A', 'T', 'G', 'C']

    if refbase != 'N':
        bases.remove(refbase)
    else:
        pass

    altbase = random.choice(bases)
    variantString = f'{gPosition.chromosome}\t{gPosition.start}\t.\t{refbase}\t{altbase}\t.\t.\t.'
    return variantString

def write_temp_vcf(genomic_position, reference_fasta, temp_dir='.'):
    '''write a temporary vcf file that will be used as an input for VEP command'''
    temp_vcf = os.path.join(temp_dir, str(genomic_position) + '.vcf')
    with open(temp_vcf, 'w') as f:
        f.write('##fileformat=VCFv4.1\n##reference=file:///path/to/human_g1k_v37.fasta\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        variantString = make_temp_variant(genomic_position, reference_fasta)
        f.write(variantString)

    return temp_vcf


def sv_annotation(bp1_position, bp2_position, reference_fasta):
    try:
        bp1_canonical = get_vep_annotation(bp1_position, reference_fasta)
        bp2_canonical = get_vep_annotation(bp2_position, reference_fasta)

        bp1_intersect_gene = bp1_canonical.split('|')[3]
        bp1_nearest_gene =  bp1_canonical.split('|')[35]
        bp1_exon = bp1_canonical.split('|')[8]
        bp1_intron =  bp1_canonical.split('|')[9]
        bp1_consequence = annotation_consequence_adjust(bp1_canonical.split('|')[1])
        bp1_strand = bp1_canonical.split('|')[20]

        bp2_intersect_gene = bp2_canonical.split('|')[3]
        bp2_nearest_gene =  bp2_canonical.split('|')[35]
        bp2_exon = bp2_canonical.split('|')[8]
        bp2_intron =  bp2_canonical.split('|')[9]
        bp2_consequence = annotation_consequence_adjust(bp2_canonical.split('|')[1])
        bp2_strand = bp2_canonical.split('|')[20]
        return ( '|'.join([bp1_intersect_gene, bp1_strand, bp1_exon, bp1_intron, bp1_consequence, bp1_nearest_gene]), '|'.join([bp2_intersect_gene, bp2_strand, bp2_exon, bp2_intron, bp2_consequence, bp2_nearest_gene]))

    except IndexError:
        print(bp1_canonical)
        print(bp1_canonical.split('|'))
        print(len(bp1_canonical.split('|')))
        print(bp2_canonical)
        print(bp2_canonical.split('|'))
        print(len(bp2_canonical.split('|')))
        print('exiting...')
        sys.exit()



def annotation_consequence_adjust(consequence):
    '''since we introduced an arbitrary point mutation to get the SV annotation from VEP, 
    consequences such as missense, stop_lost, stop_gain, are arbitrary as well'''
    if re.search(r'intergenic_variant', consequence):
        return 'intergenic_variant'
    elif re.search(r'intron_variant', consequence):
        return 'intron_variant'
    elif re.search(r'upstream_gene_variant', consequence):
        return 'upstream_gene_variant'
    elif re.search(r'downstream_gene_variant', consequence):
        return 'downstream_gene_variant'
    elif re.search(r'missense_variant|synonymous_variant|stop', consequence):
        return 'exonic_variant'
    else:
        return consequence

def argument_parser():
    parser = argparse.ArgumentParser(description='finds the gene overlapping/closest to the breakpoints involved in structural variations')
    parser.add_argument('breakpoints', nargs=2, help='two breakpoints separated by a space')
    args = vars(parser.parse_args())
    return args['breakpoints'][0], args['breakpoints'][1]

def main():
    bp1, bp2 = argument_parser()
    if Position.check_format(bp1) and Position.check_format(bp2):
        svanno = sv_annotation(bp1, bp2, '/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta')
        print(svanno)
    else:
        sys.exit()


if __name__=='__main__':
    main()


