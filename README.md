# Annotate VCF with VEP
```
$ python annotate_vcf.py -h
usage: annotate_vcf.py [-h] -i INPUT_VCF

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_VCF, --input_vcf INPUT_VCF
                        Input vcf or vcf.gz file to be annotated with VEP
```

# Annotate VCF with Panel of Normal VAFs

```
$ python pon_annotate.py -h
usage: pon_annotate.py [-h] -i INPUT [-r REFERENCE [REFERENCE ...]]
                       [-n NORMAL_BAMS [NORMAL_BAMS ...]] [-o OUTPUT_DIR]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input vcf to be annotated with panel of normal
  -r REFERENCE [REFERENCE ...], --reference REFERENCE [REFERENCE ...]
                        Reference fasta
  -n NORMAL_BAMS [NORMAL_BAMS ...], --normal_bams NORMAL_BAMS [NORMAL_BAMS ...]
                        Normal Bams to calculate panel of normal VAFS
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory
```

