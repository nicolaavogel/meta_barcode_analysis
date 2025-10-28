#!/bin/bash

module load bowtie2
module load samtools

FASTQ=$1
OUTPRE=$2
DB=$3

## Mapping with bowtie2

bowtie2 --threads 20 -k 1000 -t -x /projects/caeg/people/bfj994/barcodeMiner/DBs/viridiplantae_chloroplast/${DB} -U $FASTQ --no-unal -t | samtools view -bS - > $OUTPRE.bam

## bamfilter 


## metaDMG 

samtools sort -n -o $OUTPRE.sort.bam $OUTPRE.bam

/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp lca \
  --names /datasets/caeg_dataset/references/ncbi/20250530/taxonomy/ncbi/taxdump/names.dmp \
  --nodes /datasets/caeg_dataset/references/ncbi/20250530/taxonomy/ncbi/taxdump/nodes.dmp \
  --acc2tax /datasets/caeg_dataset/references/ncbi/20250530/taxonomy/ncbi/refseq_plastid.genomic.acc2taxid.gz \
  --sim_score_low 0.95 --sim_score_high 1.0 --how_many 15 --weight_type 0 \
  --fix_ncbi 0 --threads 10 \
  --bam $OUTPRE.sort.bam --out_prefix $OUTPRE
