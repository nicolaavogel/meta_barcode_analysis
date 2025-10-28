#!/usr/bin/env bash
set -euo pipefail

# usage: ./subsample_fastqs.sh fastq_list.txt output_dir
fastq_list=$1
outdir=$2

mkdir -p "$outdir"

# Subsample sizes
sizes=(10000 50000 100000 250000 500000 1000000)

# Copy and subsample
new_fastq_list="${outdir}/fastq_list.txt"
> "$new_fastq_list"

while read -r fq; do
    # Copy original file
    base=$(basename "$fq")
    cp "$fq" "$outdir/$base"
    echo "$outdir/$base" >> "$new_fastq_list"

    for n in "${sizes[@]}"; do
        out="$outdir/${base%.fastq}_sub${n}.fastq"
        seqtk sample -s100 "$fq" "$n" > "$out"
        echo "$out" >> "$new_fastq_list"
    done
done < "$fastq_list"

