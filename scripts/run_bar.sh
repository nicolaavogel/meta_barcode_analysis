#!/bin/bash

# Set file size identifiers and their human-readable suffixes
sizes=("10000" "100000" "1000000" "50000" "500000" "250000")
labels=("10k" "100k" "1m" "50k" "500k" "250k")

# Loop over 10 iterations
for i in {1..10}; do
    echo "=== Iteration $i ==="
    for idx in "${!sizes[@]}"; do
        size="${sizes[$idx]}"
        label="${labels[$idx]}"
        fq="BN_$i_${size}.fq.gz"

        # Barcode pipeline
        ./bar_miner_pip2.sh "$fq" "../results/BN_${label}_barcode_$i" barcode/viri_barcode.fa

        # trnL pipeline
        ./bar_miner_pip2.sh "$fq" "../results/BN_${label}_trnl_$i" trnl/viri_trnl.fa

        # Full chloroplast pipeline
        ./bar_miner_pip2.sh "$fq" "../results/BN_${label}_full_$i" fasta/viri_full_chlo.fa
    done
done
