#!/bin/bash

# Define your base commands (species, prefix, number of reads)
commands=(
    "Betula_nana BN 1000000"
    "Betula_nana BN 100000"
    "Betula_nana BN 10000"
    "Betula_nana BN 50000"
    "Betula_nana BN 500000"
    "Betula_nana BN 250000"
)

# Run each command 10 times
for i in {1..10}; do
    echo "=== Iteration $i ==="
    for cmd in "${commands[@]}"; do
        # Split the command into parts
        species=$(echo $cmd | awk '{print $1}')
        base_prefix=$(echo $cmd | awk '{print $2}')
        num_reads=$(echo $cmd | awk '{print $3}')
        
        # Create a unique prefix
        prefix="${base_prefix}_$i"
        
        echo "Running: ./rungar.sh $species $prefix $num_reads"
        ./rungar.sh "$species" "$prefix" "$num_reads"
    done
done

