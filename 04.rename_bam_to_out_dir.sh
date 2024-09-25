#!/bin/bash

# Step 1: Create the directory 'out' in the current working directory if it doesn't already exist
mkdir -p ./out

# Step 2: Move and rename *.bam files
for file in ./temp/*_output.sorted.reheaded.bam; do
    # Construct the new file path in the 'out' directory with the new name
    new_file="./out/$(basename "$file" _output.sorted.reheaded.bam).bam"
    # Move and rename the file
    mv "$file" "$new_file"
done

# Step 3: Move and rename *.bai files
for file in ./temp/*_output.sorted.reheaded.bam.bai; do
    # Construct the new file path in the 'out' directory with the new name
    new_file="./out/$(basename "$file" _output.sorted.reheaded.bam.bai).bam.bai"
    # Move and rename the file
    mv "$file" "$new_file"
done
