#!/bin/bash
#SBATCH --job-name=bowtie2-index
#SBATCH --output=./out/bwt.%j.out
#SBATCH --error=./out/bwt%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00
#SBATCH --partition=XXX
#SBATCH --output=./bowtie2build_%j.out
#SBATCH --error=./bowtie2build_%j.err

module load bowtie2

# Path to the genome FASTA file
GENOME_FA="XXX/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
# Base name for the index files
INDEX_BASE="XXX/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/Bowtie2_Index"

bowtie2-build $GENOME_FA $INDEX_BASE #--threads=16