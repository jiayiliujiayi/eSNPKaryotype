import os
import argparse

def generate_slurm_script(input_dir, output_dir, output_filename, read_group_prefix):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Find all .fq.gz files in the input directory
    fq_files = sorted([f for f in os.listdir(input_dir) if f.endswith('.fq.gz')])

    # Assuming paired-end reads
    pairs = {}
    for file in fq_files:
        prefix = '_'.join(file.split('_')[:-1])  # Modify this if file naming convention changes
        if prefix not in pairs:
            pairs[prefix] = []
        pairs[prefix].append(file)

    # Check for any unpaired files
    for key, files in pairs.items():
        if len(files) != 2:
            raise ValueError(f"Unpaired files found for prefix {key}: {files}")

    # Generate commands for each pair
    commands = []
    for index, (prefix, files) in enumerate(pairs.items()):
        bam_path = os.path.join(output_dir, f"{prefix}_output.bam")
        command = f"""
echo "Processing {prefix}"

bowtie2 \\
    -x XXX/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/Bowtie2_Index \\
    -1 {os.path.join(input_dir, files[0])} \\
    -2 {os.path.join(input_dir, files[1])} \\
    --rg-id {read_group_prefix}{prefix} \\
    --rg LB:lib{prefix} --rg PL:ILLUMINA --rg SM:{prefix} --rg PU:{read_group_prefix}{prefix} \\
    -p 8 \\
    -S {bam_path.replace('.bam', '.sam')} && \\
    samtools view -bS {bam_path.replace('.bam', '.sam')} > {bam_path} 
"""
        commands.append(command)

    # Create the SLURM script content
    slurm_content = f"""#!/bin/bash
#SBATCH --job-name=bowtie2_alignment
#SBATCH --partition=XXX
#SBATCH --output={output_dir}/bowtie2_alignment_%j.out
#SBATCH --error={output_dir}/bowtie2_alignment_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=1-12:00:00

module load samtools
module load bowtie2

# Execute each command sequentially
{''.join(commands)}
    """

    # Write the script to a file
    with open(output_filename, 'w') as file:
        file.write(slurm_content)

    print(f"SLURM script '{output_filename}' has been created.")

def main():
    parser = argparse.ArgumentParser(description='Generate a SLURM script for fq to bam conversion using Bowtie2.')
    parser.add_argument('--input_dir', type=str, required=True, help='Directory containing input fq.gz files')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to output the BAM files')
    parser.add_argument('--output_filename', type=str, default='01.bowtie2.slurm', help='Filename for the generated SLURM script')
    parser.add_argument('--read_group_prefix', type=str, default='group', help='Prefix for read group IDs and PUs')

    args = parser.parse_args()

    generate_slurm_script(args.input_dir, args.output_dir, args.output_filename, args.read_group_prefix)

if __name__ == '__main__':
    main()
