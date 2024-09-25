import os
import argparse

def generate_slurm_script(input_dir, output_dir, output_filename):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Find all .bam files that follow the specific naming convention
    #bam_files = sorted([f for f in os.listdir(input_dir) if f.endswith('_output.sorted.reheaded.bam')])
    bam_files = sorted([f for f in os.listdir(input_dir) if f.endswith('.bam')])

    # Generate commands for each bam file
    commands = []
    for bam_file in bam_files:
        # Extract the prefix 'ipsc179' from the filename
        #prefix = bam_file.split('_output.sorted.reheaded.bam')[0]
        prefix = bam_file.split('.bam')[0]

        variant_path = os.path.join(output_dir, f"{prefix}.vcf")
        command = f"""
singularity run --cleanenv --nv \\
    --bind /home/$USER:/home/$USER,/scratch/$USER:/scratch/$USER,/projectsp/XXX:/projectsp/XXX,/projects:/projects,/cache:/cache \\
    XXXX/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \\
    --ref XXX/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \\
    --in-bam {os.path.join(input_dir, bam_file)} \\
    --out-variants {variant_path} \\
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10
"""
        commands.append(command)

    # Create the SLURM script content
    slurm_content = f"""#!/bin/bash
#SBATCH --job-name=parabricks_haplotypecaller
#SBATCH --partition=gpu
#SBATCH --output={output_dir}/parabricks_haplotypecaller_%j.out
#SBATCH --error={output_dir}/parabricks_haplotypecaller_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:8
#SBATCH --mem=80G
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --exclude='cuda[001-008]'

#srun -p gpu --gres=gpu:8 --mem=80G --cpus-per-task=12 --time=4:00:00 --exclude='cuda[001-008]' --pty zsh

module load singularity/3.8.3

# Set environment variables for CUDA
export CUDA_VISIBLE_DEVICES=0,1,2,3

# Execute each command sequentially
{''.join(commands)}
    """

    # Write the script to a file
    with open(output_filename, 'w') as file:
        file.write(slurm_content)

    print(f"SLURM script '{output_filename}' has been created.")

def main():
    parser = argparse.ArgumentParser(description='Generate a SLURM script for running HaplotypeCaller using Parabricks.')
    parser.add_argument('--input_dir', type=str, required=True, help='Directory containing input BAM files')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to output the VCF files')
    parser.add_argument('--output_filename', type=str, default='99.parabricks_haplotypecaller_job.slurm', help='Filename for the generated SLURM script')

    args = parser.parse_args()

    generate_slurm_script(args.input_dir, args.output_dir, args.output_filename)

if __name__ == '__main__':
    main()
