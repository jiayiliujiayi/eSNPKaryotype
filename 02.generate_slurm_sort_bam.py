import os

def generate_slurm_script(out_dir="./temp", job_name="sort_bam", partition="XXX", cpus=8, mem="50G", time="2-00:00:00"):
    # Create a list of BAM files that need to be sorted
    bam_files = [f for f in os.listdir(out_dir) if f.endswith("_output.bam")]
    
    # Begin SLURM script
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --output={out_dir}/{job_name}_%j.out
#SBATCH --error={out_dir}/{job_name}_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}

module load samtools
module load bowtie2

# Execute each command sequentially
"""
    
    # Add sorting commands for each BAM file
    for bam_file in bam_files:
        sorted_bam_file = bam_file.replace("_output.bam", "_output.sorted.bam")
        slurm_script += f'\necho "Processing {bam_file}"\n'
        slurm_script += f'samtools sort -o {out_dir}/{sorted_bam_file} {out_dir}/{bam_file}\n'
    
    return slurm_script

# Usage example
slurm_script = generate_slurm_script()
with open("02.sort_bam.slurm", "w") as f:
    f.write(slurm_script)

print("SLURM file generated successfully!")
