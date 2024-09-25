# Load necessary libraries
library(data.table)
library(optparse)

# Parse command line options
option_list = list(
  make_option(c("--bamdir"), type="character", default=NULL, help="Directory containing BAM files", metavar="directory")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check if bamdir argument was provided
if (is.null(opt$bamdir)) {
  stop("Error: Please specify the BAM directory using --bamdir", call. = FALSE)
} else {
    bamdir = opt$bamdir
}
temp_dir = "./temp"

# Define the path to samtools binary (update with your actual path)
bin = "XXXX"  # Set your bin directory
CMD_samtools = sprintf("%s/samtools", bin)

# List BAM files in the provided directory
bam_file_list <- list.files(path = opt$bamdir, full.names = F, pattern = "\\_output.sorted.bam$")
base_names = gsub("\\_output.sorted.bam", "", bam_file_list)

# Prepare sed for creating new header file using chromosome aliases
aliases <- fread("source/hg38.chromAlias.txt")

# Prepare sed replacements for alias to sequenceName and sequenceName to V1
aliases[, alias_to_seq := paste0("s/SN:", `alias names`, "/SN:", `# sequenceName`, "/g")]
aliases[, seq_to_v1 := paste0("s/SN:", `# sequenceName`, "/SN:", V1, "/g")]

# Specific replacements for X, Y, and other chromosomes
specific_replacements <- c(
  "s/SN:X/SN:chrX/g",
  "s/SN:Y/SN:chrY/g", 
  "s/SN:KI270752.1/SN:chrUn_KI270752v1/g"
)

# Combine all sed commands
all_sed_commands <- c(aliases$alias_to_seq, aliases$seq_to_v1, specific_replacements)

# Write the sed commands to a file
writeLines(all_sed_commands, "./temp/update_header.sed")

# Loop through BAM files to rehead and index
for (k in 1:length(bam_file_list)) {
  print(sprintf("Started reheading %s at %s", base_names[k], date()))
  
  # Extract the current BAM header into a file
  header_file = sprintf("%s/%s_header.sam", temp_dir, base_names[k])
  system(sprintf("%s view -H %s/%s > %s", CMD_samtools, bamdir, bam_file_list[k], header_file))
  
  # Edit the header file using sed
  system(sprintf("sed -i -f ./%s/update_header.sed %s", temp_dir, header_file))
  
  # Rehead the BAM file using the edited header
    system(sprintf("%s reheader %s %s/%s > %s/%s_output.sorted.reheaded.bam", CMD_samtools, header_file, bamdir, bam_file_list[k], bamdir, base_names[k]))
  
  # Index the reheaded BAM file
  system(sprintf("%s index -@ 20 %s/%s_output.sorted.reheaded.bam", CMD_samtools, bamdir, base_names[k]))

    ## Print finished message
    print(sprintf("Finished reheading %s at %s", base_names[k], date()))

    
}
