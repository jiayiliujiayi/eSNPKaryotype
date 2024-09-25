# library

library(zoo)
library(gplots)
library(stringr)
library(data.table)
library(optparse)

# Parse options
## Set up command line argument parsing
option_list = list(
  make_option(c("--vcfdir"), type="character", default=NULL, help="Directory containing vcf files", metavar="directory")
)

## Parse the command line arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## Check if vcfdir argument was provided
if (is.null(opt$vcfdir)) {
  stop("Error: Please specify the vcf directory using --vcfdir", call. = FALSE)
} else {
    vcfdir = opt$vcfdir
}

## Print the provided directory
print(paste("VCF directory:", opt$vcfdir))

## Here you could add code to process vcf files found in the directory
# For instance, listing vcf files:
vcf_files = list.files(path = opt$vcfdir, pattern = "\\.vcf$", full.names = TRUE)
print("VCF files found:")
print(vcf_files)

# function sources

function_list_to_source <-
  list.files(path = "./source",
             pattern = "*.R",
             full.names = T)
for (i in 1:length(function_list_to_source)) {
  source(function_list_to_source[i])
}

# PATHS

bin = "XXX" ## set your bin
CMD_samtools = sprintf("%s/samtools", bin)
CMD_java = sprintf("%s/java", bin)
Picard_Path <- "XXX/picard-2.27.5-0"
CMD_gatk <- sprintf("%s/gatk", bin)

Genome_Fa <- "XXX/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
DIR_Genome_DICT = "XXX/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.dict"
Organism <- "Human"

# DIRs

temp_dir <- "./temp_dir"
# test if temp_dir exists
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir)
}

plot_output_directory <- "./plot_output"
if (!dir.exists(plot_output_directory)) {
  dir.create(plot_output_directory)
}

### DELETION TABLE function needs a BAMDIR!

# make a eSNPKaryotyping-compatible dbSNP file (one-off) #####
# ! Use $Databases/Genomes/hg38/common_snp150_by_chr ! #####
#Edit_dbSNP_Files(Directory = "snp",
#                 File_Name = "snp151Common",
#                 Organism = "Human")

# loop through vcf list

vcf_file_list <-
  list.files(path = opt$vcfdir,
             full.names = F, 
             recursive = F,
             pattern = "\\.vcf$")
base_names = gsub("\\.vcf", "", vcf_file_list)

k <- 1
for (k in 1:length(vcf_file_list))
    {
    print(sprintf("Started %s at %s ", base_names[k], date()))

    ###########################
    ###      Edit VCF       ###
    ###########################
    
    EditVCF(vcfdir = opt$vcfdir, 
            base_name = base_names[k], 
            Organism = "Human", 
            out_dir = "./temp_dir") 
    
    ###########################
    ###  Read VCF table in ###
    ###########################
    
    VCF_table <-
    read.delim(file = sprintf('%s/%s.csv',temp_dir, base_names[k]), 
               header = T, sep = "\t",
               quote = "",
               dec = ".")
    VCF_table$chr <- as.numeric(VCF_table$chr)
    VCF_table <- VCF_table[order(VCF_table$chr,
                                 VCF_table$position), ]
    VCF_table <- VCF_table[VCF_table$chr > 0, ]
    
    #######################################
    ###  Return MajorMinorCalc results  ###
    #######################################
    
    MajorMinorCal_results <-
    MajorMinorCalc(Table = VCF_table,
                   minDP = 20, # ! change back to 20 later !
                   maxDP = 1000,
                   minAF = 0.2)
    
    ###########################
    ###  set the file name  ###
    ### and title for the PDF##
    ###########################
    
    pdf_series_name <- base_names[k]
    
    ###########################################
    #####   Plot Allelic ratio along the  #####
    #####  genome for duplication detection ###
    ###########################################
    ## noxy
    pdf(file = sprintf("%s/%s_genome_noxy.pdf", 
                       plot_output_directory, base_names[k]), 
        width=12, height=5, 
        paper = "USr",
        title = base_names[k])
    try({
        PlotGenome_noxy(orderedTable = MajorMinorCal_results,
                   Window = 151,
                   Ylim = 3,
                   PValue = TRUE,
                   Organism = Organism)
        dev.off()
    })
    ## withxy
    pdf(file = sprintf("%s/%s_genome.pdf", 
                       plot_output_directory, base_names[k]), 
        width=12, height=5, 
        paper = "USr",
        title = base_names[k])
    try({
        PlotGenome(orderedTable = MajorMinorCal_results,
                   Window = 151,
                   Ylim = 3,
                   PValue = TRUE,
                   Organism = Organism)
        dev.off()
    })
    
    tbl_DeletionTable_output <-
    DeletionTable(Directory = "./temp_dir",
                  Table = MajorMinorCal_results,
                  base_name = base_names[k], 
                  dbSNP_Data_Directory = "./snp",
                  dbSNP_File_Name = "Edited_snp151Common_chr",
                  Genome_Fa_dict = Genome_Fa,
                  Organism = "Human",
                  source_bam = sprintf("./out/%s.bam", 
                                       base_names[k] 
                                      ),
                  temp_dir = temp_dir)
    # PATH_tbl_deletion = sprintf("%s/%s.txt", "./temp_dir", 
    #                     base_names[k])
    # fwrite(tbl_DeletionTable_output, 
    #        PATH_tbl_deletion, 
    #        quote = F, col.names = T, row.names = F, sep = "\t")
    
    #######################
    #####   Plot LOH  #####
    #######################
    ## noxy
    pdf(file = paste(plot_output_directory,
                     '/',
                   base_names[k],
                     "_zygosity_blocks_noxy.pdf", sep = ""),
        title = base_names[k])
    Plot_Zygosity_Blocks_noxy(tbl_DeletionTable_output, 
                         Window = 1500000, Max = 6, Max2 = 60, Organism = "Human")
    dev.off()
    
    ## withxy
    pdf(file = paste(plot_output_directory,
                     '/',
                   base_names[k],
                     "_zygosity_blocks.pdf", sep = ""),
        title = base_names[k])
    Plot_Zygosity_Blocks(tbl_DeletionTable_output, 
                         Window = 1500000, Max = 6, Max2 = 60, Organism = "Human")
    dev.off()
    # tbl_DeletionTable_output$snp = 
    # write.table(
    #     sprintf("%s", ), 
    #     quote = F, 
    #     col.names = T, row.names = F, 
    #     sep = "\t"
    # )
    print(sprintf("Finished %s at %s ", base_names[k], date()))
    
    }
