# This script aims to move CellRanger output files from 'outs' directories in the
# directory of each sample to a common directory, e.g., 'Cellranger_outs'


library(tidyverse) # The tidyverse is an opinionated collection of R packages
# designed for data science.

working_dir <- "/fs/scratch/PAS1475/liyang/Katherine"
setwd(working_dir)
all_dirs <- list.dirs(recursive = F) # list all the subdirectories in the current
# directory, i.e., the directories for all the samples.
dir.create("cellranger_outs", showWarnings = F) # the directory to save all the
# files output by CellRanger outputs

dir.create("bam_outs", showWarnings = F) # the directory to save the bam files

# i = 1
for (i in seq_along(all_dirs)) {
  this_dir <- all_dirs[i]
  if (any(str_detect(list.dirs(this_dir, recursive = F), "outs"))) { # find the
  # 'outs' directory
    destination_dir <- paste("cellranger_outs", this_dir, sep = "/")
    bam_dir <- paste("bam_outs", this_dir, sep = "/")
    dir.create(destination_dir, showWarnings = F) # create the subdirectory
    web_summary_file <-
      paste(this_dir, "outs", "web_summary.html", sep = "/") # web summary report
    
    h5_file <-
      paste(this_dir, "outs", "filtered_feature_bc_matrix.h5", sep = "/")
      # the feature-bc matrix in binary encoding format
      
    loupe_file <- paste(this_dir, "outs", "cloupe.cloupe", sep = "/")
    # Loupe Browser is a desktop application that provides interactive visualization
    # functionality to analyze data from different 10x Genomics solutions.
    
    bc_matrix_file <-
      paste(this_dir, "outs", "filtered_feature_bc_matrix", sep = "/")
      # the feature-bc matrix file
      
    bam_file <-
      paste(this_dir, "outs", "possorted_genome_bam.bam", sep = "/")
      # the alignment result file in binary encoding
      
    bam_index_file <-
      paste(this_dir, "outs", "possorted_genome_bam.bam.bai", sep = "/")
      # IGV requires that both SAM and BAM files be sorted by position and indexed,
      # and that the index files follow a specific naming convention
      
    # Copy the files for Seurat
    file.copy(web_summary_file, destination_dir)
    file.copy(h5_file, destination_dir)
    file.copy(loupe_file, destination_dir)
    file.copy(bc_matrix_file, destination_dir, recursive = T)
    
    # Copy BAM files
    file.copy(bam_file, paste0("bam_outs/", this_dir, ".bam"))
    file.copy(bam_index_file, paste0("bam_outs/", this_dir, ".bam.bai"))
  }
}
