# this script aims to perform alignment based on fastq files to generate a
# series of files acceptable by Seurat, including feature-bc matrices, h5 files,
# and bam files, etc.
# 
# The script submit_count.sh calls this script.


#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=16:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks=16
#SBATCH --mem=128GB

wd=/fs/scratch/PAS1475/liyang/Katherine # the current working directory

CellRanger=/fs/ess/PCON0022/liyang/scRNA-Seq_tools/cellranger
# the location of CellRanger
FastqFolder=/fs/scratch/PAS1475/liyang/Katherine/fastq
# the directory of all the fastq files, which may consist of a number of
# subdirectories (see the following example)
# 
#################################################################
# 
# S28_CKDL210002750-1a-SI_TT_B2_H2NV5DSX2_S1_L002_I1_001.fastq.gz
# S28_CKDL210002750-1a-SI_TT_B2_H2NV5DSX2_S1_L002_I2_001.fastq.gz
# S28_CKDL210002750-1a-SI_TT_B2_H2NV5DSX2_S1_L002_R1_001.fastq.gz
# S28_CKDL210002750-1a-SI_TT_B2_H2NV5DSX2_S1_L002_R2_001.fastq.gz
# 
#################################################################

# GRCh38_Refer=/fs/project/PCON0005/BISR_Tools/refdata-cellranger-GRCh38-3.0.0
MM10_Refer=/fs/project/PAS1475/liyang/tools/genome/refdata-gex-mm10-2020-A
# the location of the reference files for an orrganism, e.g., genome and genes

cd $wd

${CellRanger} count --id=$NAME --transcriptome=${MM10_Refer} --fastqs=${FastqFolder}/$NAME\
--sample=$NAME --localcores=16 --localmem=128
# The cellranger count pipeline aligns sequencing reads in FASTQ files to a reference
# transcriptome and generates a .cloupe file for visualization and analysis in
# Loupe Browser, along with a number of other outputs compatible with other
# publicly-available tools for further analysis.
# 
# --id            : the name of a ample
# --transcriptome : Path of folder containing 10x-compatible transcriptome reference
# --fastqs        : Path to input FASTQ data
# --sample        : Prefix of the filenames of FASTQs to select, i.e., sample name
# --localcores    : Set max cores the pipeline may request at one time. Only applies
#                   to local jobs
# --localmem      : Set max GB the pipeline may request at one time. Only applies to local jobs
