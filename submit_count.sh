# This script aims to submit the jobs of generating the files using CellRanger 
# (script: cellranger_count.sh), including feature-bc matrices, h5 files, 
# and bam files.


mkdir log # used to save the running log
while read NAME # the name of a sample saved in the 'sample.txt' file
do
   echo $NAME
   #qsub -v "NAME=$FASTQ1" -N "$FASTQ1" run_chardo.sh 
   sbatch --job-name=$NAME.run --output=./log/$NAME.out --export=NAME=$NAME cellranger_count.sh
   # submit the job to system
   
   ##############################################################
   # 
   # -J, --job-name=<jobname>
   # Specify a name for the job allocation
   # 
   # -o, --output=<filename pattern>
   # Instruct Slurm to connect the batch script's standard output 
   # directly to the file name specified in the "filename pattern"
   # 
   # --export=<[ALL,]environment variables|ALL|NONE>
   # Identify which environment variables from the submission 
   # environment are propagated to the launched application
   # 
   ##############################################################

   sleep 0.1s # delay for a specified amount of time
done < sample.txt
# a txt file of the list of sample names, each of which corresponds to a directory
# containing fastq files for pair-ended sequencing, e.g., 
#
#################################################################
# 
# S28_CKDL210002750-1a-SI_TT_B2_H2NV5DSX2_S1_L002_I1_001.fastq.gz
# S28_CKDL210002750-1a-SI_TT_B2_H2NV5DSX2_S1_L002_I2_001.fastq.gz
# S28_CKDL210002750-1a-SI_TT_B2_H2NV5DSX2_S1_L002_R1_001.fastq.gz
# S28_CKDL210002750-1a-SI_TT_B2_H2NV5DSX2_S1_L002_R2_001.fastq.gz
# 
#################################################################
