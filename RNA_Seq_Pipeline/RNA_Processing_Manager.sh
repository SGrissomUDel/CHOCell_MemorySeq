#!/bin/bash
#SBATCH --job-name=RNASEQProcess
#SBATCH --ntasks=16
#SBATCH --mem=120000
#SBATCH --output=/work/sgrissom/Scripts/Output/RNASEQ_%A_%a.o
#SBATCH --error=/work/sgrissom/Scripts/Error/RNASEQ_%A_%a.e
#SBATCH --mail-user=sgrissom@udel.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-40


#To submit the job array, use sbatch --array=1-x

#sample_file="path/to/sample_file.txt"
sample_file="/work/sgrissom/SampleFiles/FileNames.txt" 
working_directory="/work/sgrissom/WorkingDirectory"
output_directory="/work/sgrissom/OutputDirectory"
genome_fasta_file="/work/sgrissom/ReferenceFiles/CHO_fasta/Cricetulus_griseus_picr.CriGri-PICR.dna.nonchromosomal.mAb.fa"
genome_STAR_file_base="/work/sgrissom/ReferenceFiles/STAR_Index_Folder/mAb_STAR"
genome_gtf_file="/work/sgrissom/ReferenceFiles/CHO_GTF/Cricetulus_griseus_picr.CriGri-PICR.104.mAb.gtf"


#-------------------------Sample List--------------------------------------
#---Input raw read files are saved into text file list
#---Text file list consists of /path/filename for R1 and R2 fasta files...
#---...for paired end reads
#---Reads are saved on /zfs/blennermeep1/CHO_Projects/128.120.88.251/C202SC19060126/raw_data/
#NEW:
#/path/ZN_D0_1_1.fq.gz
#/path/ZN_D0_1_2.fq.gz
#OLD:
#/path/MW#(#)_S##_R1_001.fastq.gz
#/path/MW#(#)_S##_R2_001.fastq.gz

#---Calculate line number for R1 and R2 /path/filename in sample_list.txt
#---expr does math; PBS_ARRAY_Index is number of job loop in array 
#---"\*" indicates multiplication and - is minus
R1_line_number=$(expr $SLURM_ARRAY_TASK_ID "*" 2 - 1)
R2_line_number=$(expr $SLURM_ARRAY_TASK_ID "*" 2)
#---Extract correct line number from sample_list.txt file
#---aka if on first loop of job - running sample 1 in list this...
#---...will pull out line 1 for R1 file and line 2 for R2 file
#---"head -1" indicates display one line (-n indicates number of lines)
R1_file=$(tail -n+$R1_line_number $sample_file | head -1)
R2_file=$(tail -n+$R2_line_number $sample_file | head -1)
echo ${SLURM_ARRAY_TASK_ID}
echo $R1_file
echo $R2_file


#Command to run the script
chmod u+r+x /work/sgrissom/Scripts/Trim_to_HTSeq.sh
#while getopts "i:j:g:h:t:w:o:" opt;
#do
#  case $opt in
#    i ) R1_file=$OPTARG ;;
#    j ) R2_file=$OPTARG ;;
#    g ) genome_fasta_file=$OPTARG ;;
#    h ) genome_STAR_file_base=$OPTARG ;;
#    t ) genome_gtf_file=$OPTARG ;;
#    w ) working_directory=$OPTARG ;;
#    o ) output_directory=$OPTARG ;;
#  esac
#done
/work/sgrissom/Scripts/Trim_to_HTSeq.sh -i $R1_file -j $R2_file -g $genome_fasta_file -h $genome_STAR_file_base -t $genome_gtf_file -w $working_directory -o $output_directory
