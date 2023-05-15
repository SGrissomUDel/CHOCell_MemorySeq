#!/bin/bash
PATH=$PATH:/usr/local/FastQC

while getopts "i:j:g:h:t:w:o:" opt;
do
  case $opt in
    i ) R1_file=$OPTARG ;;
    j ) R2_file=$OPTARG ;;
    g ) genome_fasta_file=$OPTARG ;;
    h ) genome_STAR_file_base=$OPTARG ;;
    t ) genome_gtf_file=$OPTARG ;;
    w ) working_directory=$OPTARG ;;
    o ) output_directory=$OPTARG ;;
  esac
done

##Extract sample names from R1/R2 file names
tmp_name=$(basename $R1_file)
sample_name="${tmp_name%_*}"

#Set output of trimgalore to the working directory input
trimgalore_output="$working_directory"

#Set file basename to sample_name (i.e. MW## from line 49-50)
#Trimgalore uses this as the "basename" for all file names
basename_tg="$sample_name"

echo "Starting Trim"
#trim_galore commands using the following flags:
trim_galore --paired --illumina --cores 3 --fastqc -o $trimgalore_output \
--basename $basename_tg \
$R1_file $R2_file
echo "Finished Trim"

#Prepare output filenames - use trimgalore defaults   
Trim_output_1="$working_directory/${sample_name}_val_1.fq.gz"
Trim_output_2="$working_directory/${sample_name}_val_2.fq.gz"

#Make directory to copy all output gtf files into so it's easier to access
fastqc_dir="$output_directory/fastqc_files" 
mkdir -p $fastqc_dir

#Give fastqc html files a variable name
fastqc_1="$working_directory/${sample_name}_val_1_fastqc.html"
fastqc_2="$working_directory/${sample_name}_val_2_fastqc.html"

#Copy all fastqc reports to single directory
cp $fastqc_1 $fastqc_dir 
cp $fastqc_2 $fastqc_dir

#Prepare STAR output variables
STAR_output_dir="${working_directory}/STAR_${sample_name}"
mkdir -p $STAR_output_dir
STAR_output_prefix="${STAR_output_dir}/${sample_name}_"
STAR_align="Aligned.out.bam"
STAR_output_prefix_bam="${STAR_output_prefix}${STAR_align}"
#This needs to go for --outFileNamePrefix line
#Software uses a directory output to dump all output files generated into
#User supplies a file prefix name aka "sample_name"
#Run STAR alignment using defined variables
echo "Running STAR: Aligning..."
STAR --runThreadN 12 --genomeDir $genome_STAR_file_base \
--readFilesIn <(gunzip -c $Trim_output_1 $Trim_output_2) \
--outFileNamePrefix $STAR_output_prefix \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outSAMattrIHstart 0 \
--outSAMstrandField intronMotif
#--readFilesIn $Trim_output_1 $Trim_output_2
#--readFilesIn <(gunzip -c $Trim_output_1 $Trim_output_2)
#--readFilesCommand zcat
echo "Finished alignment"
rm $Trim_output_1
rm $Trim_output_2


samtools_sorted_bam_output="$working_directory/$sample_name.sorted.bam"
echo "Starting Samtools"
#Run samtools using "samtools view" and "samtools sort" 
samtools sort -n $STAR_output_prefix_bam -o $samtools_sorted_bam_output
echo "Finished Samtools"
rm -R $STAR_output_dir

##Prepare HTseq output file
htseq_output_dir="$working_directory/htseq_counts"
mkdir -p $htseq_output_dir
htseq_count_output="$htseq_output_dir/${sample_name}_count.txt"
#Note: must put {} around "sample_name" since it is followed by "_count" which...
#...is considered a valid variable name character so it will try to look for...
#...the variable "sample_name_count" which does not exist
#Note: variable names must start with letters or _ not numbers 
#Run HTseq
echo "running htseq"
htseq-count -f bam -r name --stranded=yes -t exon -i gene_id \
$samtools_sorted_bam_output $genome_gtf_file > $htseq_count_output
echo "counts generated"
#Edit htseq output count table to add on header row with "geneID" and...
#...$sample_name labels for use in DESeq2
#...Also removes "GeneID:" from each row name in htseq_count_output...
#...to work with DESeq2; saves in new file "htseq_count_output_mod"

#Set headers to geneID and sample name & drop GeneID:
headers="geneID\t$sample_name"
htseq_count_output_tmp="$htseq_output_dir/${sample_name}_count_tmp.txt"
htseq_count_output_mod="$htseq_output_dir/${sample_name}_count_mod.txt"
cut -f 2- -d ':' $htseq_count_output > $htseq_count_output_tmp
echo -e $headers | cat - $htseq_count_output_tmp > $htseq_count_output_mod
rm $samtools_sorted_bam_output
rm $htseq_count_output_tmp
rm $htseq_count_output
