#!/bin/bash
#SBATCH --job-name=CombineTable
#SBATCH --ntasks=16
#SBATCH --mem=120000
#SBATCH --output=/work/sgrissom/Scripts/Output/CombineTable_%A.o
#SBATCH --error=/work/sgrissom/Scripts/Error/CombineTable_%A.e
#SBATCH --mail-user=sgrissom@udel.edu
#SBATCH --mail-type=ALL


R1_file="/work/sgrissom/SampleFiles/Bulk-Population-2_R1_001.fastq" 
R2_file="/work/sgrissom/SampleFiles/Bulk-Population-2_R2_001.fastq"
working_directory="/work/sgrissom/WorkingDirectory"
output_directory="/work/sgrissom/OutputDirectory"
genome_fasta_file="/work/sgrissom/ReferenceFiles/CHO_fasta/Cricetulus_griseus_picr.CriGri-PICR.dna.nonchromosomal.fa"
genome_STAR_file_base="/work/sgrissom/ReferenceFiles/starIndex"
genome_gtf_file="/work/sgrissom/ReferenceFiles/CHO_GTF/Cricetulus_griseus_picr.CriGri-PICR.104.gtf"

#Set headers to geneID and sample name & drop GeneID:
#headers="geneID\t$sample_name"
#htseq_count_output_tmp="$htseq_output_dir/${sample_name}_count_tmp.txt"
#htseq_count_output_mod="$htseq_output_dir/${sample_name}_count_mod.txt"
#cut -f 2- -d ':' $htseq_count_output > $htseq_count_output_tmp
#echo -e $headers | cat - $htseq_count_output_tmp > $htseq_count_output_mod



for i in {1..40}
  do 
  j=$(expr $i - 1)
  FileName="${working_directory}/htseq_counts/Single-Cell-${i}_R1_count_mod.txt"
  tmpName="${working_directory}/htseq_counts/tmp_${i}.txt"
  tmpNameOld="${working_directory}/htseq_counts/tmp_${j}.txt"
   if (( $i==1))
    then  
     cp ${FileName} ${tmpName} 
   else
     cut -f 2 ${FileName} | paste ${tmpNameOld} - > ${tmpName}
   fi
 done
FinalTable="${working_directory}/htseq_counts/SingleSample_FinalTable.txt"
mv $tmpName $FinalTable
rm ${working_directory}/htseq_counts/tmp*
