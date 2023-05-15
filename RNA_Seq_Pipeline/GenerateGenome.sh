#!/bin/bash
#SBATCH --job-name=CHO_Genome_Index
#SBATCH --ntasks=16
#SBATCH --mem=120000
#SBATCH --output=STAR_Index.o
#SBATCH --error=STAR_Index.e
#SBATCH --mail-user=sgrissom@udel.edu
#SBATCH --mail-type=ALL


echo "Running STAR: Generating Genome Index ..."
STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir /work/sgrissom/ReferenceFiles/STAR_Index_Folder/mAb_STAR \
--genomeFastaFiles /work/sgrissom/ReferenceFiles/CHO_fasta/Cricetulus_griseus_picr.CriGri-PICR.dna.nonchromosomal.mAb.fa \
--sjdbGTFfile /work/sgrissom/ReferenceFiles/CHO_GTF/Cricetulus_griseus_picr.CriGri-PICR.104.mAb.gtf \
--sjdbOverhang 149
echo "Finished Generating Genome Index"
