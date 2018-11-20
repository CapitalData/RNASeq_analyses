#!bin/bash




mkdir RNASEQ_data

stardir='RNASEQ_data/star_REP'
rsemdir='RNASEQ_data/rsem_REP'
mkdir $stardir
mkdir $rsemdir

# get the list of FASTQ files.


ARRAY=(`find data/RAW_DATA/ -name "*fastq.gz"  -printf %f!`)
ParentPATH="data/RAW_DATA/Project_Docker_works/"
touch ProcessedRSEM.txt

ARRAY=(`ls $ParentPATH`)

#For each path, run STAR and RSEM
for i in ${ARRAY[*]}; do

  echo $ParentPATH$i
  ARRAY2=(`find $ParentPATH$i -name "*fastq.gz"`)
  NAME=`basename ${ARRAY2[0]} .fastq.gz`
  echo "making directory: RNASEQ_data/star_REP"

  # RUN THE STAR ALIGNMENT
  opt/STAR-2.5.3a/bin/Linux_x86_64/STAR --genomeDir data/genome_ref --readFilesCommand zcat \
  --readFilesIn ${ARRAY2[*]} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
  --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM \
  --runThreadN 20 --outFileNamePrefix "RNASEQ_data/star_REP/"$NAME;

   # RUN THE RSEM ALIGNMENT	
  opt/RSEM-1.3.0/rsem-calculate-expression --bam --no-bam-output -p 20 \
  --paired-end --forward-prob 0 \
  "RNASEQ_data/star_REP/"$NAME"Aligned.toTranscriptome.out.bam" \
  "RNASEQ_data/rsem_REP/"$NAME"_rsem.log"

echo $NAME >> ProcessedRSEM.txt

done
