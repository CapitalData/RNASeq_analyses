FS=!
#ARRAY=(`find /data/RAW_DATA/ -name "*fastq.gz"  -printf %f!`)
ParentPATH="/data/RAW_DATA/Project_16494_FUID1031578/"

ARRAY=(`ls $ParentPATH`)

for i in ${ARRAY[*]}; do
    
    echo $ParentPATH$i
    ARRAY2=(`find $ParentPATH$i -name "*fastq.gz"`)
    NAME=`basename ${ARRAY2[0]} .fastq.gz`
   echo "making directory: RNASEQ_data/star_REP"
   #echo "${ARRAY2}"
   mkdir "RNASEQ_data/star_REP"
   printf -v var "%s\n" "${ARRAY2[@]}"
   
  # RUN THE STAR ALIGNEMNT
   /home/jabust/STAR-2.6.0a/bin/Linux_x86_64/STAR --genomeDir /home/jabust/genome --readFilesCommand zcat \
  --readFilesIn ${ARRAY2[*]} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
  --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM \
  --runThreadN 20 --outFileNamePrefix "RNASEQ_data/star_REP/"$NAME;

done
