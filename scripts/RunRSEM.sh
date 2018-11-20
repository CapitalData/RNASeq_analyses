#!/bin/bash

readonly usage="
Usage: RunRSEM.sh -p project name
  
  Runs first STAR, then RSEM Alignment on all the .fastq.gz files 
  in the ./data/RAW_DATA/project_name directory.

    -p         project name: name of the project you're working on.
"

# parsing provided arguments
while getopts p: opt; do
    case "$opt" in
        p)   project_name="$OPTARG";;
        '?') echo "$usage" >&2; exit 1;;
        ':') echo "option $OPTARG requires an argument" >&2; exit 1;;
        '') echo 
    esac
done
shift "$(($OPTIND-1))"

if [[ -z $project_name ]]; then
  echo "project_name argument must be provided."
  echo "$usage"
  exit 1
else
  project_path="./data/RAW_DATA/${project_name}"
fi

star_dir="./data/RNASEQ_data/star_REP"
rsem_dir="./data/RNASEQ_data/rsem_REP"
genome_dir="./data/genome_dir"
mkdir -p "$star_dir" "$rsem_dir"

# get the list of FASTQ filepaths.
fastq_filepaths=(`find "$project_path" -name "*.fastq.gz"`)

# for each FASTQ file, run STAR and RSEM
for filename in "${fastq_filepaths[@]}"; do
  # using regex to parse the name of the analysis from the filename
  analysis_name="$(expr "$filename" : '.*/\(.*\).fastq.gz')"
  echo "Processing analysis $analysis_name"

  # RUN THE STAR ALIGNMENT
  # TODO: pass genome dir as an argument to script?
  /docker_home/opt/STAR-2.5.3a/bin/Linux_x86_64/STAR --genomeDir "${genome_dir}" --readFilesCommand zcat \
  --readFilesIn "${filename}" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
  --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM \
  --runThreadN 20 --outFileNamePrefix "${star_dir}/${analysis_name}";

  # RUN THE RSEM ALIGNMENT
  # TODO: pass HomoSapStar_ref dir as an argument to script?
  # TODO: is that Aligned filename right? Should there be a '.' before 'Aligned'?
  /docker_home/opt/RSEM-1.3.0/rsem-calculate-expression --bam --no-bam-output -p 20 \
  --paired-end --forward-prob 0 \
  "${star_dir}/${analysis_name}Aligned.toTranscriptome.out.bam" \
  "ref/HomoSapStar_ref" "${rsem_dir}/${analysis_name}_rsem" >& \
  "${rsem_dir}/${analysis_name}_rsem.log"

  # writing the analyzed filename to another file
  # TODO: figure out why this is done
  echo "${analysis_name}" >> ProcessedRSEM.txt

done
