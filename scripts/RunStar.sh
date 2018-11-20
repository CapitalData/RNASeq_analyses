#!/bin/bash

# TODO: I left this here because I don't know what it does. Figure that out.
FS=!
readonly usage="
Usage: RunRSEM.sh -p project_name
  
  Runs STAR alignment on all the .fastq.gz files in the ./data/RAW_DATA/project_name directory
    
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
mkdir -p "$star_dir"

# get the list of FASTQ filepaths.
fastq_filepaths=(`find "$project_path" -name "*.fastq.gz"`)

# for each FASTQ file, run STAR and RSEM
for filename in "${fastq_filepaths[@]}"; do
  # using regex to parse the name of the analysis from the filename
  analysis_name="$(expr "$filename" : '.*/\(.*\).fastq.gz')"
  echo "Processing analysis $analysis_name"

  # RUN THE STAR ALIGNMENT
  # TODO: pass genome dir as an argument to script?
  /home/jabust/STAR-2.6.0a/bin/Linux_x86_64/STAR --genomeDir /home/jabust/genome --readFilesCommand zcat \
  --readFilesIn "${filename}" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
  --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM \
  --runThreadN 20 --outFileNamePrefix "${star_dir}/${analysis_name}";

done
