#!/bin/bash

# go to dir with .gtf files
genome_dir=$(find /docker_home -name "genome_ref")
cd $genome_dir

# creating a temp file as a copy of the original gtf file in case the conversion fails
gtf_filename=$(ls | grep "gtf")
temp_filename="temp_gtf_file.gtf"
cp "$gtf_filename" "$temp_filename"

grep -e "\tgene\t" Homo_sapiens.GRCh38.86.gtf | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$6,".",$4,$10,$12,$14 }' | \
sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.86.gene.bed

grep -e "\ttranscript\t" Homo_sapiens.GRCh38.86.gtf | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$10,".",$4,$14,$16,$18 }' | \
sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.86.transcript.bed

# checks for conversion failure by checking the size of the bed file, then recreates the original gtf file from the temp file
minimumsize=90000
actualsize=$(wc -c <"Homo_sapiens.GRCh38.86.transcript.bed")

if [ $actualsize -ge $minimumsize ]; then
    rm "$temp_filename"
else
	mv "$temp_filename" "$gtf_filename"
    echo "Conversion failed. Fix script and rerun."
fi