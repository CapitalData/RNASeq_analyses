#!/bin/bash

source ./dir_locations.sh

$(find /docker_home -name "rsem-prepare-reference") --gtf "${genome_ref_dir}/${genome_name}.86.gtf" \
    --star \
    --star-path "${star_bin_dir}" \
    -p 16 \
    "${genome_ref_dir}/${genome_name}.dna.primary_assembly.fa" \
    "${genome_ref_dir}/HomoSapStar_ref" &&\
  echo RSEM reference prepped successfully.
