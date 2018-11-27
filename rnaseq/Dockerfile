# adapted from https://github.com/CapitalData/gtex-pipeline/blob/master/rnaseq/Dockerfile

# RNAseq Dockerfile 

# use an official parent image
FROM ubuntu:16.04

RUN mkdir -p /docker_home/opt

# set the working directory
WORKDIR /docker_home

# install ubuntu tooling
RUN apt-get update &&\
  apt-get install -y software-properties-common &&\
  add-apt-repository -y ppa:openjdk-r/ppa &&\
  apt-get update &&\
  apt-get install -y \
  build-essential \
  cmake \
  curl \
  libboost-all-dev \
  libbz2-dev \
  libcurl3-dev \
  liblzma-dev \
  libncurses5-dev \
  libssl-dev \
  openjdk-8-jdk \
  r-base-core=3.2.3-4 \
  unzip \
  vim-common \
  wget \
  zlib1g-dev &&\
  rm -rf /var/lib/apt/lists/*

# install R packages
RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("qvalue", "sva", "edgeR", "cluster", "Biobase"));'
RUN wget https://cran.r-project.org/src/base/R-3/R-3.2.5.tar.gz &&\
  tar -xf R-3.2.5.tar.gz && rm R-3.2.5.tar.gz &&\
  cd R-3.2.5 &&\
  ./configure --with-x=no &&\
  cd src/nmath/standalone &&\
  make

# create RMATH env variable
ENV RMATH /R-3.2.5/src

#-----------------------------
# Pipeline components
#-----------------------------

WORKDIR /docker_home/opt
# htslib
RUN wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2 &&\
  tar -xf htslib-1.8.tar.bz2 && rm htslib-1.8.tar.bz2 &&\
  cd htslib-1.8 &&\
  ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs &&\
  make && make install && make clean

# samtools
RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 &&\
  tar -xf samtools-1.8.tar.bz2 && rm samtools-1.8.tar.bz2 &&\
  cd samtools-1.8 &&\
  ./configure --with-htslib=/opt/htslib-1.8 &&\
  make && make install && make clean

# STAR v2.5.3a
RUN wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz &&\
  tar -xf 2.5.3a.tar.gz && rm 2.5.3a.tar.gz &&\
  cd STAR-2.5.3a/source &&\
  make STAR && make STARlong &&\
  mv STAR* ../bin/Linux_x86_64/

ENV PATH /docker_home/opt/STAR-2.5.3a/bin/Linux_x86_64:$PATH

# RSEM v1.3.0
RUN wget --no-check-certificate https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz &&\
  tar -xf v1.3.0.tar.gz && rm v1.3.0.tar.gz &&\
  cd RSEM-1.3.0 &&\
  make

ENV PATH /docker_home/opt/RSEM-1.3.0:$PATH

# RNA-SeQC
RUN wget --no-check-certificate https://github.com/francois-a/rnaseqc/releases/download/v1.1.9/RNA-SeQC_1.1.9.zip &&\
  unzip RNA-SeQC_1.1.9.zip -d RNA-SeQC_1.1.9 && rm RNA-SeQC_1.1.9.zip

CMD tail -f "/docker_home/RunRSEM.sh"

