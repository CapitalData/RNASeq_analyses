# adapted from https://github.com/CapitalData/gtex-pipeline/blob/master/rnaseq/Dockerfile

#RNAseq Dockerfile 

#use an official parent image

FROM ubuntu:16.04


#set the working directory

WORKDIR /docker_home


#copy the current directory contents into the container at /docker_home
COPY . /docker_home


RUN apt-get update && apt-get install -y software-properties-common && add-apt-repository -y ppa:openjdk-r/ppa && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        curl \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        openjdk-7-jdk \
        openjdk-8-jdk \
        #python3 \
        #python3-pip \
        r-base-core=3.2.3-4 \
        unzip \
        vim-common \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*



# R

#RUN wget https://github.com/downloads/PMBio/peer/R_peer_source_1.3.tgz && \
#    R CMD INSTALL R_peer_source_1.3.tgz && \
#    rm R_peer_source_1.3.tgz && \
#    echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
#    Rscript -e "install.packages(c('argparser’), dependencies=TRUE)" && \
#    Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("qvalue"); biocLite("sva"); biocLite("edgeR");'


    #below was original from gtex
    #Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("qvalue"); biocLite("sva"); biocLite("edgeR");'
    
    #here is what I added (last two packages) so we can follow the directions here: https://ycl6.gitbooks.io/rna-seq-data-analysis/install_r_libraries.html
    #RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("qvalue"); biocLite("sva"); biocLite("edgeR"); biocLite("cluster"); biocLite("Biobase");'

#found this in MaayanLab/Zika-RNAseq-Pipeline/blob/master/Dockerfile

RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("qvalue"); biocLite("sva"); biocLite("edgeR"); biocLite("cluster"); biocLite("Biobase");'



# standalone R math library
RUN wget https://cran.r-project.org/src/base/R-3/R-3.2.5.tar.gz && tar -xf R-3.2.5.tar.gz && rm R-3.2.5.tar.gz \
    && cd R-3.2.5 && ./configure --with-x=no && cd src/nmath/standalone && make
ENV RMATH /R-3.2.5/src

#I put this part - - is it right???

RUN R -e "install.packages(c('gapminder'), repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages(c('argparser’), repos = 'http://cran.us.r-project.org')"


#-----------------------------
# Pipeline components
#-----------------------------

RUN mkdir /docker_home/opt

# htslib
RUN cd /docker_home/opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2 && \
    tar -xf htslib-1.8.tar.bz2 && rm htslib-1.8.tar.bz2 && cd htslib-1.8 && \
    ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
    make && make install && make clean

# samtools
RUN cd /docker_home/opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 && \
    tar -xf samtools-1.8.tar.bz2 && rm samtools-1.8.tar.bz2 && cd samtools-1.8 && \
    ./configure --with-htslib=/opt/htslib-1.8 && make && make install && make clean


# STAR v2.5.3a
RUN cd /docker_home/opt && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz && \
    tar -xf 2.5.3a.tar.gz && rm 2.5.3a.tar.gz && \
    make STAR -C STAR-2.5.3a/source && make STARlong -C STAR-2.5.3a/source && \
    mv STAR-2.5.3a/source/STAR* STAR-2.5.3a/bin/Linux_x86_64/
ENV PATH /docker_home/opt/STAR-2.5.3a/bin/Linux_x86_64:$PATH

# RSEM v1.3.0
RUN cd /docker_home/opt && \
    wget --no-check-certificate https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz && \
    tar -xvf v1.3.0.tar.gz && rm v1.3.0.tar.gz && cd RSEM-1.3.0 && make
ENV PATH /docker_home/opt/RSEM-1.3.0:$PATH

# RNA-SeQC
RUN cd /docker_home/opt && \
    wget --no-check-certificate https://github.com/francois-a/rnaseqc/releases/download/v1.1.9/RNA-SeQC_1.1.9.zip && \
    unzip RNA-SeQC_1.1.9.zip -d RNA-SeQC_1.1.9 && rm RNA-SeQC_1.1.9.zip


#run runRSEM.sh

CMD /docker_home/RunRSEM.sh

