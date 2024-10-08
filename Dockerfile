# syntax=docker/dockerfile:1

# INSTALLED: 

FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y upgrade && \
apt-get install -y python3-pip && \
apt-get install -y wget tabix unzip && \
cd / && mkdir -p tools

# static files

RUN cd /tools && \
wget ftp://ftp.ebi.ac.uk/1000g/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# updates / packages
ENV LANGUAGE=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# long downloads are up, worth updating again
RUN apt-get update 

RUN apt-get install -y build-essential gfortran xorg-dev libpcre3-dev \
        libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev git fort77 libreadline-dev \
        cmake curl libboost-all-dev libgd-dev default-jre nano libncurses5 bc locales bsdmainutils gawk && \
        locale-gen en_US.UTF-8

# necessary toolkits

# anaconda

RUN cd /tools && \
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh && \
bash Anaconda3-2021.05-Linux-x86_64.sh -b -p /tools/anaconda && \
rm Anaconda3-2021.05-Linux-x86_64.sh


RUN apt-get install pip
SHELL ["/bin/bash", "-c"]

# samtools

RUN cd /tools && \
    wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && \
    tar xjf samtools-1.12.tar.bz2 && \
    rm samtools-1.12.tar.bz2 && \
    cd samtools-1.12 && \
    ./configure --prefix $(pwd) && \
    make

RUN apt-get install minimap2
RUN apt-get install -y bcftools
	
RUN cd /tools && \ 
wget https://netcologne.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2 && \ 
tar xjf bwa-0.7.17.tar.bz2 && \ 
rm bwa-0.7.17.tar.bz2 && \
cd bwa-0.7.17 && \ 
make

# bwakit

RUN cd /tools && \
wget https://netcologne.dl.sourceforge.net/project/bio-bwa/bwakit/bwakit-0.7.15_x64-linux.tar.bz2 && \ 
tar xjf bwakit-0.7.15_x64-linux.tar.bz2 && \ 
rm bwakit-0.7.15_x64-linux.tar.bz2 && \
cp bwa.kit/resource-GRCh38/hs38DH.fa.alt .


ENV PATH=${PATH}:/tools/samtools-1.12:/tools/bwa-0.7.17:/tools/anaconda/bin:/tools/minimap2-master/minimap2

RUN conda init bash

RUN conda create --name workspace python=3.8

RUN conda run -n workspace pip install pysam 
RUN conda run -n workspace pip install truvari
RUN conda run -n workspace pip install --upgrade pip setuptools==57.5.0
RUN conda run -n workspace pip install PyVCF3
RUN conda run -n workspace pip install opencv-python
RUN conda run -n workspace pip install tensorflow==2.12.0
RUN conda run -n workspace pip install tqdm

# bedtools2
RUN apt-get install bedtools

# SNIFFLES
RUN pip install sniffles==2.0.7

# CUTESV
RUN pip install cuteSV==2.0.2

# SVIM
RUN pip install svim==2.0.0

# DYSGU
RUN pip install numpy dysgu==1.3.16

# NANOVAR
RUN conda run -n workspace pip install nanovar==1.4.1

# install makeblastdb and windowmasker
RUN cd /tools && \
	wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/ncbi-blast-2.3.0+-x64-linux.tar.gz && \
	tar zxf ncbi-blast-2.3.0+-x64-linux.tar.gz && \
	cp ncbi-blast-2.3.0+/bin/makeblastdb ~/bin && cp ncbi-blast-2.3.0+/bin/windowmasker ~/bin

# install hs-blastn
RUN cd /tools && \
	git clone https://github.com/chenying2016/queries.git && \
	sed -i 's/isnan/std::isnan/g' /tools/queries/hs-blastn-src/v0.0.5/sources/utility.h && \
	cd queries/hs-blastn-src/v0.0.5 && \
	make && \
	cp hs-blastn ~/bin

RUN curl -s https://get.nextflow.io | bash && \
mv nextflow /bin

RUN conda install -y -c bioconda pbmm2
RUN conda install -y -c bioconda pbsv
RUN apt install tabix

RUN cd /tools && \
        git clone https://github.com/SFGLab/ConsensusSV-ONT-pipeline.git

WORKDIR /tools
