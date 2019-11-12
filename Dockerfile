# SV-STAT (standalone) -- container with configured environment and all prerequisites
# 
# VERSION       0.1
#
#
# No copyright, license, warranty, or support.
#
# Build: 
#  docker build --rm -t cfdavis/docker-sv-stat .
#
# Run: 
#  docker run --rm -v /path/to/reference:/work/ref -v /path/to/svstat/repo:/work/sv-stat -it cfdavis/docker-sv-stat /work/sv-stat/src/svstat.sh -m /work/sv-stat/test/bam/hs1011.bam /work/sv-stat/test/metadata/ABori.txt

FROM centos:7

WORKDIR /work

#### BASE OS ####
RUN yum update -y
RUN yum install -y epel-release

RUN yum install -y git
RUN yum install -y unzip
RUN yum install -y bzip2
RUN yum install -y wget
RUN yum install -y make
RUN yum install -y zlib-devel
RUN yum install -y ncurses-devel
RUN yum install -y gcc-c++
RUN yum install -y gcc

#### TOOLS ####

# samtools
RUN wget https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
RUN tar -xvjf samtools-0.1.19.tar.bz2 
WORKDIR /work/samtools-0.1.19
RUN make CFLAGS=-fPIC
RUN ln -s /work/samtools-0.1.19/samtools /usr/local/bin/samtools
WORKDIR /work

# cdbfasta
RUN git clone https://github.com/gpertea/cdbfasta.git
WORKDIR /work/cdbfasta
RUN make
RUN ln -s /work/cdbfasta/cdbfasta /usr/local/bin/cdbfasta
RUN ln -s /work/cdbfasta/cdbyank /usr/local/bin/cdbyank
WORKDIR /work
# CDBFASTA_BIN_DIR=/work/cdbfasta

# bwa
RUN wget https://github.com/lh3/bwa/archive/0.5.9.tar.gz
RUN tar -xvzf 0.5.9.tar.gz
WORKDIR /work/bwa-0.5.9
RUN make
RUN ln -s /work/bwa-0.5.9/bwa /usr/local/bin/bwa
WORKDIR /work
# BWA_BIN_DIR=/work/bwa-0.5.9/bwa

# java
RUN yum install -y java-1.8.0-openjdk

# picard
RUN wget https://sourceforge.net/projects/picard/files/picard-tools/1.40/picard-tools-1.40.zip
RUN unzip picard-tools-1.40.zip
# PICARD_JAR_DIR=/work/picard-tools-1.40

# bedtools
RUN wget https://github.com/arq5x/bedtools/archive/Version-2.12.0.tar.gz
RUN tar -xvzf Version-2.12.0.tar.gz
WORKDIR /work/bedtools-Version-2.12.0
RUN make
WORKDIR /work
# BEDTOOLS_BIN_DIR=/work/bedtools-Version-2.12.0/bin

# bioperl
RUN yum install -y perl-devel
RUN yum install -y perl-App-cpanminus
RUN cpanm Test::More
RUN cpanm Env
RUN cpanm File::Which
RUN cpanm Fatal

RUN yum install -y libxml2-devel

RUN cpanm TAP::Harness
RUN cpanm Archive::Tar
RUN cpanm Pod::Readme
RUN cpanm inc::latest
RUN cpanm MLDBM

RUN yum install -y gd-devel

RUN cpanm JSON
RUN cpanm Tk

RUN yum install -y libxslt-devel

RUN cpanm XML::LibXSLT
RUN cpanm Graph::Directed
RUN cpanm GD
RUN cpanm Bio::Root::Version

ENV SAMTOOLS=/work/samtools-0.1.19
RUN cpanm Bio::DB::Sam

#RUN echo ${SAMTOOLS}


#### CONFIGURE ENVIRONMENT ####

ENV SVSTAT_SRC_DIR=/work/sv-stat/src
ENV CDBFASTA_BIN_DIR=/work/cdbfasta
ENV BWA_BIN_DIR=/work/bwa-0.5.9
ENV PICARD_JAR_DIR=/work/picard-tools-1.40
ENV BEDTOOLS_BIN_DIR=/work/bedtools-Version-2.12.0/bin
ENV REF_DIR=/work/ref
