FROM ubuntu:18.04

RUN apt-get -y update \
	&& apt-get -y install make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree \
    liblzma-dev wget curl less pkg-config git libncurses5-dev \
    bedtools python3 \
	&& apt-get clean

RUN ln -s /usr/bin/python3 /usr/bin/python

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 \
    && tar xvf samtools-1.20.tar.bz2 \
    && rm samtools-1.20.tar.bz2 \
    && cd samtools-1.20 \
    && ./configure \
    && make \
    && make install

# install latest bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 \
    && tar xvf bcftools-1.20.tar.bz2 \
    && rm bcftools-1.20.tar.bz2 \
    && cd bcftools-1.20 \
    && make \
    && make install

# install latest htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 \
    && tar xvf htslib-1.20.tar.bz2 \
    && rm htslib-1.20.tar.bz2 \
    && cd htslib-1.20 \
    && ./configure \
    && make \
    && make install

# install seqkit
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.9.0/seqkit_linux_amd64.tar.gz \
    && tar xvf seqkit_linux_amd64.tar.gz \
    && rm seqkit_linux_amd64.tar.gz \
    && mv seqkit /usr/bin/
