FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    unzip \
    zlib1g \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

RUN wget -O bowtie.zip \
    "https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-src.zip/download" \
    && unzip bowtie.zip \
    && cd bowtie-0.12.7 \
    && make \
    && cp bowtie bowtie-build bowtie-inspect /usr/local/bin/ \
    && cd /opt \
    && rm -rf bowtie.zip bowtie-0.12.7

CMD ["bowtie", "--version"]
