FROM ubuntu:22.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    sra-toolkit \
    ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /work