#!/bin/bash

#script pour lancer les installations necessaires pour quand vous lancez une VM (avoir snakemake etc)
set -e

sudo apt-get update
sudo apt-get install -y python3.10-venv curl

python3 -m venv "$HOME/snakemake-env"

"$HOME/snakemake-env/bin/python" -m pip install \
    "snakemake==7.32.4" \
    "pulp==2.7.0"



cd "$(dirname "$0")"

#Construire les images Docker nécessaires
docker build -t bowtie-0127 \
    -f containers/Dockerfile.bowtie0127 .

docker build -t featurecounts-146p3 \
    -f containers/Dockerfile.featurecounts .

docker build -t workflow-utils \
    -f containers/Dockerfile.utils .

#image contenant R et DESeq2
docker build -t deseq2:1.16.1 \
    -f containers/Dockerfile.deseq2 .


echo "Pour utiliser Snakemake :"
echo "source ~/snakemake-env/bin/activate"