#!/bin/bash

#script pour lancer les installations necessaires pour quand vous lancez une VM (avoir snakemake etc)
set -e

sudo apt-get update
sudo apt-get install -y python3.10-venv

python3 -m venv "$HOME/snakemake-env"

"$HOME/snakemake-env/bin/python" -m pip install \
    "snakemake==7.32.4" \
    "pulp==2.7.0"

echo "Pour utiliser Snakemake :"
echo "source ~/snakemake-env/bin/activate"