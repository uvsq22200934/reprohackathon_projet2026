#!/bin/bash
set -eo pipefail

cd "$(dirname "$0")"

echo "========================================"
echo " RNA-seq reproducibility workflow"
echo "========================================"

# Installation automatique si nécessaire
bash setup.sh

# Activation de l'environnement Damona
source "$HOME/.config/damona/damona.sh"
if damona activate reprohackathon_p1; then
    :
else
    echo "ERROR: Damona activation failed."
    exit 1
fi

# Activation de l'environnement Python
source "$HOME/snakemake-env/bin/activate"

echo ""
echo "Checking tool versions..."

fasterq-dump --version
samtools --version | head -n 2

echo ""
echo "Checking Docker images..."

docker image inspect bowtie-0127 >/dev/null
docker image inspect cutadapt-1.11 >/dev/null
docker image inspect featurecounts-146p3 >/dev/null
docker image inspect deseq2:1.16.1 >/dev/null
docker image inspect workflow-utils >/dev/null

echo "All Docker images are available."

echo ""
echo "Checking Snakemake workflow..."
snakemake --snakefile Snakefile --dry-run

echo ""
echo "Starting RNA-seq analysis..."
snakemake --snakefile Snakefile --cores 4 --rerun-incomplete

echo ""
echo "RNA-seq workflow completed successfully."
