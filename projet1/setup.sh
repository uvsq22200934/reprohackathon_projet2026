#!/bin/bash
set -eo pipefail

cd "$(dirname "$0")"

echo "========================================"
echo " Setting up RNA-seq reproducibility"
echo "========================================"

# Initialisation de Damona
if [ ! -f "$HOME/.config/damona/damona.sh" ]; then
    echo "ERROR: Damona configuration not found."
    exit 1
fi

source "$HOME/.config/damona/damona.sh"

: "${DAMONA_EXE:=}"

# Vérification des prérequis système
for cmd in docker apptainer python3 damona; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd is not installed."
        exit 1
    fi
done

# Vérifier que Docker fonctionne
docker info >/dev/null

# Installer les dépendances Python si nécessaire
if ! python3 -m venv --help >/dev/null 2>&1; then
    echo "Installing python3-venv..."
    sudo apt-get update
    sudo apt-get install -y python3-venv
fi

# Environnement Snakemake
SNAKEMAKE_ENV="$HOME/snakemake-env"

if [ ! -x "$SNAKEMAKE_ENV/bin/snakemake" ]; then
    echo "Creating Snakemake environment..."
    python3 -m venv "$SNAKEMAKE_ENV"
fi

if ! "$SNAKEMAKE_ENV/bin/snakemake" --version 2>/dev/null | grep -qx "7.32.4" \
   || ! "$SNAKEMAKE_ENV/bin/python" -c \
       'import pulp; assert pulp.__version__ == "2.7.0"' 2>/dev/null; then

    echo "Installing Snakemake dependencies..."
    "$SNAKEMAKE_ENV/bin/python" -m pip install \
        "snakemake==7.32.4" \
        "pulp==2.7.0"
fi

# Initialisation de Damona
DAMONA_ENV="reprohackathon_p1"

# Encadrer les commandes Damona pour eviter les problemes avec set -e
run_damona() {
    if damona "$@"; then
        return 0
    else
        return $?
    fi
}


if ! "$DAMONA_EXE_INTERN" info "$DAMONA_ENV" >/dev/null 2>&1; then
    echo "Creating Damona environment..."
    run_damona create "$DAMONA_ENV"
fi

run_damona activate "$DAMONA_ENV"

# Installation des outils Damona si absents
if ! "$DAMONA_EXE_INTERN" info "$DAMONA_ENV" 2>/dev/null | grep "sratoolkit_3.1.0.img" >/dev/null; then
    echo "Installing SRA Toolkit 3.1.0..."
    run_damona install sratoolkit:3.1.0
fi

if ! "$DAMONA_EXE_INTERN" info "$DAMONA_ENV" 2>/dev/null | grep "samtools_1.15.0.img" >/dev/null; then
    echo "Installing Samtools 1.15.0..."
    run_damona install samtools:1.15.0
fi

# Construction des images Docker uniquement si elles sont absentes
build_if_missing() {
    IMAGE="$1"
    DOCKERFILE="$2"

    if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
        echo "Building Docker image: $IMAGE"
        docker build -t "$IMAGE" -f "$DOCKERFILE" .
    else
        echo "Docker image already available: $IMAGE"
    fi
}

build_if_missing bowtie-0127 containers/Dockerfile.bowtie0127
build_if_missing featurecounts-146p3 containers/Dockerfile.featurecounts
build_if_missing workflow-utils containers/Dockerfile.utils
build_if_missing deseq2:1.16.1 containers/Dockerfile.deseq2
build_if_missing cutadapt-1.11 containers/Dockerfile.cutadapt

echo "========================================"
echo " Setup completed successfully"
echo "========================================"
