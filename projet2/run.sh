
#!/usr/bin/env bash

# ============================================================
# PROJET 2 - RNA-seq Toxoplasma gondii
#
# Script de lancement du workflow Snakemake.
#
# Utilisation :
#   ./run.sh --dry-run    Vérifier le workflow sans l'exécuter
#   ./run.sh              Exécuter le workflow complet
#
# ============================================================

set -euo pipefail

# ------------------------------------------------------------
# SE PLACER DANS LE REPERTOIRE DU PROJET
# ------------------------------------------------------------

# Permet dexécuter le script depuis nimporte quel dossier.
cd "$(dirname "$0")"

# Nombre maximal de cœurs utilisés par Snakemake.
CORES=2

# ------------------------------------------------------------
# VERIFIER LES DEPENDANCES
# ------------------------------------------------------------

echo "Vérification des dépendances..."

for command in docker snakemake; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "Erreur : $command n'est pas installé ou accessible."
        exit 1
    fi
done

# Vérifier que Docker est opérationnel.
docker info >/dev/null 2>&1 || {
    echo "Erreur : Docker n'est pas accessible."
    exit 1
}

# Vérifier la présence des fichiers nécessaires.
for file in Snakefile config/samples.tsv \
    containers/Dockerfile.wget \
    containers/Dockerfile.gffread \
    containers/Dockerfile.R; do

    if [[ ! -f "$file" ]]; then
        echo "Erreur : fichier manquant : $file"
        exit 1
    fi
done

echo "Dépendances vérifiées."

# ------------------------------------------------------------
# CONSTRUCTION DES IMAGES DOCKER
# ------------------------------------------------------------

# Construire les images personnalisées à partir des Dockerfiles.
# Les images seront reproductibles sur les machines des collègues.

echo "Construction de l'image wget..."
docker build \
    -t reprohackathon-wget:1.0 \
    -f containers/Dockerfile.wget .

echo "Construction de l'image gffread..."
docker build \
    -t reprohackathon-gffread:1.0 \
    -f containers/Dockerfile.gffread .

# Construire l'image contenant R, DESeq2 et tximport.
echo "Construction de l'image R..."
docker build \
    -t reprohackathon-r:1.0 \
    -f containers/Dockerfile.R .

# Télécharger les images nécessaires aux autres étapes.
echo "Récupération des images bioinformatiques..."

docker pull quay.io/biocontainers/sra-tools:3.0.3--h87f3376_0

docker pull quay.io/biocontainers/trimmomatic:0.39--1

docker pull quay.io/biocontainers/salmon:1.2.0--hf69c8f4_0

# ------------------------------------------------------------
# CREATION DES REPERTOIRES
# ------------------------------------------------------------

mkdir -p \
    data/fastq \
    data/tmp \
    data/reference \
    results/trimmed \
    results/salmon \
    logs

# ------------------------------------------------------------
# EXECUTION DU WORKFLOW
# ------------------------------------------------------------

if [[ "${1:-}" == "--dry-run" ]]; then

    echo "Lancement du dry-run Snakemake..."

    snakemake -n --cores "$CORES"

elif [[ $# -eq 0 ]]; then

    echo "Lancement du workflow Snakemake..."

    snakemake --cores "$CORES" --printshellcmds

else

    echo "Argument inconnu : $1"
    echo "Utilisation : ./run.sh [--dry-run]"
    exit 1

fi

echo "Fin du script."