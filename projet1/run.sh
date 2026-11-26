#!/bin/bash
set -e

cd "$(dirname "$0")"

# Trouver Snakemake dans l'environnement actif
if command -v snakemake >/dev/null 2>&1; then
    SNAKEMAKE="$(command -v snakemake)"
elif [ -x "$HOME/snakemake-env/bin/snakemake" ]; then
    SNAKEMAKE="$HOME/snakemake-env/bin/snakemake"
else
    echo "ERREUR : Snakemake est introuvable."
    echo "Installez Snakemake ou activez votre environnement Conda."
    exit 1
fi

# Vérifier Docker
if ! command -v docker >/dev/null 2>&1; then
    echo "ERREUR : Docker n'est pas installé."
    exit 1
fi

if ! docker info >/dev/null 2>&1; then
    echo "ERREUR : Docker ne semble pas fonctionner."
    exit 1
fi

echo "Version de Snakemake :"
"$SNAKEMAKE" --version

# Construire les images Docker
echo "Construction des images Docker"
"$SNAKEMAKE" --snakefile Snakefile --cores 1 build_images

# Vérifier le workflow
echo "Vérification du workflow (dry-run)"
"$SNAKEMAKE" --snakefile Snakefile --dry-run --cores 4

# Lancer l'analyse
echo "Lancement de l'analyse"
"$SNAKEMAKE" --snakefile Snakefile --cores 4
