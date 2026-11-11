#!/bin/bash
set -e #arreter le script si erreur

cd "$(dirname "$0")"

#activer env python contenant snakemake
source "$HOME/snakemake-env/bin/activate"

#test verif du pipeline Snakemake
echo "Vérification du workflow (dry-run)"
snakemake --snakefile Snakefile --dry-run 

#Lancement de l'analyse
echo "Lancement du pipeline Snakemake"
snakemake --snakefile Snakefile --cores 4