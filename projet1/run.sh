#!/bin/bash

#test verif du pipeline Snakemake
echo "Vérification du workflow (dry-run)"
snakemake --snakefile Snakefile --dry-run 

#Lancement de l'analyse
echo "Lancement du pipeline Snakemake"
snakemake --snakefile Snakefile --cores 4