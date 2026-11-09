# Image de base contenant déjà R et l'environnement Bioconductor.
FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Installation des packages R nécessaires à l'analyse différentielle :
# DESeq2 pour comparer l'expression des transcrits et tximport pour importer
# les fichiers de quantification produits par Salmon (quant.sf).
RUN R -e "BiocManager::install(c('DESeq2', 'tximport'), ask=FALSE, update=FALSE)"

# Définit /project comme répertoire de travail à l'intérieur du conteneur.
WORKDIR /project
