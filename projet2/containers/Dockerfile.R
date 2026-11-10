# Image de base contenant déjà R et l'environnement Bioconductor.
FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Installation des packages nécessaires à l'analyse différentielle :
# DESeq2 pour l'analyse d'expression, tximport pour importer les quantifications
# Salmon et VennDiagram pour générer les diagrammes de Venn.
RUN R -e "BiocManager::install(c('DESeq2', 'tximport'), ask=FALSE, update=FALSE)" \
    && R -e "install.packages('VennDiagram', repos='https://cloud.r-project.org')"

# Définit /project comme répertoire de travail à l'intérieur du conteneur.
WORKDIR /project
