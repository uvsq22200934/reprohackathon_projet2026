# reprohackathon_projet2026


#Projet 1

Peyrusson et al. (2020) — Intracellular Staphylococcus aureus persisters upon antibiotic exposure, Nature Communications.

L'étude s'intéresse au Staphylocoque doré dont les bactéries persistantes sont capables de survivre à un traitement antibiotique à l'intérieur des cellules hôtes.

Les auteurs réalisent une analyse RNA-seq pour comparer l'expression des gènes entre des bactéries persistantes exposées à un antibiotique, et des bactéries témoins.

L'objectif est de réaliser cette analyse. Il faut construire les étapes bioinformatics depuis les données brutes jusqu'à l'analyse d'expression.

Dans l'analyse, on utilise 6 échantillons RNA-seq : 
SRR10379721 (bactérie persistante, réplicat 1)
SRR10379722(bactérie persistante, réplicat 2)
SRR10379723(bactérie persistante, réplicat 3)
SRR10379724(bactérie contrôle, réplicat 1)
SRR10379725(bactérie contrôle, réplicat 2)
SRR10379726(bactérie contrôle, réplicat 3)


Les échantillons de bactéries persistantes ont survécu a l'exposition de l'oxaciline pendant 1 jour 
Les échantillons contrôlés correspondent à des bactéries témoins non exposées

Le génome de référence utilisé est celui de Staphylococcus aureus NCTC 8325, accession CP000253.1.



ORGANISATION DU PROJET

Snakefile : definit les étapes de l'analyse et les dépendances
setup.sh : prepare l'environnement pour l'exécution
run.sh : lancer l'analyse
containers/ : contient les containers pour construire les images Docker.
data/ : données d'entrée et genome de ref
resultats/ : résultats des analyses


Logiciels utilisés : 
Python 3.10.12 
Snakemake 7.32.4 
Docker 29.4.1 
SRA tools 3.0.3 Téléchargement des données RNA-seq
Cutadapt 1.11 Nettoyage des lectures
Bowtie 1.3.1 alignement sur le génome de ref
SAMtools 1.13 Tri et indexation des alignements
featureCounts 1.4.6-p3 comptage des lectures

deux images sont construites à partir de Dockerfiles :
bowtie-0127 pour bowtie
featurecounts-146p3 pour featurecounts

Les autres outils sont exécutes à partir d'images de Biocontainers



Le workflow :

1) téléchargement des données
   rules : download_genome, download_annotation

2) preparation des annotations(revoir)
  rule : convert_gff_to_saf

3) indexation du genome
   rule : index_genome
 Création d'un index du génome avec bowtie-build.

L'index permet à Bowtie de rechercher efficacement la position des lectures RNA-seq sur le génome de référence.

4) Téléchargement des données RNA-seq
   rule : download_fastq
   environ 9.2 Go compréssés

5) nettoyage des lectures
   rule : trim_reads
   Nettoyage des fichiers faste avec cutadapt

   Paramètres utilisés :

-q 20 : retrait des bases de faible qualité aux extrémités
-m 25 : conservation des lectures d'au moins 25 nucléotides
-a AGATCGGAAGAGC: recherche et retrait de la séquence d'adaptateur indiquée

>> limiter les erreurs lors de l'alignement

6) alignement sur le genome avec Bowtie
   rule :align_reads

7) Comptage des lectures
   rule : count_reads

8) regle finale
  rule : all

  snakemake identifie automatiquement les étapes nécessaires, leur ordre d'exécution et les fichiers qui ont deja été produits.



  (j'ai pas fini mais en gros)







