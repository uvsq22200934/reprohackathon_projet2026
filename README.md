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




# Projet 2 : réponse transcriptomique du cerveau de souris à l'infection par *Toxoplasma gondii* (Hu et al., Frontiers in Microbiology, 2020)

**Article scientifique :** [Hu et al. (2020) - Frontiers in Microbiology](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2020.570903/full)

## Introduction

Ce projet vise à reproduire une partie de l'analyse RNA-seq présentée dans l'article de Hu et al. (2020), qui étudie les modifications de l'expression des transcrits dans le cerveau de souris infectées par *Toxoplasma gondii*.

L'objectif est de comparer les souris infectées aux souris témoins à deux stades de l'infection : aigu (J11) et chronique (J33), afin d'identifier les transcrits différentiellement exprimés (DETs) et de reproduire le diagramme de Venn de la Figure 2B de l'article.

## Organisation du projet

```text
projet2/
│
├── config/
│   └── samples.tsv
│
├── containers/
│   ├── Dockerfile.wget
│   ├── Dockerfile.gffread
│   └── Dockerfile.R
│
├── data/
│   ├── metadata/
│   │   └── SraRunInfo.csv
│   └── reference/
│
├── scripts/
│   └── deseq2.R
│
├── results/              # Généré lors de l'exécution
│
├── Snakefile
├── run.sh
└── README.md
```

### Rôle des principaux fichiers et dossiers

- **`config/`** : contient les informations sur les échantillons (conditions expérimentales, jours d'infection et réplicats).
- **`containers/`** : contient les Dockerfiles permettant de construire les environnements logiciels nécessaires au workflow.
- **`data/`** : contient les métadonnées et les données de référence utilisées pour l'analyse.
- **`scripts/`** : contient les scripts d'analyse, notamment `deseq2.R` pour l'analyse différentielle.
- **`results/`** : accueillera les fichiers intermédiaires et les résultats générés par le pipeline.
- **`Snakefile`** : définit les différentes étapes du workflow et leurs dépendances à l'aide de Snakemake.
- **`run.sh`** : automatise la préparation des environnements Docker et le lancement du workflow.

## Étapes du workflow

### 1. Téléchargement des données RNA-seq

Téléchargement des reads RNA-seq depuis la base de données NCBI SRA à partir des identifiants des échantillons.

Les 12 échantillons sont répartis en quatre groupes :
- Infection aiguë (J11) : infectés et témoins.
- Infection chronique (J33) : infectés et témoins.

Chaque groupe comprend trois réplicats biologiques.

### 2. Nettoyage des reads (Trimmomatic)

Nettoyage des reads bruts avec Trimmomatic afin de retirer les adaptateurs et les séquences de mauvaise qualité.

### 3. Préparation du transcriptome de référence

Téléchargement du génome de souris GRCm38 et de son annotation Ensembl 100.

Utilisation de gffread pour extraire les séquences transcriptomiques à partir du génome et de son annotation.

### 4. Construction de l'index transcriptomique (Salmon)

Construction d'un index du transcriptome avec Salmon, nécessaire à la quantification des transcrits.

### 5. Quantification des transcrits (Salmon)

Quantification des reads nettoyés pour chaque échantillon avec Salmon.

Cette étape produit notamment un fichier `quant.sf` par échantillon, contenant les estimations d'abondance des transcrits (TPM et nombre estimé de reads).

### 6. Analyse différentielle (DESeq2)

Importation des fichiers `quant.sf` dans R avec tximport, puis analyse différentielle avec DESeq2.

Deux comparaisons sont réalisées :
- J11 : souris infectées vs souris témoins.
- J33 : souris infectées vs souris témoins.

L'objectif est d'identifier les transcrits différentiellement exprimés pour chaque stade de l'infection.

### 7. Reproduction du diagramme de Venn (Figure 2B)

Comparaison des listes de transcrits différentiellement exprimés obtenues à J11 et J33 afin d'identifier les transcrits communs et ceux spécifiques à chaque stade de l'infection.

Ces résultats permettront de reproduire le diagramme de Venn présenté dans la Figure 2B de l'article.

## Exécution du workflow

Le pipeline utilise Snakemake pour orchestrer les étapes et Docker pour gérer les environnements logiciels.

Depuis le dossier `projet2/` :

```bash
bash run.sh
```

Pour simuler les étapes du workflow :

```bash
snakemake -n --cores 2
```

Le mode simulation permet de vérifier le plan d'exécution sans lancer les analyses.

