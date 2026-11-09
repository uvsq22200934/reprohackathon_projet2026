# Import des packages nécessaires
library(tximport)
library(DESeq2)

# Lecture des informations des échantillons
metadata <- read.delim(
  "config/samples.tsv",
  stringsAsFactors = FALSE
)

# Création des chemins vers les fichiers quant.sf produits par Salmon
files <- file.path(
  "results/salmon",
  metadata$sample,
  "quant.sf"
)

names(files) <- metadata$sample

# Vérification de la présence des fichiers
if (!all(file.exists(files))) {
  stop("Certains fichiers quant.sf sont absents. Vérifier la quantification Salmon.")
}

# Création du dossier de sortie
dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

# Analyse séparée des deux temps d'infection
for (day_value in c(11, 33)) {

  # Sélection des échantillons du jour étudié
  metadata_day <- metadata[metadata$day == day_value, ]
  files_day <- files[metadata_day$sample]

  # La condition WT sert de référence pour la comparaison
  metadata_day$condition <- factor(
    metadata_day$condition,
    levels = c("WT", "In")
  )

  rownames(metadata_day) <- metadata_day$sample

  # Import des quantifications Salmon au niveau des transcrits
  txi <- tximport(
    files_day,
    type = "salmon",
    txOut = TRUE
  )

  # Création de l'objet DESeq2
  dds <- DESeqDataSetFromTximport(
    txi,
    colData = metadata_day,
    design = ~ condition
  )

  # Filtrage des transcrits très faiblement exprimés
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]

  # Analyse différentielle
  dds <- DESeq(dds)

  # Comparaison In versus WT
  res <- results(
    dds,
    contrast = c("condition", "In", "WT"),
    alpha = 0.05
  )

  # Conversion en tableau et tri par p-value
  res_df <- as.data.frame(res)
  res_df <- res_df[order(res_df$pvalue), ]

  # Sauvegarde de tous les résultats
  write.table(
    res_df,
    file = paste0("results/deseq2/day", day_value, "_all.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  # Sélection des transcrits différentiellement exprimés
  dets <- res_df[
    !is.na(res_df$pvalue) &
    res_df$pvalue < 0.05 &
    abs(res_df$log2FoldChange) >= 1.5,
  ]

  # Sauvegarde des DET
  write.table(
    dets,
    file = paste0("results/deseq2/day", day_value, "_DETs.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  cat(
    "Jour", day_value, ":",
    nrow(dets), "transcrits différentiellement exprimés\n"
  )
}
