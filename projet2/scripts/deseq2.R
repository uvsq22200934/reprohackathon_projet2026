# ============================================================
# PROJET 2 - Analyse différentielle DESeq2
# Hu et al. (2020) - Infection aiguë (J11) et chronique (J33)
#
# Comparaisons :
#   In versus WT à J11
#   In versus WT à J33
#
# Deux critères DET sont étudiés :
#   1. p-value brute < 0.05 et |log2FC| >= 1.5
#   2. p-value ajustée < 0.05 et |log2FC| >= 1.5
# ============================================================

library(tximport)
library(DESeq2)

# ------------------------------------------------------------
# 1. Lecture des métadonnées
# ------------------------------------------------------------

metadata <- read.delim(
  "config/samples.tsv",
  stringsAsFactors = FALSE
)

required_columns <- c("sample", "condition", "day", "replicate")

if (!all(required_columns %in% colnames(metadata))) {
  stop("Colonnes manquantes dans config/samples.tsv")
}

if (nrow(metadata) != 12) {
  stop("Le fichier samples.tsv doit contenir 12 échantillons.")
}

# ------------------------------------------------------------
# 2. Préparation des fichiers Salmon
# ------------------------------------------------------------

files <- file.path(
  "results/salmon",
  metadata$sample,
  "quant.sf"
)

names(files) <- metadata$sample

if (!all(file.exists(files))) {
  stop("Certains fichiers quant.sf sont absents.")
}

dir.create(
  "results/deseq2",
  recursive = TRUE,
  showWarnings = FALSE
)

summary_results <- data.frame()

# ------------------------------------------------------------
# 3. Analyse séparée à J11 et J33
# ------------------------------------------------------------

for (day_value in c(11, 33)) {

  cat("\nAnalyse du jour", day_value, "\n")

  metadata_day <- metadata[
    metadata$day == day_value,
  ]

  if (nrow(metadata_day) != 6) {
    stop(paste("Nombre incorrect d'échantillons à J", day_value))
  }

  metadata_day$condition <- factor(
    metadata_day$condition,
    levels = c("WT", "In")
  )

  if (any(is.na(metadata_day$condition))) {
    stop("Valeur de condition inconnue dans samples.tsv")
  }

  rownames(metadata_day) <- metadata_day$sample

  files_day <- files[metadata_day$sample]

  # Import des quantifications au niveau transcriptique
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

  # Analyse différentielle
  dds <- DESeq(dds)

  # Comparaison In versus WT
  res <- results(
    dds,
    contrast = c("condition", "In", "WT"),
    alpha = 0.05
  )

  # Conversion en tableau
  res_df <- as.data.frame(res)
  res_df$transcript_id <- rownames(res_df)

  res_df <- res_df[
    order(res_df$pvalue, na.last = TRUE),
  ]

  # Sauvegarde de tous les résultats
  write.table(
    res_df,
    file = paste0(
      "results/deseq2/day", day_value, "_all.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # 4. Définition des DET
  # ----------------------------------------------------------

  # Critère 1 : p-value brute
  dets_raw <- res_df[
    !is.na(res_df$pvalue) &
    res_df$pvalue < 0.05 &
    abs(res_df$log2FoldChange) >= 1.5,
  ]

  # Critère 2 : p-value ajustée
  dets_fdr <- res_df[
    !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 1.5,
  ]

  # Sauvegarde des listes complètes
  write.table(
    dets_raw,
    file = paste0(
      "results/deseq2/day", day_value, "_DETs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    dets_fdr,
    file = paste0(
      "results/deseq2/day", day_value, "_DETs_FDR.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # 5. Séparation surexprimés / sous-exprimés
  # ----------------------------------------------------------

  if (day_value == 11) {
    prefix <- "acute"
  } else {
    prefix <- "chronic"
  }

  # Listes basées sur la p-value brute
  up_raw <- dets_raw[
    dets_raw$log2FoldChange > 0,
  ]

  down_raw <- dets_raw[
    dets_raw$log2FoldChange < 0,
  ]

  # Listes basées sur la p-value ajustée
  up_fdr <- dets_fdr[
    dets_fdr$log2FoldChange > 0,
  ]

  down_fdr <- dets_fdr[
    dets_fdr$log2FoldChange < 0,
  ]

  write.table(
    up_raw,
    paste0("results/deseq2/", prefix, "_upregulated.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    down_raw,
    paste0("results/deseq2/", prefix, "_downregulated.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    up_fdr,
    paste0("results/deseq2/", prefix, "_upregulated_FDR.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    down_fdr,
    paste0("results/deseq2/", prefix, "_downregulated_FDR.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # 6. Résumé des effectifs
  # ----------------------------------------------------------

  summary_results <- rbind(
    summary_results,
    data.frame(
      day = day_value,
      raw_total = nrow(dets_raw),
      raw_up = nrow(up_raw),
      raw_down = nrow(down_raw),
      FDR_total = nrow(dets_fdr),
      FDR_up = nrow(up_fdr),
      FDR_down = nrow(down_fdr)
    )
  )

  cat("DET avec p-value brute :", nrow(dets_raw), "\n")
  cat("  Surexprimés :", nrow(up_raw), "\n")
  cat("  Sous-exprimés :", nrow(down_raw), "\n")

  cat("DET avec p-value ajustée :", nrow(dets_fdr), "\n")
  cat("  Surexprimés :", nrow(up_fdr), "\n")
  cat("  Sous-exprimés :", nrow(down_fdr), "\n")
}

# ------------------------------------------------------------
# 7. Sauvegarde du résumé général
# ------------------------------------------------------------

write.table(
  summary_results,
  "results/deseq2/deseq2_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nAnalyse DESeq2 terminée.\n")
EOF
