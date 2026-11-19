# Analyse différentielle avec DESeq2
# Comparaison : Persister vs Control
# Projet 1 - Reprohackathon

library(DESeq2)

# 1. Définition des échantillons et conditions

samples <- c(
    "SRR10379721", "SRR10379722", "SRR10379723",
    "SRR10379724", "SRR10379725", "SRR10379726"
)

condition <- factor(
    c("Persister", "Persister", "Persister",
      "Control", "Control", "Control"),
    levels = c("Control", "Persister")
)

coldata <- data.frame(
    condition = condition,
    row.names = samples
)

# 2. Lecture des comptages featureCounts

files <- file.path(
    "results/counts",
    paste0(samples, ".txt")
)

if (!all(file.exists(files))) {
    stop("Un ou plusieurs fichiers de comptage sont manquants.")
}

count_tables <- lapply(files, function(file) {
    read.delim(file, comment.char = "#", check.names = FALSE)
})

gene_ids <- count_tables[[1]]$Geneid

for (i in seq_along(count_tables)) {
    if (!identical(gene_ids, count_tables[[i]]$Geneid)) {
        stop("Les identifiants des gènes diffèrent entre les échantillons.")
    }
}

count_matrix <- do.call(
    cbind,
    lapply(count_tables, function(x) {
        as.integer(x[[ncol(x)]])
    })
)

rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- samples

if (any(is.na(count_matrix))) {
    stop("La matrice de comptage contient des valeurs manquantes.")
}

# 3. Création du dataset DESeq2

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = ~ condition
)

# 4. Analyse différentielle

dds <- DESeq(dds)

res <- results(
    dds,
    contrast = c("condition", "Persister", "Control")
)

res <- res[order(res$padj), ]

# 5. Export des résultats

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

write.csv(
    as.data.frame(res),
    file = "results/deseq2/differential_expression.csv"
)

write.table(
    counts(dds, normalized = TRUE),
    file = "results/deseq2/normalized_counts.tsv",
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

saveRDS(dds, "results/deseq2/dds.rds")

# 6. Fonction commune pour les graphiques MA

plot_custom_MA <- function(data, title, output_file,
                           annotations = NULL,
                           aminoacyl_genes = NULL) {

    # Retirer les valeurs non exploitables
    valid <- !is.na(data$baseMean) &
             !is.na(data$log2FoldChange) &
             is.finite(data$baseMean) &
             is.finite(data$log2FoldChange)

    data <- data[valid, ]

    # Significativité statistique
    significant <- !is.na(data$padj) & data$padj < 0.05

    # Couleur des points
    point_colors <- ifelse(significant, "red", "grey70")

    # Position des gènes
    x <- log10(data$baseMean + 1)
    y <- data$log2FoldChange

    pdf(output_file, width = 9, height = 7)

    plot(
        x, y,
        pch = 16,
        col = point_colors,
        cex = 0.7,
        xlab = "log10(mean normalized counts + 1)",
        ylab = "log2 fold change (Persister / Control)",
        main = title,
        ylim = c(-5, 5)
    )

    abline(h = 0, lty = 2, col = "black")

    # Cercles noirs autour des aminoacyl-ARNt synthétases
    if (!is.null(aminoacyl_genes)) {
        selected <- rownames(data) %in% aminoacyl_genes

        if (any(selected)) {
            points(
                x[selected], y[selected],
                pch = 1,
                col = "black",
                cex = 1.4,
                lwd = 1.2
            )
        }
    }

    # Annotation des six gènes ciblés
    if (!is.null(annotations)) {
        selected <- match(rownames(data), annotations$gene_id)
        target <- !is.na(selected)

        if (any(target)) {
            text(
                x[target], y[target],
                labels = annotations$symbol[selected[target]],
                pos = 3,
                cex = 0.85,
                col = "black"
            )
        }
    }

    legend(
        "topright",
        legend = c(
            "Significant (padj < 0.05)",
            "Not significant",
            "Aminoacyl-tRNA synthetases"
        ),
        col = c("red", "grey70", "black"),
        pch = c(16, 16, 1),
        bty = "n",
        cex = 0.8
    )

    dev.off()
}

# 7. MA-plot de tous les gènes

plot_custom_MA(
    data = as.data.frame(res),
    title = "Persister vs Control - All genes",
    output_file = "results/deseq2/MA_plot_all_genes.pdf"
)

# 8. Lecture des gènes liés à la traduction

translation_file <- "data/kegg/translation_related_genes_gff.txt"
aminoacyl_file <- "data/kegg/aminoacyl_trna_genes.txt"
annotation_file <- "data/kegg/target_gene_annotations.tsv"

if (!file.exists(translation_file) ||
    !file.exists(aminoacyl_file) ||
    !file.exists(annotation_file)) {
    stop("Un ou plusieurs fichiers d'annotation sont manquants.")
}

translation_genes <- readLines(translation_file)
aminoacyl_genes <- readLines(aminoacyl_file)

annotations <- read.delim(
    annotation_file,
    header = TRUE,
    stringsAsFactors = FALSE
)

# Vérifier les correspondances avec les résultats DESeq2

translation_genes <- intersect(translation_genes, rownames(res))
aminoacyl_genes <- intersect(aminoacyl_genes, translation_genes)

if (length(translation_genes) == 0) {
    stop("Aucun gène de traduction ne correspond aux résultats DESeq2.")
}

# 9. MA-plot des gènes liés à la traduction

res_translation <- as.data.frame(res[translation_genes, ])

plot_custom_MA(
    data = res_translation,
    title = "Persister vs Control - Translation-related genes",
    output_file = "results/deseq2/MA_plot_translation_genes.pdf",
    annotations = annotations,
    aminoacyl_genes = aminoacyl_genes
)

# 10. Informations sur l'environnement R

writeLines(
    capture.output(sessionInfo()),
    "results/deseq2/sessionInfo.txt"
)
