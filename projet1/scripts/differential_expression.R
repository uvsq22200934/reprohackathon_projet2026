
#Production des figures avec R

# Figure 3 supp et figure 3c :
# Expression différentielle avec DESeq2 : Persister vs Control

library(DESeq2)


echantillons <- c("SRR10379721", "SRR10379722", "SRR10379723",
  "SRR10379724", "SRR10379725", "SRR10379726") #echantillons

condition <- factor(c("Persister", "Persister", "Persister",
    "Control", "Control", "Control"),
    levels = c("Control", "Persister")) #on leur attribue leurs conditions

tab_data <- data.frame(condition = condition, row.names = echantillons)


#lectures des fichiers de comptage :

files <- file.path("results/counts", paste0(echantillons, ".txt"))

if (!all(file.exists(files))) {
  stop("fichiers manquants")
}

count_tables <- lapply(files, function(file) {
  read.delim(file, comment.char = "#", check.names = FALSE)
})


gene_ids <- count_tables[[1]]$Geneid #verif que tous les echs ont le mm identifiant de gene

for (i in seq_along(count_tables)) {
  if (!identical(gene_ids, count_tables[[i]]$Geneid)) {
    stop("id de gene differe selon l'ech")
  }
}


comptage_matrix <- do.call(
  cbind, lapply(count_tables, function(x) {
    as.integer(x[[ncol(x)]])
  })
)

rownames(comptage_matrix) <- gene_ids
colnames(comptage_matrix) <- echantillons


dds <- DESeqDataSetFromMatrix( #creation dataset DESeq2
  countData = comptage_matrix,
  colData = tab_data,
  design = ~ condition
)


#expression differentielle

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Persister", "Control"))

#tri selon la pvalue
res <- res[order(res$padj), ]


#exporter les donnees:

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

write.csv(as.data.frame(res),
  file = "results/deseq2/differential_expression.csv")

write.table(counts(dds, normalized = TRUE),
  file = "results/deseq2/normalized_counts.tsv",
  sep = "\t", quote = FALSE, col.names = NA)

saveRDS(dds, "results/deseq2/dds.rds")


#MA-plot global : Supplementary Figure 3

significatifs <- !is.na(res$padj) & res$padj < 0.05

pdf("results/deseq2/MA_plot.pdf", width = 8, height = 6)

plot(log10(res$baseMean + 1), res$log2FoldChange,
  pch = 16,
  col = ifelse(significatifs, "red", "grey"),
  xlab = "log10(mean expression)",
  ylab = "log2 fold change",
  main = "Persister vs Control",
  ylim = c(-5, 5))

abline(h = 0, lty = 2)

dev.off()


#MA-plot des genes impliques dans la translation : Figure 3c

translation_genes <- readLines("resources/translation_genes_kegg.txt")

res_translation <- res[rownames(res) %in% translation_genes, ]

if (nrow(res_translation) == 0) {
  stop("Aucun gene KEGG Translation ne correspond aux identifiants des comptages")
}

significatifs_translation <- !is.na(res_translation$padj) &
  res_translation$padj < 0.05

pdf("results/deseq2/MA_plot_translation.pdf", width = 8, height = 6)

plot(log10(res_translation$baseMean + 1),
  res_translation$log2FoldChange,
  pch = 16,
  col = ifelse(significatifs_translation, "red", "grey"),
  xlab = "log10(mean expression)",
  ylab = "log2 fold change",
  main = "Translation genes - Persister vs Control",
  ylim = c(-5, 5))

abline(h = 0, lty = 2)

#annotation des six genes mentionnes dans les diapositives

noms_genes <- c(
  "SAOUHSC_01236" = "frr",
  "SAOUHSC_02489" = "infA",
  "SAOUHSC_01234" = "tsf",
  "SAOUHSC_01786" = "infC",
  "SAOUHSC_01246" = "infB",
  "SAOUHSC_00475" = "pth"
)

genes_a_annoter <- intersect(names(noms_genes), rownames(res_translation))

if (length(genes_a_annoter) > 0) {
  text(
    log10(res_translation[genes_a_annoter, "baseMean"] + 1),
    res_translation[genes_a_annoter, "log2FoldChange"],
    labels = noms_genes[genes_a_annoter],
    pos = 3,
    cex = 0.8
  )
}

dev.off()


writeLines(
  capture.output(sessionInfo()),
  "results/deseq2/sessionInfo.txt"
)
