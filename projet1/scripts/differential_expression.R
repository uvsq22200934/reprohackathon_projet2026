# Analyse différentielle avec DESeq2 : Bactéries persistantes vscontrol
#figure 

library(DESeq2)

# echantillons et conditions
samples= c("SRR10379721", "SRR10379722", "SRR10379723",
    "SRR10379724","SRR10379725", "SRR10379726")

condition= factor(
    c("Persister", "Persister", "Persister",
      "Control", "Control", "Control"),
    levels = c("Control", "Persister"))

coldata= data.frame(
    condition = condition,
    row.names = samples)

#lcture des comptages
files=file.path(
    "results/counts",
    paste0(samples, ".txt"))

if (!all(file.exists(files))) {
    stop("Un ou plusieurs fichiers de comptage sont manquants.")}

count_tables= lapply(files, function(file) {
    read.delim(file,comment.char = "#", check.names = FALSE)
})

gene_ids =count_tables[[1]]$Geneid


count_matrix= do.call(
    cbind,
    lapply(count_tables, function(x) {
        as.integer(x[[ncol(x)]])
    }))

rownames(count_matrix) = gene_ids
colnames(count_matrix)= samples


# Création du dataset DESeq2

dds =DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = ~ condition
)

# Analyse différentielle

dds= DESeq(dds)

res =results(
    dds,
    contrast = c("condition", "Persister", "Control"))

res =res[order(res$padj), ]


dir.create("results/deseq2", recursive= TRUE, showWarnings = FALSE)

write.csv(
    as.data.frame(res),
    file = "results/deseq2/differential_expression.csv"
)

#nombre de gènes différentiellement exprimés :

res_df =as.data.frame(res)

significant = !is.na(res_df$padj) & res_df$padj < 0.05

upregulated =significant &
               !is.na(res_df$log2FoldChange) &
               res_df$log2FoldChange > 0

downregulated= significant &
                 !is.na(res_df$log2FoldChange) &
                 res_df$log2FoldChange < 0

summary= data.frame(
    category = c(
        "Total genes",
        "Genes tested by DESeq2",
        "Differentially expressed genes (padj < 0.05)",
        "Upregulated in Persisters",
        "Downregulated in Persisters"
    ),
    count = c(
        nrow(res_df),
        sum(!is.na(res_df$padj)),
        sum(significant),
        sum(upregulated),
        sum(downregulated)
    )
)

write.csv(
    summary,
    file = "results/deseq2/differential_expression_summary.csv",
    row.names = FALSE
)

print(summary)

# Export des comptages normalisés

write.table(
    counts(dds, normalized = TRUE),
    file = "results/deseq2/normalized_counts.tsv",
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

saveRDS(dds, "results/deseq2/dds.rds")


#MA plot de tous les gènes
data = as.data.frame(res)

#Retirer les valeurs non exploitables

valid= !is.na(data$baseMean) &
         !is.na(data$log2FoldChange) &
         is.finite(data$baseMean) &
         is.finite(data$log2FoldChange)

data = data[valid, ]

#significativité statistique
significant = !is.na(data$padj) & data$padj < 0.05

point_colors =ifelse(significant, "red", "grey70")

#position des gènes

x <- log10(data$baseMean + 1)
y <- data$log2FoldChange

pdf("results/deseq2/MA_plot.pdf", width= 9, height = 7)

plot(
    x, y,
    pch =16,
    col = point_colors,
    cex = 0.7,
    xlab = "log10(mean normalized counts + 1)",
    ylab ="log2 fold change (Persister /control)",
    main= "Persister vs Control (tous les gènes)",
    ylim = c(-5, 5))

abline(h = 0,lty = 2, col = "black")

legend(
    "topright",
    legend = c(
        "Significatif(padj < 0.05)",
        "Non significatif"
    ),
    col = c("red", "grey70"),
    pch =c(16, 16),
    bty = "n",
    cex = 0.8
)

dev.off()

writeLines(
    capture.output(sessionInfo()),
    "results/deseq2/sessionInfo.txt"
)