
# ============================================================
# Differential expression analysis with DESeq2
# Persister vs Control
# ============================================================

library(DESeq2)

# 1. Sample information --------------------------------------

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

# 2. Read featureCounts files --------------------------------

files <- file.path(
  "results/counts",
  paste0(samples, ".txt")
)

if (!all(file.exists(files))) {
  stop("Some featureCounts files are missing.")
}

count_tables <- lapply(files, function(file) {
  read.delim(file, comment.char = "#", check.names = FALSE)
})

# Check that all samples have the same gene identifiers
gene_ids <- count_tables[[1]]$Geneid

for (i in seq_along(count_tables)) {
  if (!identical(gene_ids, count_tables[[i]]$Geneid)) {
    stop("Gene identifiers differ between samples.")
  }
}

# Extract the read counts (last column of each file)
count_matrix <- do.call(
  cbind,
  lapply(count_tables, function(x) {
    as.integer(x[[ncol(x)]])
  })
)

rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- samples

# 3. Create DESeq2 dataset -----------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ condition
)

# 4. Differential expression analysis -----------------------

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "Persister", "Control")
)

# Sort by adjusted p-value
res <- res[order(res$padj), ]

# 5. Export results ------------------------------------------

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

# 6. Global MA plot ------------------------------------------

pdf("results/deseq2/MA_plot.pdf", width = 8, height = 6)

plotMA(
  res,
  ylim = c(-5, 5),
  main = "Persister vs Control"
)

dev.off()

# 7. Save session information --------------------------------

writeLines(
  capture.output(sessionInfo()),
  "results/deseq2/sessionInfo.txt"
)
