library(tximport)
library(DESeq2)
library(ggplot2)

metadata <- read.delim("config/samples.tsv")
rownames(metadata) <- metadata$sample
metadata$condition <- factor(metadata$condition)
metadata$day <- factor(metadata$day)

#quantifications Salmon
files <- file.path("results/salmon", metadata$sample, "quant.sf")
names(files) <- metadata$sample
txi <- tximport(files, type = "salmon", txOut = TRUE)

#objet des
dds <- DESeqDataSetFromTximport(txi, metadata, design = ~ day + condition)
#normaliser les objets deseq2 
vsd <- vst(dds, blind = TRUE)
p <- plotPCA(vsd, intgroup = c("condition", "day"))

dir.create("results/pca", showWarnings = FALSE)
ggsave("results/pca/PCA.png", p, width = 8, height = 6)
