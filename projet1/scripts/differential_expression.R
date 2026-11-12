#Production des figures avec R

# figure 3 supp :
# Expression différenciel avec DEseq2 :Persister vsControl



library(DESeq2)



echantillons <-c("SRR10379721", "SRR10379722", "SRR10379723",
  "SRR10379724", "SRR10379725", "SRR10379726") #echantillons

condition <- factor(c("Persister", "Persister", "Persister",
    "Control", "Control", "Control"),levels = c("Control", "Persister")) #on leur attribut leurs conditions



tab_data<- data.frame(condition =condition,row.names = samples)


#lectures des fichiers de comptage :

files<- file.path("results/counts",paste0(samples, ".txt"))

if (!all(file.exists(files))) {
  stop("fichiers manquants")
}

count_tables <- lapply(files,function(file) {
  read.delim(file, comment.char = "#", check.names = FALSE)})


gene_ids <- count_tables[[1]]$Geneid #verif que tous les echs ont le mm identififiant de gene

for (i in seq_along(count_tables)) {
  if (!identical(gene_ids, count_tables[[i]]$Geneid)) {
    stop("id de gene differe selon l'ech")
  }}

comptage_matrix <-do.call( #extraction des comptages (derniere colonne de chaque fichier)
  cbind,lapply(count_tables, function(x) {as.integer(x[[ncol(x)]])}))

rownames(comptage_matrix)<- gene_ids
colnames(comptage_matrix)<- echantillons


dds <- DESeqDataSetFromMatrix( #creation dataset DEqes2
  countData = count_matrix,
  colData = coldata,
  design =condition)



#expression differencielle

dds <- DESeq(dds)

res <-results(dds,contrast = c("condition","Persister", "Control"))

#tirage selon la pvalue
res <- res[order(res$padj), ]

#exporter les donnees:
write.csv(as.data.frame(res)file ="results/deseq2/differential_expression.csv")

write.table(counts(dds, normalized = TRUE),file = "results/deseq2/normalized_counts.tsv",sep = "\t",quote = FALSEcol.names = NA)

saveRDS(dds, "results/deseq2/dds.rds")


#MA-plot
pdf("results/deseq2/MA_plot.pdf",width = 8,height= 6)

plotMA(res,ylim = c(-5, 5),main = "Persister vs Control")

dev.off()

writeLines(
  capture.output(sessionInfo()),
  "results/deseq2/sessionInfo.txt")
