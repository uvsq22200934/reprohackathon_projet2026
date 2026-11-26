#figure 3c de l'article : MA plot des gènes impliqués dans la traduction

library(ggplot2)

res=read.csv("results/deseq2/differential_expression.csv",
  row.names =1) #chargement résultats de l'analyse différentielle qu'on a effectuée avant

res$gene_id = rownames(res)


translation_genes=scan(
  "data/kegg/translation_related_genes_gff.txt",
  what ="character", quiet = TRUE)

res =res[res$gene_id %in% translation_genes, ] 

#charger la liste des aa-ARNt synthétases
aminoacyl_genes= scan("data/kegg/aminoacyl_trna_synthetases.txt",
  what ="character",
  quiet = TRUE)

res$aminoacyl= res$gene_id %in% aminoacyl_genes


res$significatif =!is.na(res$padj) & res$padj < 0.05

res$log2BaseMean = log2(res$baseMean + 1)

#chargement des noms des gènes à annoter
annotations=read.delim(
  "data/kegg/target_gene_annotations.tsv")

res$symbol= annotations$symbol[
  match(res$gene_id, annotations$gene_id)
]


#Graphique

res$categorie= ifelse(
  res$significatif,
  "significatif",
  "Non significatif")

p =ggplot(res, aes(x = log2BaseMean, y = log2FoldChange)) +

  #points gris et rouges pour tous les gènes
  geom_point(
    aes(color = categorie),
    size = 1.8) +

  #cercles noirs autour des aa-ARNt synthétases
  geom_point(
    data = subset(res, aminoacyl),
    mapping = aes(
      x = log2BaseMean,
      y = log2FoldChange,
      shape = "Aminoacyl-ARNt synthetases"
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 2.5,
    stroke = 0.8
  ) +

  # Définir les couleurs de la légende
  scale_color_manual(
    name = NULL,
    values = c(
      "significatif" = "red",
      "Non significatif" = "grey65")) +

  scale_shape_manual(
    name = NULL,
    values = c("Aminoacyl-ARNt synthetases" = 1)
  ) +

  geom_hline(
    yintercept= 0,
    linetype = "dashed",
    color = "grey50") +

  #mettre le nom sur les 6 gènes
  geom_text(
    data = subset(res, !is.na(symbol)),
    aes(label= symbol),
    hjust=-0.2,
    vjust= -0.5,
    size = 3) +

  labs(
    title = "Expression différentielle des gènes impliqués dans la traduction",
    x = "log2 des comptages normalisés moyen (+ 1)",
    y = "log2 fold change") +

  theme_classic() +

  theme(legend.position = "right")

ggsave("results/deseq2/MA_plot_translation_genes.pdf", plot = p,
  width = 7,
  height = 5)