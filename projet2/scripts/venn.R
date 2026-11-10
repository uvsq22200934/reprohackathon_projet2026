# ============================================================
# PROJET 2 - Diagrammes de Venn
# Comparaison des transcrits différentiellement exprimés
# entre infection aiguë (J11) et chronique (J33)
# ============================================================

library(VennDiagram)
library(grid)

input_dir <- "results/deseq2"
output_dir <- "results/venn"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Lecture des identifiants des transcrits
read_transcripts <- function(filename) {
  data <- read.delim(
    file.path(input_dir, filename),
    stringsAsFactors = FALSE
  )

  if (!"transcript_id" %in% colnames(data)) {
    stop(paste("Colonne transcript_id absente dans", filename))
  }

  unique(data$transcript_id)
}

# Génération des diagrammes pour les deux critères statistiques
criteria <- c("raw", "FDR")
directions <- c("upregulated", "downregulated")

for (criterion in criteria) {
  for (direction in directions) {

    suffix <- if (criterion == "raw") "" else "_FDR"

    acute_file <- paste0("acute_", direction, suffix, ".tsv")
    chronic_file <- paste0("chronic_", direction, suffix, ".tsv")

    acute <- read_transcripts(acute_file)
    chronic <- read_transcripts(chronic_file)

    output_file <- file.path(
      output_dir,
      paste0("venn_", direction, "_", criterion, ".png")
    )

    venn.diagram(
      x = list(
        "Aigu (J11)" = acute,
        "Chronique (J33)" = chronic
      ),
      filename = output_file,
      imagetype = "png",
      height = 2000,
      width = 2000,
      resolution = 300,
      compression = "lzw",
      fill = c("skyblue", "salmon"),
      alpha = 0.5,
      cex = 1.4,
      cat.cex = 1.2,
      main = paste(
        "Transcrits", direction,
        "- critère", criterion
      ),
      main.cex = 1.4
    )

    cat("Diagramme créé :", output_file, "\n")
  }
}

cat("Génération des diagrammes de Venn terminée.\n")
