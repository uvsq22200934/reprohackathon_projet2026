library(VennDiagram)

read_transcripts <- function(filename) {
  data <- read.delim(file.path("results/deseq2", filename))
  unique(data$transcript_id)
}

for (criterion in c("raw", "FDR")) {

  suffix <- if (criterion == "raw") "" else "_FDR"

  sets <- list(
    Acute_upregulated = read_transcripts(
      paste0("acute_upregulated", suffix, ".tsv")
    ),
    Acute_downregulated = read_transcripts(
      paste0("acute_downregulated", suffix, ".tsv")
    ),
    Chronic_upregulated = read_transcripts(
      paste0("chronic_upregulated", suffix, ".tsv")
    ),
    Chronic_downregulated = read_transcripts(
      paste0("chronic_downregulated", suffix, ".tsv")
    )
  )

  venn.diagram(
    x = sets,
    filename = paste0("results/venn/figure2B_", criterion, ".png"),
    imagetype = "png",
    height = 2500,
    width = 2500,
    resolution = 300,
    fill = c("skyblue", "yellow", "salmon", "lightgreen"),
    alpha = 0.5,
    cex = 1.3,
    cat.cex = 1,
    margin = 0.1
  )
}
