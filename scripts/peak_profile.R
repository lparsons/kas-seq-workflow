#!/usr/bin/env Rscript --vanilla

library(BiocManager, quietly=TRUE)
library(ChIPseeker, quietly=TRUE)
library(GenomicFeatures, quietly=TRUE)
library(cowplot, quietly=TRUE)
library(readr, quietly=TRUE)
library(argparser, quietly=TRUE)

p <- arg_parser("Profile of peaks binding to TSS regions")
p <- add_argument(p, "--peak_files", help="Peak files to annotate and compare",
                  nargs=Inf)
p <- add_argument(p, "--annotation_file", help="GFF3 or GTF file of gene annotations")

# Add an optional arguments
p <- add_argument(p, "--names", help="Sample names for each peak file",
                  nargs=Inf)
p <- add_argument(p, "--tag_profile_plot",
                  help="Average peak profile plot filename",
                  default="avgProfilePlot.pdf")
p <- add_argument(p, "--tag_heatmap_plot",
                  help="Tag heatmap plot PDF filename",
                  default="tagHeatmapPlot.pdf")
p <- add_argument(p, "--tag_matrix_list",
                  help="List of outputs from getTagMatrix CHiPseeker function",
                  default="tagMatrixList.Rdata")

# Parse arguments (interactive, snakemake, or command line)
if (exists("snakemake")) {
  # Arguments via Snakemake
  argv <- parse_args(p, c(
    "--peak_files", snakemake@input[["peak_files"]],
    "--annotation_file", snakemake@input[["annotation_file"]],
    "--names", snakemake@params[["names"]],
    "--tag_profile_plot", snakemake@output[["tag_profile_plot"]],
    "--tag_heatmap_plot", snakemake@output[["tag_heatmap_plot"]],
    "--tag_matrix_list", snakemake@output[["tag_matrix_list"]]
  ))
} else if (interactive()) {
  # Arguments supplied inline (for debug/testing when running interactively)
  print("Running interactively...")
  input_file <- c("results_2020-12-03/macs2/D701-lane1_peaks.broadPeak",
                  "results_2020-12-03/macs2/D702-lane1_peaks.broadPeak", 
                  "results_2020-12-03/macs2/D703-lane1_peaks.broadPeak", 
                  "results_2020-12-03/macs2/D704-lane1_peaks.broadPeak", 
                  "results_2020-12-03/macs2/D705-lane1_peaks.broadPeak")
  names <- c("D701", "D702", "D703", "D704", "D705")
  annotation_file <- "genomes/hg38/annotation/Homo_sapiens.GRCh38.101.gtf"
  argv <- parse_args(p, c("--peak_files", input_file,
                          "--names", names,
                          "--annotation_file", annotation_file))
  print(argv)
} else {
  # Arguments from command line
  argv <- parse_args(p)  
  print(argv)
}

# Set names
if (!anyNA(argv$names)) {
  peakFileNames <- argv$names
} else {
  peakFileNames <- sapply(argv$peak_files, basename)
}
names(argv$peak_files) <- peakFileNames

# Create txdb object from supplied annotation file
txdb <- makeTxDbFromGFF(argv$annotation_file)

# Profile of peaks binding to TSS regions

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(argv$peak_files, getTagMatrix, windows=promoter)
saveRDS(tagMatrixList, file = argv$tag_matrix_list)

# Average tag profile
tagProfilePlot <- plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
save_plot(
  filename = argv$tag_profile_plot,
  plot = tagProfilePlot,
  base_height = 14,
  base_width = 14
)

# Tag heatmap
# Note: tagHeatMap does not return a ggplot object, we must save it differently
pdf(
  file = argv$tag_heatmap_plot,
  height = 14,
  width = 14
)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
dev.off()

