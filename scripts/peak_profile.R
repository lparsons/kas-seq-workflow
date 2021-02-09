#!/usr/bin/env Rscript --vanilla

library(BiocManager, quietly = TRUE)
library(ChIPseeker, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(readr, quietly = TRUE)
library(argparser, quietly = TRUE)

p <- arg_parser("Profile of peaks binding to TSS regions")
p <- add_argument(p, "--peak_files",
                  help = "Peak files to annotate and compare",
                  nargs = Inf)
p <- add_argument(p, "--txdb_file",
                  help = "File to load txdb using AnnotationDbi::loadDb()")
p <- add_argument(p, "--annotation_file",
                  help = paste0("GFF3 or GTF file of gene annotations used ",
                                "to build txdb"))
p <- add_argument(p, "--txdb",
                  help = "Name of txdb package to install from Bioconductor")

# Add an optional arguments
p <- add_argument(p, "--names", help = "Sample names for each peak file",
                  nargs = Inf)
p <- add_argument(p, "--tag_profile_plot",
                  help = "Average peak profile plot filename",
                  default = "avgProfilePlot.pdf")
p <- add_argument(p, "--tag_heatmap_plot",
                  help = "Tag heatmap plot PDF filename",
                  default = "tagHeatmapPlot.pdf")
p <- add_argument(p, "--tag_matrix_list",
                  help = paste0("List of outputs from getTagMatrix ",
                                "CHiPseeker function"),
                  default = "tag_matrix_list.Rdata")

# Parse arguments (interactive, snakemake, or command line)
if (exists("snakemake")) {
  # Arguments via Snakemake
  argv <- parse_args(p, c(
    "--peak_files", snakemake@input[["peak_files"]],
    "--txdb_file", snakemake@input[["txdb_file"]],
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
  txdb_file <- "txdb.db"
  argv <- parse_args(p, c("--peak_files", input_file,
                          "--names", names,
                          "--txdb_file", txdb_file))
  print(argv)
} else {
  # Arguments from command line
  argv <- parse_args(p)
  print(argv)
}

# Set names
if (!anyNA(argv$names)) {
  peak_file_names <- argv$names
} else {
  peak_file_names <- sapply(argv$peak_files, basename)
}
names(argv$peak_files) <- peak_file_names

# Get txdb object
if (!is.na(argv$txdb)) {
  # Load (install if needed) txdb from bioconductor
  library(pacman, quietly = TRUE)
  pacman::p_load(argv$txdb, character.only = TRUE)
  txdb <- eval(parse(text = argv$txdb))
} else if (!is.na(argv$txdb_file)) {
  # Load txdb
  library(AnnotationDbi, quietly = TRUE)
  txdb <- AnnotationDbi::loadDb(argv$txdb_file)
} else if (!is.na(argv$annotation_file)) {
  # Create txdb object from supplied annotation file
  library(GenomicFeatures, quietly = TRUE)
  txdb <- GenomicFeatures::makeTxDbFromGFF(argv$annotation_file)
} else {
  stop("Must specify one of --txdb, --txdb_file, or --annotation_file")
}


# Profile of peaks binding to TSS regions
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tag_matrix_list <- lapply(argv$peak_files, getTagMatrix, windows = promoter)
saveRDS(tag_matrix_list, file = argv$tag_matrix_list)

# Average tag profile
tag_profile_plot <- plotAvgProf(tag_matrix_list, xlim = c(-3000, 3000))
save_plot(
  filename = argv$tag_profile_plot,
  plot = tag_profile_plot,
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
tagHeatmap(tag_matrix_list, xlim = c(-3000, 3000), color = NULL)
dev.off()
