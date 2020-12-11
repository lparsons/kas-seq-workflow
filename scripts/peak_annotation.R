#!/usr/bin/env Rscript --vanilla

library(BiocManager, quietly=TRUE)
library(ChIPseeker, quietly=TRUE)
library(org.Hs.eg.db, quietly=TRUE)
library(GenomicFeatures, quietly=TRUE)
library(cowplot, quietly=TRUE)
library(readr, quietly=TRUE)
library(argparser, quietly=TRUE)

p <- arg_parser("KAS-Seq Peak Annotation and Comparison")
p <- add_argument(p, "--peak_files", help="Peak files to annotate and compare",
                  nargs=Inf)
p <- add_argument(p, "--annotation_file", help="GFF3 or GTF file of gene annotations")

# Add an optional arguments
p <- add_argument(p, "--names", help="Sample names for each peak file",
                  nargs=Inf)
p <- add_argument(p, "--output_dir",
                  help="Directory for output files",
                  default="peak_annotation")
p <- add_argument(p, "--annotation_distribution_plot",
                  help="Peak annotation distribution barplot filename",
                  default="annotationDistributionPlot.pdf")
p <- add_argument(p, "--peak_annotation_list_rdata",
                  help="Peak annotation list Rdata file",
                  default="peakAnnoList.Rdata")

# Parse arguments (interactive, snakemake, or command line)
if (exists("snakemake")) {
  # Arguments via Snakemake
  argv <- parse_args(p, c(
    "--peak_files", snakemake@input[["peak_files"]],
    "--annotation_file", snakemake@input[["annotation_file"]],
    "--names", snakemake@params[["names"]],
    "--output_dir", snakemake@params[["output_dir"]],
    "--annotation_distribution_plot", snakemake@output[["annotation_distribution_plot"]],
    "--peak_annotation_list_rdata", snakemake@output[["peak_annotation_list_rdata"]]
  ))
} else if (interactive()) {
  # Arguments supplied inline (for debug/testing when running interactively)
  print("Running interactively...")
  peak_files <- c("results_2020-12-03/macs2/D701-lane1_peaks.broadPeak",
                  "results_2020-12-03/macs2/D702-lane1_peaks.broadPeak", 
                  "results_2020-12-03/macs2/D703-lane1_peaks.broadPeak", 
                  "results_2020-12-03/macs2/D704-lane1_peaks.broadPeak", 
                  "results_2020-12-03/macs2/D705-lane1_peaks.broadPeak")
  names <- c("D701", "D702", "D703", "D704", "D705")
  annotation_file <- "genomes/hg38/annotation/Homo_sapiens.GRCh38.101.gtf"
  argv <- parse_args(p, c("--peak_files", peak_files,
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

# Output directory
if (!dir.exists(argv$output_dir)) {
  dir.create(argv$output_dir, recursive = TRUE)
}

# Create txdb object from supplied annotation file
txdb <- makeTxDbFromGFF(argv$annotation_file)

# Peak Annotation
# TODO Provide config parameter for annoDb
peakAnnoList <- lapply(argv$peak_files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), annoDb="org.Hs.eg.db",
                       verbose=FALSE)
lapply(names(peakAnnoList), function(name) {
  filebase = file.path(argv$output_dir, basename(argv$peak_files[[name]]))
  write_tsv(as.data.frame(peakAnnoList[[name]]),
            file=paste(filebase, ".annotated.tsv.gz", sep=""))
  sink(file = paste(filebase, ".annotated.summary.txt", sep=""))
  print(peakAnnoList[[name]])
  sink()
})
saveRDS(peakAnnoList, 
     file = argv$peak_annotation_list_rdata)


peakAnnotationDistributionPlot <- plotAnnoBar(peakAnnoList)
save_plot(
  filename = argv$annotation_distribution_plot,
  plot = peakAnnotationDistributionPlot,
  base_height = 14,
  base_width = 14
)
