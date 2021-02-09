#!/usr/bin/env Rscript --vanilla

library(BiocManager, quietly=TRUE)
library(ChIPseeker, quietly=TRUE)
library(org.Hs.eg.db, quietly=TRUE)
library(cowplot, quietly=TRUE)
library(readr, quietly=TRUE)
library(argparser, quietly=TRUE)
# library(clusterProfiler, quietly=TRUE)

p <- arg_parser("KAS-Seq Peak Annotation and Comparison")
p <- add_argument(p, "--peak_files", help="Peak files to annotate and compare",
                  nargs=Inf)
p <- add_argument(p, "--txdb_file",
                  help="File to load txdb using AnnotationDbi::loadDb()")
p <- add_argument(p, "--annotation_file",
                  help="GFF3 or GTF file of gene annotations used to build txdb")
p <- add_argument(p, "--txdb",
                  help="Name of txdb package to install from Bioconductor")

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
p <- add_argument(p, "--venn_diagram",
                  help="Venn digagram of annotated genes per sample pdf filename",
                  default="annotationVennDiagram.pdf")

# Parse arguments (interactive, snakemake, or command line)
if (exists("snakemake")) {
  # Arguments via Snakemake
  argv <- parse_args(p, c(
    "--peak_files", snakemake@input[["peak_files"]],
    "--txdb_file", snakemake@input[["txdb_file"]],
    "--names", snakemake@params[["names"]],
    "--output_dir", snakemake@params[["output_dir"]],
    "--annotation_distribution_plot", snakemake@output[["annotation_distribution_plot"]],
    "--venn_diagram", snakemake@output[["venn_diagram"]],
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
  txdb_file <- "txdb.db"
  argv <- parse_args(p, c("--peak_files", peak_files,
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
  peakFileNames <- argv$names
} else {
  peakFileNames <- sapply(argv$peak_files, basename)
}
names(argv$peak_files) <- peakFileNames

# Output directory
if (!dir.exists(argv$output_dir)) {
  dir.create(argv$output_dir, recursive = TRUE)
}

# Get txdb object
if (!is.na(argv$txdb)) {
  # Load (install if needed) txdb from bioconductor
  library(pacman, quietly=TRUE)
  pacman::p_load(argv$txdb, character.only = TRUE)
  txdb <- eval(parse(text = argv$txdb))
} else if (!is.na(argv$txdb_file)) {
  # Load txdb
  library(AnnotationDbi, quietly=TRUE)
  txdb <- AnnotationDbi::loadDb(argv$txdb_file)
} else if (!is.na(argv$annotation_file)) {
  # Create txdb object from supplied annotation file
  library(GenomicFeatures, quietly=TRUE)
  txdb <- GenomicFeatures::makeTxDbFromGFF(argv$annotation_file)
} else {
  stop("Must specify one of --txdb, --txdb_file, or --annotation_file")
}


# Peak Annotation
# TODO Provide config parameter for annoDb
peakAnnoList <- lapply(argv$peak_files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), annoDb="org.Hs.eg.db",
                       verbose=FALSE)
lapply(names(peakAnnoList), function(name) {
  filebase = file.path(argv$output_dir, basename(argv$peak_files[[name]]))
  write_tsv(as.data.frame(peakAnnoList[[name]]),
            path=paste(filebase, ".annotated.tsv.gz", sep=""))
  sink(file = paste(filebase, ".annotated.summary.txt", sep=""))
  print(peakAnnoList[[name]])
  sink()
})
saveRDS(peakAnnoList, 
     file = argv$peak_annotation_list_rdata)

# Peak annotation distribution plot
peakAnnotationDistributionPlot <- plotAnnoBar(peakAnnoList)
save_plot(
  filename = argv$annotation_distribution_plot,
  plot = peakAnnotationDistributionPlot,
  base_height = 14,
  base_width = 14
)

# Functional profiles comparison
# NOTE: Not all data return enrichment
# genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
# names(genes) = sub("_", "\n", names(genes))
# compKEGG <- compareCluster(geneCluster   = genes,
#                            fun           = "enrichKEGG",
#                            pvalueCutoff  = 0.05,
#                            pAdjustMethod = "BH")
# dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")

# Venn Diagram of gene annotated per sample
# Note: vennplot does not return a ggplot object
pdf(
  file = argv$venn_diagram,
  height = 14,
  width = 14
)
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)
dev.off()
