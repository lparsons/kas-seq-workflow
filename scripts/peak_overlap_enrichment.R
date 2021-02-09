#!/usr/bin/env Rscript --vanilla

library(BiocManager, quietly = TRUE)
library(ChIPseeker, quietly = TRUE)
library(readr, quietly = TRUE)
library(argparser, quietly = TRUE)

p <- arg_parser("KAS-Seq Peak Overlap Enrichment Analysis")
p <- add_argument(p, "--query_peak_file",
                  help = "Query peak files to compare targets against")
p <- add_argument(p, "--target_peak_files",
                  help = "Target peak files to compare against query peaks",
                  nargs = Inf)
p <- add_argument(p, "--txdb_file",
                  help = "File to load txdb using AnnotationDbi::loadDb()")
p <- add_argument(p, "--annotation_file",
                  help = paste0("GFF3 or GTF file of gene annotations used ",
                                "to build txdb"))
p <- add_argument(p, "--txdb",
                  help = "Name of txdb package to install from Bioconductor")

# Add an optional arguments
p <- add_argument(p, "--nshuffle", help = "Number of shuffle iterations",
                  default = 1000)
p <- add_argument(p, "--p_adjust_method", help = "pvalue adjustment method",
                  default = "BH")
p <- add_argument(p, "--output_tsv_file",
                  help = "Path and filename of output tsv file",
                  default = "peak_overlap_enrichment.tsv")
p <- add_argument(p, "--cores",
                  help = "number of cores (threads) to use",
                  default = NA)
p <- add_argument(p, "--chromlist",
                  help = paste0("List of chromosomes to include, others will ",
                                "be filtered out, new peaks stored in ",
                                "[peak_file].filtered"),
                  nargs = Inf, default = NA)

# Parse arguments (interactive, snakemake, or command line)
if (exists("snakemake")) {
  # Arguments via Snakemake
  argv <- parse_args(p, c(
    "--query_peak_file", snakemake@input[["query_peak_file"]],
    "--target_peak_files", snakemake@input[["target_peak_files"]],
    "--txdb_file", snakemake@input[["txdb_file"]],
    "--nshuffle", snakemake@params[["nshuffle"]],
    "--p_adjust_method", snakemake@params[["p_adjust_method"]],
    "--output_tsv_file", snakemake@output[["output_tsv_file"]],
    "--cores", snakemake@threads,
    "--chromlist", snakemake@params[["chromlist"]]
  ))
} else if (interactive()) {
  # Arguments supplied inline (for debug/testing when running interactively)
  print("Running interactively...")
  query_peak_file <- "results_2020-12-03/macs2/D701-lane1_peaks.broadPeak"
  target_peak_files <- c("results_2020-12-03/macs2/D702-lane1_peaks.broadPeak",
                         "results_2020-12-03/macs2/D703-lane1_peaks.broadPeak",
                         "results_2020-12-03/macs2/D704-lane1_peaks.broadPeak",
                         "results_2020-12-03/macs2/D705-lane1_peaks.broadPeak")
  nshuffle <- 50
  txdb_file <- "txdb.db"
  chromlist <- paste(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                       "12", "13", "14", "15", "16", "17", "18", "19", "20",
                       "21", "22", "23"), sep = ", ")
  argv <- parse_args(p, c("--query_peak_file", query_peak_file,
                          "--target_peak_files", target_peak_files,
                          "--nshuffle", nshuffle,
                          "--txdb_file", txdb_file,
                          "--chromlist", chromlist))
  print(argv)
} else {
  # Arguments from command line
  argv <- parse_args(p)
  print(argv)
}

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



# Number of cores
if (is.na(argv$cores)) {
  cores <- detectCores() - 1
} else {
  cores <- argv$cores
}

# Peak Overlap Enrichment
# Running with a single core can cause:
# Error in sample.int(length(x), size, replace, prob) :
#    cannot take a sample larger than the population when 'replace = FALSE'
# Running with multiple cores throws only a warning:
# Warning message:
# In mclapply(seq_along(idx), function(j) { :
#    all scheduled cores encountered errors in user code
# This results in incorrect output, thus we throw and error when we encounter
# that warning message
withCallingHandlers(
  peak_overlap_enrichment <- enrichPeakOverlap(
    queryPeak     = argv$query_peak_file,
    targetPeak    = unlist(argv$target_peak_files),
    TxDb          = txdb,
    pAdjustMethod = argv$p_adjust_method,
    nShuffle      = argv$nshuffle,
    chainFile     = NULL,
    verbose       = TRUE,
    mc.cores      = cores,
  ),
  warning = function(w) {
    if (grepl("encountered errors in user code", w$message)) {
      stop(paste0("Errors encounted executing enrichPeakOverlap: ", w$message))
    } else {
      message(w$message)
    }
  }
)

write_tsv(peak_overlap_enrichment, argv$output_tsv_file)
