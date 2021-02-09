#!/usr/bin/env Rscript --vanilla

library(BiocManager, quietly=TRUE)
library(AnnotationDbi, quietly=TRUE)
library(readr, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(argparser, quietly=TRUE)

p <- arg_parser("Get a TXDB object and save for use in other scripts")
p <- add_argument(p, "txdb_input",
                  help=paste0("Name of txdb package to install from Bioconductor ",
                              "or annotation file (GFF3 or GTF) to use to ",
                              "build a TXDB object"))
p <- add_argument(p, "--chrom_info",
                  help=paste0("Chromosome info file with columns: ",
                  "chrom, length, is_circular (optional)"))
p <- add_argument(p, "--output",
                  help="File to save txdb object in",
                  default="txdb.db")

# Parse arguments (interactive, snakemake, or command line)
if (exists("snakemake")) {
  # Arguments via Snakemake
  argv <- parse_args(p, c(
    snakemake@input[["txdb_input"]],
    "--chrom_info", snakemake@input[["chrom_info"]],
    "--output", snakemake@output[[1]]
  ))
} else if (interactive()) {
  # Arguments supplied inline (for debug/testing when running interactively)
  print("Running interactively...")
  # txdb_input <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
  txdb_input <- "genomes/hg38/annotation/Homo_sapiens.GRCh38.101.gtf"
  chrom_info <- "genomes/hg38/genome/fasta/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai"
  argv <- parse_args(p, c(txdb_input, 
                          "--chrom_info", chrom_info))
  print(argv)
} else {
  # Arguments from command line
  argv <- parse_args(p)  
  print(argv)
}

# Build txdb if txdb_input is file
if (file.exists(argv$txdb_input)) {
  # Create txdb object from supplied annotation file
  library(GenomicFeatures, quietly=TRUE)
  if (!is.na(argv$chrom_info)) {
    # Get chrom info
    is_circular_column <-
      count_fields(argv$chrom_info, tokenizer_tsv())[1] > 2
    if (is_circular_column) {
      col_names <- c("chrom", "length", "is_circular")
      col_spec <- cols_only(X1 = col_character(),
                            X2 = col_integer(),
                            X3 = col_guess())
    } else {
      col_names <- c("chrom", "length")
      col_spec <- cols(X1 = col_character(),
                       X2 = col_integer())
    }
    chrom_info <- read_tsv(argv$chrom_info,
                           col_names = FALSE,
                           col_types = col_spec)
    if (is.logical(chrom_info$X3)) {
      chrom_info <-
        chrom_info %>% select(c(
          chrom = X1,
          length = X2,
          is_circular = X3
        ))
    } else {
      chrom_info <- chrom_info %>% select(c(chrom = X1, length = X2))
    }
    txdb <- makeTxDbFromGFF(argv$txdb_input, chrominfo = chrom_info)
  } else {
    txdb <- makeTxDbFromGFF(argv$txdb_input)
  }
} else {
  library(pacman, quietly=TRUE)
  # Load (install if needed) txdb from bioconductor
  pacman::p_load(argv$txdb_input, character.only = TRUE)
  txdb <- eval(parse(text = argv$txdb_input))  
}

# Save txdb
AnnotationDbi::saveDb(txdb, argv$output)
