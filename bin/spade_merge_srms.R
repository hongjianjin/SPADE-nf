#!/usr/bin/env Rscript
# spade_merge_srms.R - CAB internal: merge AutoMapper deliverables of several SRMs.
#
# Uses MergeRSEMs() and MergeReports() from RNAseqERCC v1.1 (ported in spade_cli_helpers.R):
# count / TPM / FPKM matrices are inner-joined on geneID (same rows and annotation columns
# required); .report files are stacked. Unlike RNAseqERCC's MergeSRMs(), no YAML is rewritten:
# the merged files are written as <prefix>.report and <prefix>_RSEM_gene_{count,TPM,FPKM}.txt.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})
option_list <- list(
  make_option("--counts", type = "character", default = NA, help = "comma-separated RSEM count files (one per SRM)"),
  make_option("--reports", type = "character", default = "", help = "[optional] comma-separated .report files"),
  make_option("--tpm", type = "character", default = "", help = "[optional] comma-separated TPM files"),
  make_option("--fpkm", type = "character", default = "", help = "[optional] comma-separated FPKM files"),
  make_option(c("-o", "--outPrefix"), type = "character", default = NA, help = "outDir/prefix, e.g. mergedSRMs_2")
)
opt <- parse_args(OptionParser(usage = "%prog --counts a.txt,b.txt [--reports a.report,b.report --tpm ... --fpkm ...] -o prefix",
  option_list = option_list))
if (is.na(opt$counts) || is.na(opt$outPrefix)) stop("--counts and --outPrefix are required")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "spade_cli_helpers.R"))

splitList <- function(x) { x <- trimws(unlist(strsplit(x, ","))); x[nzchar(x)] }
countFiles <- splitList(opt$counts)
if (length(countFiles) < 2) stop("Need at least two count files to merge.")
prefix <- opt$outPrefix
dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)

# Sample IDs must be unique across SRMs, otherwise the join silently renames columns.
annoCols <- c("geneID", "geneSymbol", "bioType", "annotationLevel")
samples <- lapply(countFiles, function(f) setdiff(colnames(read.delim(f, nrows = 1, check.names = FALSE)), annoCols))
dups <- unique(unlist(samples)[duplicated(unlist(samples))])
if (length(dups) > 0) stop("Sample IDs appear in more than one SRM: ", paste(dups, collapse = ", "))

mergeOne <- function(files, label) {
  if (length(files) == 0) return(invisible(NULL))
  if (length(files) != length(countFiles)) stop("Got ", length(files), " ", label, " files for ", length(countFiles), " SRMs.")
  out <- paste0(prefix, "_RSEM_gene_", label, ".txt")
  ok <- MergeRSEMs(files, out)
  if (!isTRUE(ok)) stop("Failed to merge ", label, " files.")
  n_in <- vapply(files, function(f) length(readLines(f)) - 1L, integer(1))
  n_out <- length(readLines(out)) - 1L
  cat("\n", basename(out), "[saved]: ", n_out, " genes (inputs: ", paste(n_in, collapse = ", "), ")\n", sep = "")
}
mergeOne(countFiles, "count")
mergeOne(splitList(opt$tpm), "TPM")
mergeOne(splitList(opt$fpkm), "FPKM")

reports <- splitList(opt$reports)
if (length(reports) > 0) {
  if (length(reports) < 2) stop("Need a .report file for every SRM (or none).")
  ok <- MergeReports(reports, paste0(prefix, ".report"))
  if (!isTRUE(ok)) stop("Failed to merge .report files.")
  cat("\n", basename(prefix), ".report [saved]\n", sep = "")
}
