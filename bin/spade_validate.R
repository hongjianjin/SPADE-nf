#!/usr/bin/env Rscript
# spade_validate.R - SPADE-nf input checks, run before the DE step so mistakes fail fast.
#
#   * metadata: SPADE ReadMeta() rules (ID, NEWNAME, GROUP; unique; [A-Za-z0-9._-] only)
#   * count matrix: gene columns present, sample columns match metadata IDs
#   * species: --species human|mouse checked against geneID (ENSG / ENSMUSG) or gene symbol case
#   * ERCC: number of ERCC- rows (the DE step decides ERCC vs TMM normalization)
#   * comparisons: given list checked against GROUP (>= 2 samples per side), or all pairs generated
#
# Outputs: <prefix>_comparisons.lst, <prefix>_validation.tsv

suppressPackageStartupMessages({
  library(optparse)
  library(stringr)
  library(edgeR)
})

option_list <- list(
  make_option(c("-i", "--rawCount"), type = "character", default = NA, help = "count matrix"),
  make_option(c("-m", "--metaData"), type = "character", default = NA, help = "metadata file"),
  make_option(c("-L", "--compList"), type = "character", default = NULL, help = "[optional] comparisons.lst"),
  make_option(c("-s", "--species"), type = "character", default = "human", help = "human, mouse or other [default %default]"),
  make_option(c("-o", "--outPrefix"), type = "character", default = NA, help = "outDir/prefix"),
  make_option("--spade_dir", type = "character", default = NULL, help = "SPADE checkout [default: $SPADE_HOME or /opt/spade]")
)
parser <- OptionParser(usage = "%prog -i <counts.txt> -m <meta.txt> -s <human|mouse|other> -o <prefix> [-L comparisons.lst]",
  option_list = option_list)
opt <- parse_args(parser)
if (is.na(opt$rawCount) || is.na(opt$metaData) || is.na(opt$outPrefix)) {
  print_help(parser)
  quit(status = 1)
}
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "spade_cli_helpers.R"))
spade_load(opt$spade_dir)

errors <- character()
notes  <- character()
fail <- function(...) errors <<- c(errors, paste0(...))

# --- metadata -----------------------------------------------------------------
metaObj <- ReadMeta(opt$metaData)
if (!isTRUE(metaObj$valid)) {
  cat(metaObj$errors, "\n")
  stop("Invalid metadata file: ", gsub("\n", " ", metaObj$errors))
}
meta <- metaObj$used

# --- counts -------------------------------------------------------------------
raw <- read.table(opt$rawCount, header = TRUE, sep = "\t", quote = "", fill = TRUE,
                  check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
annoCols <- c("geneid", "gene", "symbol", "genesymbol", "biotype", "annotationlevel")
isAnno <- tolower(colnames(raw)) %in% annoCols
if (!any(tolower(colnames(raw)) %in% c("geneid", "gene", "symbol", "genesymbol"))) {
  fail("Count matrix needs a geneID, geneSymbol, gene or symbol column.")
}
sampleCols <- colnames(raw)[!isAnno]
matched <- intersect(meta$ID, sampleCols)
missing <- setdiff(meta$ID, sampleCols)
if (length(matched) == 0) {
  fail("No metadata ID matches a count-matrix column. First columns: ", paste(head(sampleCols, 5), collapse = ", "))
} else if (length(missing) > 0) {
  notes <- c(notes, paste0("Metadata IDs not in the count matrix (dropped by SPADE): ", paste(missing, collapse = ", ")))
}
nonNumeric <- sampleCols[!vapply(raw[sampleCols], is.numeric, logical(1))]
if (length(intersect(nonNumeric, matched)) > 0) {
  fail("Non-numeric sample columns: ", paste(intersect(nonNumeric, matched), collapse = ", "))
}

symCol <- which(tolower(colnames(raw)) %in% c("genesymbol", "gene", "symbol"))[1]
idCol  <- which(tolower(colnames(raw)) == "geneid")[1]
symbols <- if (!is.na(symCol)) as.character(raw[[symCol]]) else if (!is.na(idCol)) as.character(raw[[idCol]]) else character()
nErcc <- sum(grepl("^ERCC-", symbols))
erccReads <- if (nErcc > 0 && length(matched) > 0) colSums(raw[grepl("^ERCC-", symbols), matched, drop = FALSE]) else numeric()

# --- species ------------------------------------------------------------------
species <- tolower(opt$species)
if (!species %in% c("human", "mouse", "other")) fail("--species must be human, mouse or other, not ", opt$species)
detected <- "unknown"
if (!is.na(idCol)) {
  ids <- as.character(raw[[idCol]])
  pHs <- mean(grepl("^ENSG[0-9]", ids)); pMm <- mean(grepl("^ENSMUSG[0-9]", ids))
  if (pHs > 0.5) detected <- "human" else if (pMm > 0.5) detected <- "mouse"
}
if (detected == "unknown" && length(symbols) > 0) {
  s <- symbols[!grepl("^ERCC-", symbols) & grepl("^[A-Za-z]", symbols)]
  s <- s[!grepl("^(LOC|LINC|MIR|SNOR)", s, ignore.case = TRUE)]
  if (length(s) > 100) {
    pUpper <- mean(s == toupper(s) & grepl("[A-Z]", s))
    pTitle <- mean(grepl("^[A-Z][a-z0-9]", s))
    if (pUpper > 0.8) detected <- "human" else if (pTitle > 0.8) detected <- "mouse"
  }
}
if (species %in% c("human", "mouse") && detected %in% c("human", "mouse") && detected != species) {
  fail("--species ", species, " but the gene IDs/symbols look ", detected, ". Fix --species (or use other with --custom_gmt).")
}
if (species %in% c("human", "mouse") && detected == "unknown") {
  notes <- c(notes, paste0("Could not confirm species ", species, " from gene IDs or symbols."))
}

# --- comparisons --------------------------------------------------------------
if (is.null(opt$compList) || !nzchar(opt$compList)) {
  comps <- GetComparisons(meta[meta$ID %in% matched, ])
  notes <- c(notes, paste0("Generated ", length(comps), " comparison(s) from GROUP."))
} else {
  comps <- gsub(" ", "", scan(opt$compList, character(), quote = "", quiet = TRUE))
  comps <- comps[nzchar(comps)]
}
grpN <- table(meta$GROUP[meta$ID %in% matched])
for (cmp in comps) {
  if (grepl(".*=\\(.*-.*\\)", cmp)) {
    g <- ExtractGroups(cmp)
  } else {
    g <- unlist(strsplit(cmp, "_[Vv][Ss]_"))
    if (length(g) != 2 || g[1] == g[2]) { fail("Invalid comparison '", cmp, "' (expected Exp_vs_Ctl)."); next }
  }
  absent <- g[!toupper(g) %in% toupper(names(grpN))]
  if (length(absent) > 0) { fail("Comparison '", cmp, "': group(s) not in metadata GROUP: ", paste(absent, collapse = ", ")); next }
  n <- grpN[match(toupper(g), toupper(names(grpN)))]
  if (any(n < 2)) fail("Comparison '", cmp, "': each group needs >= 2 samples (", paste(g, n, sep = "=", collapse = ", "), ").")
}

# --- write --------------------------------------------------------------------
prefix <- opt$outPrefix
dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
writeLines(comps, paste0(prefix, "_comparisons.lst"))
summ <- data.frame(
  key = c("samples_in_meta", "samples_matched", "genes", "groups", "comparisons", "species_requested",
          "species_detected", "ercc_rows", "ercc_samples_with_reads", "errors", "notes"),
  value = c(nrow(meta), length(matched), nrow(raw), paste(names(grpN), grpN, sep = "=", collapse = ";"),
            length(comps), species, detected, nErcc, sum(erccReads > 0), length(errors),
            paste(notes, collapse = " | ")))
spade_write_tsv(summ, paste0(prefix, "_validation.tsv"))
for (n in notes) cat("NOTE:", n, "\n")
if (length(errors) > 0) {
  cat("\nInput validation failed:\n", paste0("  - ", errors, collapse = "\n"), "\n", sep = "")
  quit(status = 1)
}
cat("\nInput validation passed:", length(matched), "samples,", length(comps), "comparison(s), species", species,
    "(detected:", detected, "), ERCC rows:", nErcc, "\n")
