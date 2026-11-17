#!/usr/bin/env Rscript
# spade_de.R - SPADE-nf differential expression step (QC, ERCC QC, limma-voom DE).
#
# Command line and outputs follow RNAseqERCC v1.1 stdRNAseq.ERCC.Main.R
# (https://github.com/hongjianjin/RNAseqERCC); the analysis uses the SPADE v3
# functions (https://github.com/hongjianjin/SPADE), loaded by spade_cli_helpers.R.
#
# Normalization:
#   --ERCC TRUE  and ERCC- rows present -> ERCC (spike-in) normalization
#   --ERCC TRUE  and no ERCC- rows      -> library-size (TMM) normalization, with a note in log and report
#   --ERCC FALSE                        -> library-size (TMM) normalization

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(stringr)
  library(patchwork)
  library(ggpubr)
  library(limma)
  library(edgeR)
  library(reshape2)
  library(ggrepel)
  library(parallel)
})

option_list <- list(
  make_option(c("-i", "--rawCount"), type = "character", default = NA,
    help = "count matrix (geneID, geneSymbol, [bioType, annotationLevel], one column per sample)"),
  make_option(c("-r", "--report"), type = "character", default = NULL,
    help = "[optional] AutoMapper ${SRM}.report file for mapping QC"),
  make_option(c("-m", "--metaData"), type = "character", default = NA,
    help = "metadata file with ID, NEWNAME and GROUP columns [COLOR, BATCH, GENDER, SHAPE]"),
  make_option(c("-L", "--compList"), type = "character", default = NULL,
    help = "[optional] comparisons.lst, one Exp_vs_Ctl per line; all pairs of GROUP if omitted"),
  make_option(c("-C", "--countCutoff"), type = "numeric", default = 0,
    help = "raw count cutoff; 0 = auto [default %default]"),
  make_option(c("-c", "--cpmCutoff"), type = "numeric", default = 0,
    help = "CPM cutoff; 0 = auto, -1 = no low-expression filter [default %default]"),
  make_option(c("-o", "--outPrefix"), type = "character", default = NA,
    help = "outDir/prefix for output files"),
  make_option(c("-f", "--filter"), type = "character", default = "TRUE",
    help = "filter non-coding genes by bioType/annotationLevel [default %default]"),
  make_option(c("-E", "--ERCC"), type = "character", default = "TRUE",
    help = "use ERCC normalization when ERCC- rows are present [default %default]"),
  make_option(c("-F", "--log2FCcutoff"), type = "numeric", default = 1,
    help = "log2FC cutoff to define DEGs [default %default]"),
  make_option(c("-P", "--PvalCutoff"), type = "numeric", default = 0.05,
    help = "P.Value cutoff to define DEGs [default %default]"),
  make_option(c("-q", "--qCutoff"), type = "numeric", default = 0.05,
    help = "FDR cutoff to define DEGs [default %default]"),
  make_option(c("-n", "--topn"), type = "integer", default = 20,
    help = "number of top up and down DEGs [default %default]"),
  make_option("--numGenePCA", type = "integer", default = 3000,
    help = "number of most variable genes for QC PCA [default %default]"),
  make_option("--ncores", type = "integer", default = 1,
    help = "CPUs for running comparisons in parallel [default %default]"),
  make_option("--spade_dir", type = "character", default = NULL,
    help = "SPADE checkout with module/ and data/ [default: $SPADE_HOME or /opt/spade]")
)
parser <- OptionParser(
  usage = "%prog -i <counts.txt> -m <meta.txt> -o <outDir/prefix> [-L comparisons.lst -r SRM.report -E TRUE -f TRUE ...]",
  description = "SPADE differential expression with optional ERCC spike-in normalization.",
  option_list = option_list
)
opt <- parse_args(parser)
if (is.na(opt$rawCount) || is.na(opt$metaData) || is.na(opt$outPrefix)) {
  print_help(parser)
  quit(status = 1)
}

#--------------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "spade_cli_helpers.R"))
spade_load(opt$spade_dir)

countFile    <- normalizePath(opt$rawCount, mustWork = TRUE)
metaFile     <- normalizePath(opt$metaData, mustWork = TRUE)
reportFile   <- if (!is.null(opt$report) && nzchar(opt$report)) normalizePath(opt$report, mustWork = TRUE) else NULL
filterFlag   <- spade_flag(opt$filter)
flag_ERCC    <- spade_flag(opt$ERCC)
cutoff_count <- opt$countCutoff
cutoff_cpm   <- opt$cpmCutoff
cutoff_FDR   <- opt$qCutoff
cutoff_Pval  <- opt$PvalCutoff
cutoff_logFC <- opt$log2FCcutoff
topn         <- opt$topn

outDir <- dirname(opt$outPrefix)
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
outDir  <- normalizePath(outDir)
prefix0 <- basename(opt$outPrefix)
prefix  <- paste0(outDir, "/", prefix0, "_")   # global: used by SPADE ParallelDE2() workers

# SPADE's ErccQC() inside the per-comparison DE looks for ERCC_Controls_Analysis.txt
# in the working directory when no file is passed.
erccControlFile <- spade_ercc_control_file()
if (!file.exists("ERCC_Controls_Analysis.txt")) {
  file.copy(erccControlFile, "ERCC_Controls_Analysis.txt")
}

pipeline_version <- paste0("SPADE v", spade_version())
logFile <- paste0(prefix, "log.txt")
showMsgAndLog(logFile, msg = paste0("Run started: ", date()), append = FALSE)
showMsgAndLog(logFile, msg = paste0("R version: ", R.version$major, ".", R.version$minor), nStart = 1)
showMsgAndLog(logFile, msg = pipeline_version, nStart = 1)
showMsgAndLog(logFile, msg = sprintf("Exclude non-coding genes: %s", filterFlag), nStart = 1)

#--------------------------------------------------------------------------------
showMsgAndLog(logFile, msg = "[Parse AutoMapper ${SRM}.report file]", nStart = 2, timeStamp = "none")
if (!is.null(reportFile)) {
  outFile <- paste0(prefix, "MappingReport.tsv")
  mapStat <- ParseMappingReport(reportFile, metaFile, outFile = outFile)
  showMsgAndLog(logFile, msg = paste0(basename(outFile), "[saved]"), nStart = 1)
} else {
  mapStat <- NA
  showMsgAndLog(logFile, msg = "Mapping report not provided [skipped].", nStart = 1, tStart = 1)
}
objQC <- list(mapStat = mapStat)

#--------------------------------------------------------------------------------
showMsgAndLog(logFile, msg = "[read meta data and count file]", nStart = 2, timeStamp = "none")
datSet      <- ReadMetaAndCountsERCC(metaFile, countFile, filterFlag = FALSE)
counts      <- datSet[["counts"]]
counts_ERCC <- datSet[["counts_ERCC"]]
annotation  <- datSet[["annotation"]]
meta        <- datSet[["meta"]]

# ERCC detection: SPADE needs >= 5 ERCC- rows with reads in at least half of the samples.
ercc_vals  <- if (!is.null(counts_ERCC) && nrow(counts_ERCC) > 0) as.matrix(counts_ERCC[, meta$ID, drop = FALSE]) else NULL
n_ercc     <- if (is.null(ercc_vals)) 0 else nrow(ercc_vals)
ercc_ok    <- n_ercc >= 5 && sum(colSums(ercc_vals) == 0) <= nrow(meta) / 2
ercc_requested <- flag_ERCC
if (flag_ERCC && !ercc_ok) {
  flag_ERCC <- FALSE
  if (n_ercc == 0) {
    ercc_note <- "No ERCC- spike-in rows found in the count table; using library-size (TMM) normalization instead of ERCC normalization."
  } else {
    ercc_note <- sprintf("Found %d ERCC- rows but too few ERCC reads (need >= 5 rows and reads in at least half of the samples); using library-size (TMM) normalization instead of ERCC normalization.", n_ercc)
  }
} else if (flag_ERCC) {
  ercc_note <- sprintf("Found %d ERCC- spike-in rows; using ERCC normalization.", n_ercc)
} else if (n_ercc > 0) {
  ercc_note <- sprintf("Found %d ERCC- spike-in rows, but ERCC normalization was disabled (--ERCC FALSE); using library-size (TMM) normalization.", n_ercc)
} else {
  ercc_note <- "ERCC normalization disabled and no ERCC- rows present; using library-size (TMM) normalization."
}
normalization <- if (flag_ERCC) "ERCC" else "TMM"
showMsgAndLog(logFile, msg = paste0("Normalization: ", normalization, ". ", ercc_note), nStart = 1)
showMsgAndLog(logFile, msg = paste0("used meta data: ", nrow(meta), " x ", ncol(meta)), nStart = 1)
showMsgAndLog(logFile, msg = paste0("used count data: ", nrow(counts), " x ", ncol(counts)), nStart = 1)

#--------------------------------------------------------------------------------
showMsgAndLog(logFile, msg = "[check comparisons.lst]", nStart = 2)
if (is.null(opt$compList) || !nzchar(opt$compList)) {
  outFile <- paste0(prefix, "comparisons.lst")
  comps <- GetComparisons(meta, outFile = outFile)
  showMsgAndLog(logFile, msg = paste0(basename(outFile), "[saved]"), nStart = 1, tStart = 2)
} else {
  comps <- scan(file = normalizePath(opt$compList), character(), quote = "", quiet = TRUE)
}
comps <- gsub(" ", "", comps)
comps <- comps[nzchar(comps)]

#--------------------------------------------------------------------------------
msg_1 <- if (flag_ERCC) " based on ERCC normalization." else " based on library size normalization."
showMsgAndLog(logFile, msg = sprintf("[Standard RNAseq Differential Analysis] %s", msg_1), nStart = 2, timeStamp = "none")
if (flag_ERCC) {
  showMsgAndLog(logFile, msg = "step1. Perform ERCC QC", nStart = 2)
  VizErccConc(counts_ERCC, meta, paste0(prefix, "ERCC_QC"), ercc_control_file = erccControlFile)
  correlations <- ErccQC(counts_ERCC, meta, ercc_control_file = erccControlFile)
  ErccStat <- CalcERCCnormFactor(counts, counts_ERCC, meta, txtFile = NULL)
  # RNA_index as in RNAseqERCC v1.1 (endogenous RNA relative to ERCC)
  ErccStat$RNA_index <- ErccStat$gene_libSize / ErccStat$ERCC_libSize
  ErccStat$ERCC_COR  <- correlations[ErccStat$ID]
  spade_write_tsv(ErccStat, paste0(prefix, "ERCC_QC_stat.txt"))
  Ercc_CV <- as.data.frame(ErccCV(ErccStat))
  VizErccStat(ErccStat, outFile = paste0(prefix, "ERCC_QC_stat.pdf"))
  shift_test <- lapply(comps, function(x) ErccShiftTest(stat = ErccStat, comp = x, ratio = "pct_ERCC", group = "GROUP"))
  shift_test <- as.data.frame(do.call(rbind, lapply(shift_test, as.data.frame)))
  shift_test <- shift_test[, c(setdiff(colnames(shift_test), "SHIFT"), "SHIFT")]
  rownames(shift_test) <- shift_test$comp
  objQC <- list(mapStat = mapStat, ErccStat = ErccStat, ErccCV = Ercc_CV, ErccShift = shift_test)
} else {
  Ercc_CV <- NULL; ErccStat <- NULL; shift_test <- NULL
}
objQC$normalization <- normalization
objQC$ercc_note     <- ercc_note
objQC$n_ercc        <- n_ercc

counts_all <- counts   # before gene filtering, for the QC bundle
if (filterFlag) {  # refresh counts, meta and annotation with the non-coding filter
  datSet     <- ReadMetaAndCountsERCC(metaFile, countFile, filterFlag = filterFlag)
  counts     <- datSet[["counts"]]
  meta       <- datSet[["meta"]]
  annotation <- datSet[["annotation"]]
  counts_ERCC <- datSet[["counts_ERCC"]]
}

#--------------------------------------------------------------------------------
isFormula <- any(grepl(".*=\\(.*-.*\\)", comps, ignore.case = TRUE))
if (!isFormula) {
  invalid <- vapply(comps, function(x) {
    k <- unlist(strsplit(x, "_[Vv][Ss]_"))
    length(k) != 2 || k[1] == k[2]
  }, logical(1))
  for (x in comps[invalid]) showMsgAndLog(logFile, msg = paste0("**Error: Invalid comparison! ", x, " [skipped]**"), nStart = 2, tStart = 1, EOL = 1)
  comps <- comps[!invalid]
}
if (length(comps) == 0) stop("No valid comparison found.")
names(comps) <- paste0("C", seq_along(comps))
showMsgAndLog(logFile, msg = paste0("Found ", length(comps), " valid comparisons."), nStart = 1, tStart = 1)

#--------------------------------------------------------------------------------
showMsgAndLog(logFile, msg = "step.2  perform differential analyses..\n", nStart = 2)
objDE <- ParallelDE2(counts = counts, annotation = annotation, meta = meta, comps = comps,
  cutoff_count = cutoff_count, cutoff_cpm = cutoff_cpm, logFile = logFile,
  ncores = max(1, opt$ncores), counts_ERCC = if (flag_ERCC) counts_ERCC else NULL)
if (!is.list(objDE) || !all(names(comps) %in% names(objDE))) stop("Differential analysis failed.")

statSizes <- NULL; statTopn <- NULL
numDEGs <- list(); topnDEGs <- list()
sizes <- table(meta$GROUP)
for (n in seq_along(comps)) {
  objStr     <- names(comps)[n]
  diff_table <- objDE[[objStr]]$DE
  meta1      <- objDE[[objStr]]$meta
  comp <- comps[n]
  if (!isFormula) {
    comp <- gsub("_[vV][Ss]_", "_vs_", comp)
    keywords <- unlist(strsplit(comp, "_vs_"))
    sizes1 <- c(sum(toupper(meta1$GROUP) == toupper(keywords[1])), sum(toupper(meta1$GROUP) == toupper(keywords[2])))
    strSizes <- paste0(sizes1, collapse = " vs ")
  } else {
    allGroups <- ExtractGroups(comp)
    strSizes <- paste(sizes[allGroups], collapse = " , ")
  }
  statSizes <- c(statSizes, strSizes)
  resDEGs <- GetDEGs(diff_table, cutoff_logFC = cutoff_logFC, cutoff_Pval = cutoff_Pval, cutoff_FDR = cutoff_FDR, topn = topn)
  numDEGs[[objStr]]  <- resDEGs$numDEGs
  topnDEGs[[objStr]] <- resDEGs$topnDEGs
  statTopn <- c(statTopn, paste0(length(resDEGs$topnDEGs$topUpGenes), " up, ", length(resDEGs$topnDEGs$topDownGenes), " down"))
}
statDEGs <- data.frame(ID = names(comps),
  comparison = if (isFormula) comps else gsub("_vs_", " vs ", comps, ignore.case = TRUE),
  sample_size = statSizes, top_DEGs = statTopn, row.names = NULL)
# DEG count tables are written by spade_report.R with the report cutoffs.

dat <- if (tolower(colnames(counts)[1]) == "gene") counts[, -1] else counts
log2CPM <- cpm(dat, log = TRUE)

settings <- list(countCutoff = cutoff_count, cpmCutoff = cutoff_cpm, filter = filterFlag,
  ERCC_requested = ercc_requested, normalization = normalization,
  log2FCcutoff = cutoff_logFC, PvalCutoff = cutoff_Pval, qCutoff = cutoff_FDR, topn = topn,
  pipeline = pipeline_version, countFile = basename(countFile), metaFile = basename(metaFile))
objDEs <- list(QC = objQC, data = datSet, log2CPM = log2CPM, comps = comps, statDEGs = statDEGs,
  numDEGs = numDEGs, topnDEGs = topnDEGs, allDE = objDE, annotation = annotation, settings = settings)
rdsFile <- paste0(prefix, "SPADE_Results.rds")
saveRDS(objDEs, file = rdsFile)
showMsgAndLog(logFile, msg = paste0(basename(rdsFile), "[saved]"), nStart = 1, tStart = 1)

#--------------------------------------------------------------------------------
# Extras after the results object is saved: an error here is logged, not fatal.
showMsgAndLog(logFile, msg = "Saving rnk files and GSEA inputs...", nStart = 2, tStart = 1)
for (n in seq_along(comps)) {
  tryCatch({
    rnkInput <- objDE[[names(comps)[n]]]$DE
    pv <- as.numeric(rnkInput$P.Value)
    pv[!is.finite(pv) | pv <= 0] <- .Machine$double.xmin
    rnkInput$P.Value <- pv
    compName <- if (isFormula) unlist(strsplit(comps[n], "="))[1] else comps[n]
    GenerateRank4GSEA(rnkInput, annotation = annotation, prefix = paste0(prefix, compName))
  }, error = function(e) showMsgAndLog(logFile, msg = paste0("Warning: rnk files failed for ", comps[n], ": ", conditionMessage(e)), nStart = 1))
}
tryCatch(GenerateTMMnCLS4GSEA(counts, annotation = annotation, meta = meta, prefix = paste0(prefix, "all")),
  error = function(e) showMsgAndLog(logFile, msg = paste0("Warning: GSEA TMM/cls files failed: ", conditionMessage(e)), nStart = 1))

showMsgAndLog(logFile, msg = "Saving QC bundle...", nStart = 2, tStart = 1)
tryCatch(saveStandardQcBundle(counts = counts, logCPM = log2CPM, meta = meta, outPrefix = paste0(outDir, "/", prefix0)),
  error = function(e) showMsgAndLog(logFile, msg = paste0("Warning: QC bundle failed: ", conditionMessage(e)), nStart = 1))
if (flag_ERCC) {
  tryCatch({
    erccLog2CPM <- calcErccQcLog2CPM(counts, counts_ERCC, meta)
    saveErccQcBundle(stat = ErccStat, cv = Ercc_CV, shift = shift_test, tmmLogCPM = log2CPM,
      erccLogCPM = erccLog2CPM, meta = meta, outPrefix = paste0(outDir, "/", prefix0))
  }, error = function(e) showMsgAndLog(logFile, msg = paste0("Warning: ERCC QC bundle failed: ", conditionMessage(e)), nStart = 1))
}

showMsgAndLog(logFile, msg = "Done.", nStart = 2)
cat("\n\n")
invisible(gc())
