#!/usr/bin/env Rscript
# spade_report.R - SPADE-nf HTML report.
#
# Follows RNAseqERCC v1.1 stdRNAseq.ERCC.Reporter.R: reads the results object written by
# spade_de.R, builds PCA / volcano / MA / heatmap plots per comparison and renders the Rmd.
# Additions for SPADE: library-size and RLE QC plots, fgsea and Enrichr summaries, and two
# brandings: --branding spade (default, external users) or cab (CAB internal, as RNAseqERCC).
# Report cutoffs come from the config file (CUTOFF_*, NUM_GENE_PCA, NUM_DEG_TOPN), so changing
# them only re-runs this step.

suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = NA, help = "config file (KEY: value lines) with project info and report cutoffs"),
  make_option(c("-i", "--input"), type = "character", default = NA, help = "results object (*_SPADE_Results.rds) from spade_de.R"),
  make_option(c("-o", "--outDir"), type = "character", default = ".", help = "output directory [default %default]"),
  make_option(c("-V", "--volcano"), type = "character", default = "FALSE", help = "also save one volcano plot PDF per comparison [default %default]"),
  make_option("--enrichment_dir", type = "character", default = NULL, help = "[optional] directory with *_fgsea.tsv and *_enrichr.tsv files"),
  make_option("--branding", type = "character", default = "spade", help = "spade or cab [default %default]"),
  make_option("--assets_dir", type = "character", default = NULL, help = "directory with the Rmd, bib and logos [default: ../assets/report next to this script]"),
  make_option(c("-r", "--Rmd"), type = "character", default = NULL, help = "[optional] use your own Rmd"),
  make_option("--spade_dir", type = "character", default = NULL, help = "SPADE checkout [default: $SPADE_HOME or /opt/spade]")
)
opt <- parse_args(OptionParser(usage = "%prog -c <config.yaml> -i <prefix_SPADE_Results.rds> [-o outDir --enrichment_dir dir --branding spade|cab]",
  option_list = option_list))
if (is.na(opt$config) || is.na(opt$input)) stop("--config and --input are required")

scriptDir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(scriptDir, "spade_cli_helpers.R"))
spade_load(opt$spade_dir)
suppressPackageStartupMessages({
  library(pheatmap); library(knitr); library(plotly); library(grid); library(gridExtra); library(dplyr)
})
options(bitmapType = "cairo")

branding <- tolower(opt$branding)
if (!branding %in% c("spade", "cab")) stop("--branding must be spade or cab")
assetsDir <- if (!is.null(opt$assets_dir)) opt$assets_dir else file.path(dirname(scriptDir), "assets", "report")
assetsDir <- normalizePath(assetsDir, winslash = "/", mustWork = TRUE)
RmdFile <- if (!is.null(opt$Rmd)) normalizePath(opt$Rmd, mustWork = TRUE) else file.path(assetsDir, "SPADE_DE_Report.Rmd")
bibFile <- file.path(assetsDir, "SPADE_report.bib")
logoFile <- file.path(assetsDir, if (branding == "cab") "CAB_logo.png" else "SPADE_logo.png")
for (f in c(RmdFile, bibFile, logoFile)) if (!file.exists(f)) stop("Report asset not found: ", f)

yaml <- ReadYaml(opt$config)
num <- function(key, default) { v <- suppressWarnings(as.numeric(yaml[[key]])); if (length(v) == 0 || is.na(v)) default else v }
cutoff_FDR   <- num("CUTOFF_FDR", 0.05)
cutoff_Pval  <- num("CUTOFF_PVALUE", 0.05)
cutoff_logFC <- num("CUTOFF_LOG2FC", 1)
topn_pca     <- num("NUM_GENE_PCA", 3000)
topn_deg     <- num("NUM_DEG_TOPN", 10)
volcanoFlag  <- spade_flag(opt$volcano)

outDir <- opt$outDir
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
outDir <- normalizePath(outDir)
if (is.null(yaml$OUTPUT_DIR) || !nzchar(yaml$OUTPUT_DIR)) yaml$OUTPUT_DIR <- outDir
if (is.null(yaml$SRM) || !nzchar(yaml$SRM)) yaml$SRM <- "SPADE"
yaml$BRANDING <- branding
logFile <- file.path(outDir, paste0(yaml$SRM, "_report_log.txt"))
showMsgAndLog(logFile, msg = paste0("Report started: ", date()), append = FALSE)

objDEs <- readRDS(normalizePath(opt$input, mustWork = TRUE))
if (!all(c("QC", "log2CPM", "data", "comps", "statDEGs", "allDE") %in% names(objDEs))) stop("Invalid SPADE results object: ", opt$input)
yaml$PIPELINE_RUN <- objDEs$settings$pipeline
if (is.null(yaml$PIPELINE) || !nzchar(yaml$PIPELINE)) yaml$PIPELINE <- objDEs$settings$pipeline
# numDEGs are recomputed with the report cutoffs so the report and its tables agree
log2CPM <- objDEs[["log2CPM"]]
meta    <- objDEs[["data"]][["meta"]]
comps   <- objDEs[["comps"]]

#--------------------------------------------------------------------------------
showMsgAndLog(logFile, msg = "[Preparing plots]...", nStart = 2, timeStamp = "none")
colors4DEGs <- c(fill_up = "#F9ABAB", fill_down = "#ADADF9", color_up = "#F92A2A", color_down = "#2727FA", GOI = "purple")
pdf(NULL)  # avoid Rplots.pdf
pcaPlot0 <- PlotPCA(as.data.frame(log2CPM), meta, scale = TRUE, topn = topn_pca, outFile = NA, label = TRUE)
libSizePlot <- tryCatch(PlotLibSize(counts = objDEs$data$counts, meta = meta), error = function(e) NULL)
rlePlot <- tryCatch(PlotRLE(logCPM = log2CPM, meta = meta), error = function(e) NULL)
if (volcanoFlag) { vpDir <- file.path(outDir, "volcanoPlots"); dir.create(vpDir, showWarnings = FALSE) }

pcaPlots <- list(); volcanoPlots <- list(); maPlots <- list(); heatmaps <- list(); numDEGs <- list(); topnDEGs <- list()
for (n in seq_along(comps)) {
  objStr <- names(comps)[n]
  showMsgAndLog(logFile, msg = paste(n, comps[n]), nStart = 1)
  objDE1 <- objDEs[["allDE"]][[objStr]]
  diff_table <- objDE1$DE
  meta1 <- objDE1$meta
  title <- gsub("_vs_", " vs ", comps[objStr], ignore.case = TRUE)
  pcaPlots[[objStr]] <- PlotPCA(as.data.frame(log2CPM), meta1, scale = TRUE, topn = topn_pca, outFile = NA, label = TRUE, title = title)
  resDEGs <- GetDEGs(diff_table, cutoff_logFC = cutoff_logFC, cutoff_Pval = cutoff_Pval, cutoff_FDR = cutoff_FDR, topn = topn_deg)
  numDEGs[[objStr]] <- resDEGs$numDEGs
  usedDEGs <- resDEGs$topnDEGs
  topnDEGs[[objStr]] <- usedDEGs
  # P.Value of exactly 0 (global shifts) -> smallest non-zero P.Value, as in RNAseqERCC
  pv <- diff_table$P.Value
  if (any(pv == 0, na.rm = TRUE)) diff_table$P.Value[pv == 0] <- min(pv[pv > 0], na.rm = TRUE)
  outFileV <- if (volcanoFlag) file.path(vpDir, paste0("volcanoPlot", sprintf("%03d", n), "_", comps[objStr], ".pdf")) else NA
  volcanoPlots[[objStr]] <- PlotVolcano(diff_table, cutoff_logFC = cutoff_logFC, cutoff_FDR = cutoff_FDR,
    cols = NA, topnDEGs = usedDEGs, outFile = outFileV)  # SPADE default colors
  maPlots[[objStr]] <- PlotMA(diff_table, cutoff_logFC = cutoff_logFC, cols = NA, topnDEGs = usedDEGs, outFile = NA)
  genes <- unlist(usedDEGs)
  heatmaps[objStr] <- list(if (length(genes) > 0) tryCatch(PlotHeatmap(diff_table, meta1, outFile = NA, title = title, genes = genes),
    error = function(e) NULL) else NULL)
}
objPlots <- list(pcaPlot0 = pcaPlot0, pcaPlots = pcaPlots, volcanoPlots = volcanoPlots, heatmaps = heatmaps,
  maPlots = maPlots, libSizePlot = libSizePlot, rlePlot = rlePlot)
objDEs$numDEGs <- numDEGs
objDEs$statDEGs$top_DEGs <- vapply(names(comps), function(x)
  paste0(length(topnDEGs[[x]]$topUpGenes), " up, ", length(topnDEGs[[x]]$topDownGenes), " down"), character(1))
spade_write_tsv(objDEs$statDEGs, file.path(outDir, paste0(yaml$SRM, "_DEG_summary.tsv")))
numTable <- do.call(rbind, lapply(names(numDEGs), function(x) cbind(ID = x, comparison = comps[[x]], numDEGs[[x]])))
spade_write_tsv(numTable, file.path(outDir, paste0(yaml$SRM, "_DEG_counts.tsv")))

#--------------------------------------------------------------------------------
# Enrichment summaries (optional)
readAll <- function(pattern) {
  if (is.null(opt$enrichment_dir) || !dir.exists(opt$enrichment_dir)) return(NULL)
  fs <- list.files(opt$enrichment_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(fs) == 0) return(NULL)
  tabs <- lapply(fs, function(f) tryCatch(read.delim(f, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL))
  tabs <- Filter(function(x) !is.null(x) && nrow(x) > 0, tabs)
  if (length(tabs) == 0) return(data.frame())
  dplyr::bind_rows(tabs)
}
fgseaTab <- readAll("_fgsea\\.tsv$")
enrichrTab <- readAll("_enrichr\\.tsv$")
if (!is.null(fgseaTab) && nrow(fgseaTab) > 0) {
  for (cmp in unique(fgseaTab$comparison)) {
    wb <- openxlsx::createWorkbook()
    sub <- fgseaTab[fgseaTab$comparison == cmp, ]
    for (gs in unique(sub$gene_set)) {
      sheet <- substr(gsub("[\\[\\]:*?/\\\\]", "_", gs, perl = TRUE), 1, 31)
      openxlsx::addWorksheet(wb, sheet)
      openxlsx::writeData(wb, sheet, sub[sub$gene_set == gs, ])
    }
    f <- file.path(outDir, paste0(yaml$SRM, "_", cmp, "_fgsea_results.xlsx"))
    openxlsx::saveWorkbook(wb, f, overwrite = TRUE)
    showMsgAndLog(logFile, msg = paste0(basename(f), "[saved]"), nStart = 1, tStart = 1)
  }
}
enrichment <- list(fgsea = fgseaTab, enrichr = enrichrTab, top_n = topn_deg)

#--------------------------------------------------------------------------------
showMsgAndLog(logFile, msg = "[Generating report]...", nStart = 2, timeStamp = "none")
while (!is.null(dev.list())) invisible(dev.off())
workDir <- file.path(outDir, ".report_tmp")
dir.create(workDir, showWarnings = FALSE)
invisible(file.copy(RmdFile, file.path(workDir, "SPADE_DE_Report.Rmd"), overwrite = TRUE))
rmarkdown::render(file.path(workDir, "SPADE_DE_Report.Rmd"),
  params = list(objDEs = objDEs, yaml = yaml, objPlots = objPlots, bib_path = bibFile, logo = logoFile, enrichment = enrichment),
  output_format = "html_document", output_dir = outDir, intermediates_dir = workDir, quiet = TRUE)
unlink(workDir, recursive = TRUE)
htmlFile <- file.path(outDir, "SPADE_DE_Report.html")
htmlNew <- file.path(outDir, if (branding == "cab") paste0("stdRNAseqERCC_DE_Report_", yaml$SRM, ".html") else paste0("SPADE_DE_Report_", yaml$SRM, ".html"))
if (!file.exists(htmlFile) || !file.rename(htmlFile, htmlNew)) stop("Report rendering failed: ", htmlFile)
showMsgAndLog(logFile, msg = paste0(basename(htmlNew), "[saved]"), nStart = 1, tStart = 1)
