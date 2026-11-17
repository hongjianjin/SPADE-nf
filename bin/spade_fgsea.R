#!/usr/bin/env Rscript
# spade_fgsea.R - SPADE-nf gene set enrichment for one comparison and one gene set (GMT).
#
# Uses SPADE's run_fgsea() with its custom-GMT path (fgsea::gmtPathways), so the MSigDB
# release is whatever GMT is given (the pipeline downloads v2026.2 Hs/Mm symbol GMTs).
# Ranks come from the SPADE DE table (*_log2CPM.tsv): gene symbols, metric
#   t       : moderated t statistic (t.Statistic)
#   logPval : sign(log2FC) * -log10(P.Value)
#   log2FC  : log2 fold change
# Duplicated symbols keep the value with the largest |rank|. Results with P < --pval_cutoff
# are kept (SPADE default 0.05); plots follow SPADE's bar and dot plots (top N terms).

suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option("--de", type = "character", default = NA, help = "SPADE DE table (*_log2CPM.tsv) for one comparison"),
  make_option("--gmt", type = "character", default = NA, help = "gene set GMT file (gene symbols)"),
  make_option("--gene_set", type = "character", default = NULL, help = "label for the gene set [default: GMT file name]"),
  make_option("--comparison", type = "character", default = NULL, help = "comparison label [default: from --de file name]"),
  make_option("--rank_metric", type = "character", default = "t", help = "t, logPval or log2FC [default %default]"),
  make_option("--pval_cutoff", type = "numeric", default = 0.05, help = "keep pathways with pval < cutoff [default %default]"),
  make_option("--min_size", type = "integer", default = 15, help = "minimum gene set size [default %default]"),
  make_option("--top_n", type = "integer", default = 10, help = "pathways per plot (5-25) [default %default]"),
  make_option("--seed", type = "integer", default = 42, help = "random seed for fgsea [default %default]"),
  make_option(c("-o", "--outPrefix"), type = "character", default = NA, help = "outDir/prefix"),
  make_option("--spade_dir", type = "character", default = NULL, help = "SPADE checkout [default: $SPADE_HOME or /opt/spade]")
)
opt <- parse_args(OptionParser(usage = "%prog --de <comp_log2CPM.tsv> --gmt <set.gmt> -o <prefix>", option_list = option_list))
if (is.na(opt$de) || is.na(opt$gmt) || is.na(opt$outPrefix)) stop("--de, --gmt and --outPrefix are required")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "spade_cli_helpers.R"))
spade_load(opt$spade_dir, enrichment = TRUE)

geneSet <- if (!is.null(opt$gene_set)) opt$gene_set else sub("\\.gmt$", "", basename(opt$gmt))
comp    <- if (!is.null(opt$comparison)) opt$comparison else sub("_log2CPM\\.tsv$", "", basename(opt$de))
topN    <- max(5L, min(25L, opt$top_n))
prefix  <- opt$outPrefix
dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)

de <- read.delim(opt$de, check.names = FALSE, stringsAsFactors = FALSE)
symCol <- intersect(c("geneSymbol", "gene", "symbol"), colnames(de))[1]
if (is.na(symCol)) stop("DE table has no geneSymbol/gene column")
genes <- sub("_ENS.*$", "", as.character(de[[symCol]]))  # SPADE suffixes duplicated symbols with _<geneID>
metric <- switch(opt$rank_metric,
  t       = as.numeric(de[["t.Statistic"]]),
  logPval = { p <- as.numeric(de[["P.Value"]]); p[!is.finite(p) | p <= 0] <- .Machine$double.xmin
              sign(as.numeric(de[["log2FC"]])) * -log10(p) },
  log2FC  = as.numeric(de[["log2FC"]]),
  stop("--rank_metric must be t, logPval or log2FC"))
ok <- is.finite(metric) & !is.na(genes) & nzchar(genes)
genes <- genes[ok]; metric <- metric[ok]
o <- order(abs(metric), decreasing = TRUE)
genes <- genes[o]; metric <- metric[o]
keep <- !duplicated(genes)
ranks <- sort(setNames(metric[keep], genes[keep]), decreasing = TRUE)
cat("Ranks:", length(ranks), "genes; metric", opt$rank_metric, "; gene set", geneSet, "\n")

set.seed(opt$seed)
res <- run_fgsea(ranks = ranks, custom_gmt = opt$gmt, pval_cutoff = opt$pval_cutoff,
                 sname = paste(comp, geneSet), run_only = TRUE)
res <- as.data.frame(res)
if (nrow(res) > 0 && "leadingEdge" %in% colnames(res)) {
  res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
}
res <- cbind(comparison = rep(comp, nrow(res)), gene_set = rep(geneSet, nrow(res)), res)
spade_write_tsv(res, paste0(prefix, "_fgsea.tsv"))
cat("Pathways with pval <", opt$pval_cutoff, ":", nrow(res), "\n")

# Plots (SPADE styles). run_fgsea() expects the list column back for plotting.
saveIfData <- function(p, file, height) {
  if (is.null(p) || isTRUE(attr(p, "shinyRNAseq_empty_plot"))) return(invisible(FALSE))
  ggplot2::ggsave(file, plot = p, width = 10, height = height, units = "in")
  cat(basename(file), "[saved]\n")
}
if (nrow(res) > 0) {
  plotTab <- tibble::as_tibble(res)
  plotTab$leadingEdge <- strsplit(plotTab$leadingEdge, ";", fixed = TRUE)
  h <- max(5, ceiling(topN * 0.4) + 2)
  tryCatch({
    pb <- run_fgsea(fgsea_res_tab = plotTab, pval_cutoff = opt$pval_cutoff, top_n = topN, plot_type = "bar",
                    sname = paste0(comp, ": ", geneSet))
    saveIfData(pb, paste0(prefix, "_fgsea_bar.pdf"), h)
    pd <- run_fgsea(fgsea_res_tab = plotTab, pval_cutoff = opt$pval_cutoff, top_n = topN, plot_type = "dot",
                    sname = paste0(comp, ": ", geneSet))
    saveIfData(pd$up, paste0(prefix, "_fgsea_dot_up.pdf"), h)
    saveIfData(pd$down, paste0(prefix, "_fgsea_dot_down.pdf"), h)
  }, error = function(e) cat("Warning: fgsea plots failed:", conditionMessage(e), "\n"))
}
