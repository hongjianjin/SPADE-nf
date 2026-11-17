#!/usr/bin/env Rscript
# spade_enrichr.R - SPADE-nf Enrichr over-representation analysis for one comparison.
#
# Needs outbound HTTPS to maayanlab.cloud (Enrichr web API). Gene selection follows the SPADE
# app: genes ranked by --rank_metric (t, or sign(log2FC)*-log10(P)); top N up (log2FC > 0) and
# top N down (log2FC < 0) are submitted separately. Databases: --dbs standard picks the newest
# Enrichr library for Hallmark, GO BP/MF/CC, KEGG and Reactome (SPADE's "Standard" picker);
# or give a comma-separated list of Enrichr library names.
# Terms with P < --pval_cutoff are kept; SPADE bar and dot plots of the top terms are saved.

suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option("--de", type = "character", default = NA, help = "SPADE DE table (*_log2CPM.tsv)"),
  make_option("--comparison", type = "character", default = NULL, help = "comparison label"),
  make_option("--dbs", type = "character", default = "standard", help = "'standard' or comma-separated Enrichr libraries [default %default]"),
  make_option("--species", type = "character", default = "human", help = "human or mouse (picks KEGG library) [default %default]"),
  make_option("--rank_metric", type = "character", default = "t", help = "t or logPval [default %default]"),
  make_option("--top_genes", type = "integer", default = 250, help = "genes per direction: 100, 250 or 500 [default %default]"),
  make_option("--pval_cutoff", type = "numeric", default = 0.05, help = "keep terms with P < cutoff [default %default]"),
  make_option("--top_n", type = "integer", default = 10, help = "terms per plot (5-25) [default %default]"),
  make_option(c("-o", "--outPrefix"), type = "character", default = NA, help = "outDir/prefix"),
  make_option("--spade_dir", type = "character", default = NULL, help = "SPADE checkout [default: $SPADE_HOME or /opt/spade]")
)
opt <- parse_args(OptionParser(usage = "%prog --de <comp_log2CPM.tsv> -o <prefix> [--dbs standard]", option_list = option_list))
if (is.na(opt$de) || is.na(opt$outPrefix)) stop("--de and --outPrefix are required")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "spade_cli_helpers.R"))
spade_load(opt$spade_dir, enrichment = TRUE)

comp   <- if (!is.null(opt$comparison)) opt$comparison else sub("_log2CPM\\.tsv$", "", basename(opt$de))
prefix <- opt$outPrefix
dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
topN   <- max(5L, min(25L, opt$top_n))

# --- databases ----------------------------------------------------------------
suppressPackageStartupMessages(library(enrichR))  # attaching sets the Enrichr server options
invisible(tryCatch(enrichR::setEnrichrSite("Enrichr"), error = function(e) NULL))
allDbs <- tryCatch(enrichR::listEnrichrDbs(), error = function(e) NULL)
if (is.null(allDbs) || !is.data.frame(allDbs) || nrow(allDbs) == 0) {
  stop("Cannot reach the Enrichr server (needs internet access). Run without --run_enrichr or on a node with internet.")
}
libs <- as.character(allDbs$libraryName)
pickLatest <- function(pattern) {
  hits <- grep(pattern, libs, value = TRUE, ignore.case = TRUE)
  if (grepl("KEGG", pattern)) {
    sp <- grep(if (tolower(opt$species) == "mouse") "Mouse" else "Human", hits, value = TRUE)
    if (length(sp) > 0) hits <- sp
  }
  if (length(hits) == 0) return(NULL)
  yr <- suppressWarnings(as.integer(sub(".*_(\\d{4}).*", "\\1", hits)))
  hits[order(yr, decreasing = TRUE, na.last = TRUE)][1]
}
if (tolower(opt$dbs) == "standard") {
  pats <- c("^MSigDB_Hallmark", "^GO_Biological_Process", "^GO_Molecular_Function", "^GO_Cellular_Component", "^KEGG_", "^Reactome")
  dbs <- unique(unlist(lapply(pats, pickLatest)))
} else {
  dbs <- trimws(unlist(strsplit(opt$dbs, ",")))
  unknown <- setdiff(dbs, libs)
  if (length(unknown) > 0) stop("Unknown Enrichr libraries: ", paste(unknown, collapse = ", "))
}
cat("Enrichr libraries:", paste(dbs, collapse = ", "), "\n")

# --- genes (SPADE app selection) ----------------------------------------------
de <- read.delim(opt$de, check.names = FALSE, stringsAsFactors = FALSE)
genes <- as.character(de[[intersect(c("geneSymbol", "gene", "symbol"), colnames(de))[1]]])
lfc <- as.numeric(de[["log2FC"]])
score <- if (opt$rank_metric == "t") as.numeric(de[["t.Statistic"]]) else {
  p <- as.numeric(de[["P.Value"]]); p[!is.finite(p) | p <= 0] <- .Machine$double.xmin
  sign(lfc) * -log10(p)
}
topGenes <- if (opt$top_genes %in% c(100L, 250L, 500L)) opt$top_genes else 250L
upIdx <- which(is.finite(lfc) & lfc > 0 & is.finite(score)); upIdx <- upIdx[order(score[upIdx], decreasing = TRUE)]
dnIdx <- which(is.finite(lfc) & lfc < 0 & is.finite(score)); dnIdx <- dnIdx[order(score[dnIdx], decreasing = FALSE)]
clean <- function(x) { x <- unique(x); x[!is.na(x) & nzchar(x)] }
genesUp <- clean(genes[head(upIdx, topGenes)])
genesDn <- clean(genes[head(dnIdx, topGenes)])
cat("Submitting", length(genesUp), "up and", length(genesDn), "down genes\n")
writeLines(genesUp, paste0(prefix, "_enrichr_input_up.txt"))
writeLines(genesDn, paste0(prefix, "_enrichr_input_down.txt"))

res <- run_enrichr(de, selected_dbs = dbs, enrichment_pval_cutoff = opt$pval_cutoff,
  top_up = topGenes, top_down = topGenes, genes_up = genesUp, genes_down = genesDn,
  file_prefix = basename(prefix), comp_label = comp, out_dir = dirname(prefix))

# --- tables and plots ---------------------------------------------------------
rows <- list()
for (key in names(res)) {
  e <- res[[key]]
  for (dir in c("up", "down")) {
    t <- e[[paste0(dir, "_result")]]
    if (!is.null(t) && is.data.frame(t) && nrow(t) > 0) {
      rows[[length(rows) + 1]] <- cbind(comparison = comp, direction = dir, library = e$db_name, t)
    }
  }
  stub <- paste0(prefix, "_", gsub("[^[:alnum:]]+", "_", e$db_name))
  for (pt in c("bar", "dot")) {
    p <- tryCatch(build_enrichr_plot_bundle(e, plot_type = pt, top_up = topN, top_down = topN), error = function(err) NULL)
    if (is.null(p)) next
    for (dir in c("up", "down")) {
      pp <- p[[dir]]
      if (!is.null(pp) && inherits(pp, "ggplot") && !isTRUE(attr(pp, "shinyRNAseq_empty_plot"))) {
        f <- paste0(stub, "_", pt, "_", dir, ".pdf")
        ggplot2::ggsave(f, pp, width = 10, height = max(5, ceiling(topN * 0.4) + 2), units = "in")
        cat(basename(f), "[saved]\n")
      }
    }
  }
}
tab <- if (length(rows) > 0) do.call(rbind, rows) else
  data.frame(comparison = comp, direction = NA, library = NA, Term = NA, Message = "No terms passed the P-value cutoff.")
spade_write_tsv(tab, paste0(prefix, "_enrichr.tsv"))
