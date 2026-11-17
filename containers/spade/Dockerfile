# Container for the SPADE-nf pipeline: R + Bioconductor + the SPADE analysis functions
#
# Built by .github/workflows/build-container.yml and pushed to
#   ghcr.io/hongjianjin/spade:<SPADE version>
#
# Manual build:
#   docker build --build-arg SPADE_REF=3.0.1 -t ghcr.io/hongjianjin/spade:3.0.1 containers/spade
#
# Singularity on the St. Jude HPC:
#   singularity pull spade_3.0.1.sif docker://ghcr.io/hongjianjin/spade:3.0.1
#
# SPADE (https://github.com/hongjianjin/SPADE) is a Shiny app, not an R package: its analysis
# functions (module/RNAseq*.R) and ERCC concentration table are copied to /opt/spade, where the
# pipeline's bin/ scripts source them. The Shiny UI packages are not installed.

# R 4.3 / Bioconductor 3.18, the environment SPADE 3.0 was developed with
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# git tag (or commit) of hongjianjin/SPADE to copy
ARG SPADE_REF=3.0.1

LABEL org.opencontainers.image.title="spade" \
      org.opencontainers.image.description="R environment with the SPADE analysis functions for the SPADE-nf Nextflow pipeline" \
      org.opencontainers.image.source="https://github.com/hongjianjin/SPADE-nf" \
      org.opencontainers.image.licenses="Apache-2.0"

# Use only the container's R libraries. Without this, an R_LIBS/R_LIBS_USER exported on the host
# (e.g. a personal R library in $HOME) leaks in through Singularity and R loads packages
# compiled for another R version ("libRlapack.so: cannot open shared object file").
ENV R_LIBS=/usr/local/lib/R/site-library:/usr/local/lib/R/library \
    R_LIBS_USER=/dev/null \
    R_LIBS_SITE=/usr/local/lib/R/site-library \
    SPADE_HOME=/opt/spade

# procps: `ps` for Nextflow task metrics; pandoc: HTML report
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps pandoc git \
    && rm -rf /var/lib/apt/lists/*

# CRAN packages come from a dated Posit Package Manager snapshot (Ubuntu 22.04 binaries) from the
# Bioconductor 3.18 period: current CRAN has packages that need R >= 4.4 (e.g. Deriv, needed by
# ggpubr via car/rstatix), and install.packages() only warns when one fails.
ARG CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/__linux__/jammy/2024-04-24
ENV SPADE_R_PKGS="limma edgeR fgsea enrichplot optparse ggplot2 ggpubr ggrepel ggtext patchwork reshape2 pheatmap plotly dplyr tidyr tibble stringr tidyverse scales data.table openxlsx gridExtra htmltools knitr rmarkdown formattable xfun enrichR msigdbr"

RUN Rscript -e "options(repos = c(CRAN = '${CRAN_SNAPSHOT}')); \
        BiocManager::install(c('limma', 'edgeR', 'fgsea', 'enrichplot'), ask = FALSE, update = FALSE)" \
    && Rscript -e "pkgs <- strsplit(Sys.getenv('SPADE_R_PKGS'), ' ')[[1]]; \
        install.packages(setdiff(pkgs, rownames(installed.packages())), repos = '${CRAN_SNAPSHOT}'); \
        bad <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; \
        if (length(bad)) stop('R packages not installed: ', paste(bad, collapse = ', '))"

RUN git clone https://github.com/hongjianjin/SPADE.git /tmp/SPADE \
    && cd /tmp/SPADE && git checkout "${SPADE_REF}" \
    && mkdir -p /opt/spade/module /opt/spade/data \
    && cp module/RNAseq.Funcs.R module/RNAseq.ERCC.Func.R module/RNAseq.fgsea.Func.R module/RNAseq.enrichr.Func.R /opt/spade/module/ \
    && cp data/ERCC_Controls_Analysis.txt /opt/spade/data/ \
    && cp DESCRIPTION LICENSE /opt/spade/ \
    && echo "${SPADE_REF#v}" > /opt/spade/VERSION \
    && git rev-parse HEAD > /opt/spade/COMMIT \
    && rm -rf /tmp/SPADE

# Fail the build if any package does not attach or any SPADE function file cannot be sourced
RUN Rscript -e 'for (p in strsplit(Sys.getenv("SPADE_R_PKGS"), " ")[[1]]) suppressPackageStartupMessages(library(p, character.only = TRUE)); \
        for (f in list.files("/opt/spade/module", full.names = TRUE)) sys.source(f, envir = globalenv()); \
        stopifnot(exists("ParallelDE2"), exists("run_fgsea"), exists("run_enrichr")); \
        stopifnot(rmarkdown::pandoc_available()); cat("SPADE", readLines("/opt/spade/VERSION"), "OK\n")'

CMD ["R"]
