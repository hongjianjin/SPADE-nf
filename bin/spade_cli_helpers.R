# spade_cli_helpers.R
# Shared code for the SPADE-nf command-line scripts in bin/.
#
#   1. spade_load(): sources the SPADE v3 analysis functions (module/*.R of
#      https://github.com/hongjianjin/SPADE), which the container installs in
#      /opt/spade. Override with --spade_dir or the SPADE_HOME variable.
#   2. Helpers the CLI needs that SPADE does not provide, ported unchanged from
#      RNAseqERCC v1.1 (stdRNAseq.Func.R, https://github.com/hongjianjin/RNAseqERCC),
#      except ParseMappingReport(), adapted to SPADE's ReadMeta() return value.
#   3. Small CLI utilities (spade_*).

spade_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) == 0) return(getwd())
  dirname(normalizePath(f[1]))
}

spade_find_home <- function(spade_dir = NULL) {
  cands <- c(spade_dir, Sys.getenv("SPADE_HOME", unset = ""), "/opt/spade")
  cands <- cands[!is.na(cands) & nzchar(cands)]
  for (d in cands) {
    if (file.exists(file.path(d, "module", "RNAseq.Funcs.R"))) return(normalizePath(d))
  }
  stop("SPADE function files not found. Looked in: ", paste(cands, collapse = ", "),
       ". Set --spade_dir or SPADE_HOME to a SPADE checkout.")
}

# Source SPADE's analysis functions (no Shiny needed) into the global environment.
spade_load <- function(spade_dir = NULL, enrichment = FALSE) {
  home <- spade_find_home(spade_dir)
  # Packages SPADE's functions expect on the search path (the Shiny app attaches them in app.R).
  for (p in c("stringr", "limma", "edgeR", "ggplot2", "ggpubr", "ggrepel", "reshape2", "patchwork", "dplyr", "parallel")) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  files <- c("RNAseq.Funcs.R", "RNAseq.ERCC.Func.R")
  if (enrichment) files <- c(files, "RNAseq.fgsea.Func.R", "RNAseq.enrichr.Func.R")
  for (f in files) {
    suppressPackageStartupMessages(sys.source(file.path(home, "module", f), envir = globalenv()))
  }
  options(spade.home = home)
  invisible(home)
}

spade_ercc_control_file <- function() {
  f <- file.path(getOption("spade.home", "/opt/spade"), "data", "ERCC_Controls_Analysis.txt")
  if (!file.exists(f)) stop("ERCC_Controls_Analysis.txt not found: ", f)
  f
}

spade_version <- function() {
  home <- getOption("spade.home", "/opt/spade")
  vf <- file.path(home, "VERSION")
  if (file.exists(vf)) return(trimws(readLines(vf, n = 1, warn = FALSE)))
  d <- file.path(home, "DESCRIPTION")
  if (file.exists(d)) return(as.character(read.dcf(d, fields = "Version")[1, 1]))
  "unknown"
}

# Parse "TRUE"/"T"/"yes"/"1" style flags from the command line.
spade_flag <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  toupper(as.character(x)) %in% c("T", "TRUE", "YES", "Y", "1")
}

spade_write_tsv <- function(x, file) {
  write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  cat("\n", basename(file), "[saved]\n")
  invisible(file)
}

######################################################
# showMsgAndLog (from RNAseqERCC v1.1 stdRNAseq.Func.R)
showMsgAndLog <- function(logFile, msg = msg,
                nStart=0, tStart=0, EOL=0,
                verbose = TRUE, timeStamp = "left",append=TRUE)
{
  # enhanced version, 3/11/2020 4:59 PM
  if (tStart>0){
    msg<- paste(paste(rep("\t",tStart),collapse=""),msg, sep="")
  }
  if (grepl("left",timeStamp,ignore.case=T)){
   msg<- paste0("[",format(Sys.time(), '%H:%M:%S'),"] ",msg)
  }else  if (grepl("right",timeStamp,ignore.case=T)){
   msg<- paste0(msg, " [",format(Sys.time(), '%H:%M:%S'),"]")
  }

  if (EOL>0){ # end of line
     msg<- paste(msg, paste(rep("\n",EOL),collapse=""), sep="")
 }
  if (nStart>0){
    msg<- paste(paste(rep("\n",nStart),collapse=""),msg, sep="")
  }
  if (verbose) {
    cat(msg)
  }
  if (!file.exists(logFile)){
    append <- FALSE
  }
  cat(file = logFile, text = msg, append = append)
}

######################################################
# ReadYaml (from RNAseqERCC v1.1 stdRNAseq.Func.R)
ReadYaml <- function(yamlFile){
    suppressMessages(library("stringr"))
    if (!file.exists(yamlFile)){
      cat("\n Error:yaml file not found:", yamlFile)
      stop("Exit...")
    }else{
        yamlFile <- normalizePath(yamlFile)
    }
    con <- file(yamlFile, "r")
    tmp<- readLines(con = yamlFile, warn=FALSE)
    close(con)
    df<- str_split_fixed(tmp,":", 2)
    df[,1]<- trimws(df[,1])
    df[,2]<- trimws(df[,2])
    yaml <- setNames(as.list(df[,2]), gsub("[[:space:]]","_", df[,1]))
    if ("SRM" %in% names(yaml)){
        yaml$SRM <- gsub("[[:space:]]","_",yaml$SRM)
    }
    return(yaml)
}

######################################################
# SaveYaml (from RNAseqERCC v1.1 stdRNAseq.Func.R)
SaveYaml <- function(yaml, outFile=NULL){
  if(!is.list(yaml)){
    cat("\n Error: Invalid yaml list!\n")
    return (FALSE)
  }
  if(!is.list(yaml)){
    cat("\n Error: undefined outFile!\n")
    return (FALSE)
  }
  res<- paste0(paste(paste(names(yaml),": ",yaml, sep=""),collapse="\n"),"\n")
  cat(file=outFile,text=res,append=F)
  return (TRUE)
}

######################################################
# GetComparisons (from RNAseqERCC v1.1 stdRNAseq.Func.R)
GetComparisons <- function(meta, outFile=NA){
  if (! "GROUP" %in% colnames(meta) | class(meta) != "data.frame"){
       cat("\n Error: invalid meta data file!\n")
       stop("Exiting...")
  }
  groups<- as.character(gsub(" ","", unique(meta$GROUP)))
  if (any(grepl("_vs_",groups))){
       cat("\n Error: invalid meta data file! \"_vs_\" string is not allowed in meta$GROUP!\n")
       stop("Exiting...")
  }
  if(length(groups)==1) {	# invalid metadata file
		cat("\n Error: Only one GROUP provided in your metadata file! \n")
        stop("Exiting...")
  }
  myPairs <- combn(groups,2)
  comps<-paste(myPairs[1,],"_vs_",myPairs[2,],sep="")
  if (!is.na(outFile)){
    writeLines(comps,outFile)
  }
  return(comps)
}

######################################################
# GetContrastMatrix (from RNAseqERCC v1.1 stdRNAseq.Func.R)
GetContrastMatrix <- function(comparisons,design){
  # make limma contrasts
  comparisons <- gsub(" ","",unlist(comparisons))
  if(!all(grepl(".*=\\(.*-.*\\)",comparisons)) & !all(grepl(".*_[vV][sS]_.*",comparisons))){
     cat("\n Error: invalid comparisons!")
     stop("Terminated run.")
  }else if (all(grepl(".*=\\(.*\\)",comparisons))){
    x<- paste(comparisons, collapse=",")
    prestr="makeContrasts("
    poststr=",levels=design)"
    commandstr=paste0(prestr,x,poststr)
    res<- eval(parse(text=commandstr))
  }else{
    mat <- NULL
    for(n in 1:length(comparisons)) {
        keywords <- unlist(strsplit(comparisons[n], "_[Vv][Ss]_"))
        if (keywords[1] == keywords[2]){
            invalid_comp_inds <- c(invalid_comp_inds, n)
            cat(paste0("** Invalid comparison! ", comparisons[n], " [skipped]**\n" ))
            next
        }
        cfNames <- paste0(keywords[1],"vs",keywords[2])
        mat <- c(mat, paste0( cfNames,'=(',keywords[1],'-',keywords[2],')'))
    }
    x<- paste(mat, collapse=",")
    prestr="makeContrasts("
    poststr=",levels=design)"
    commandstr=paste0(prestr,x,poststr)
    res<- eval(parse(text=commandstr))
  }
   res
}

######################################################
# CalAvgByGrp (from RNAseqERCC v1.1 stdRNAseq.Func.R)
CalAvgByGrp <- function(cpm, meta, id='ID'){
    # calculate average CPM by group label
	avg<- NULL; 	cn<- NULL
    colnames(meta)<- toupper(colnames(meta))
	for(group in unique(meta$GROUP)) {
		cols_grp<- colnames(cpm)[colnames(cpm) %in% meta[meta$GROUP==group,id]]
		if (length(cols_grp) == 1){
			grp.mean <- cpm[,cols_grp]
		}else{
			grp.mean <- rowMeans(cpm[,cols_grp])
		}
		tmpDF<-data.frame(grp.mean)
		if (is.null(avg)){
			avg<- tmpDF
			cn<-  group
		}else{
			cn<- c(cn,group)
			avg<- cbind(avg,tmpDF)
		}
	}
	colnames(avg)<-paste(cn,".AveExpr",sep='')
    avg<- round(avg,3)
	return(avg)
}

######################################################
# ReadReport (from RNAseqERCC v1.1 stdRNAseq.Func.R)
ReadReport <- function(reportFile){
    con <- file(reportFile,"r") #to disable waring "incomplete final line"
    items <- readLines(con,warn=FALSE)
    close(con)
    headerLine <- items[3] # read h
    headerLine<- unlist(strsplit(gsub(" ","",headerLine),"[|]"))
    if (headerLine[1] == ""){
        headerLine[1] <- "NA"
    }
    headerLine <- gsub('[(%)]',"",headerLine)
    # to avoid erorr: incomplete final line found by readTableHeader
	dat <- read.table(text=items,sep="|",fill=TRUE,skip=4,check.names=F)
	colnames(dat) <- gsub(" ","",headerLine)
    dat$SAMPLE <- gsub(" ","",dat[,"SAMPLE"] )
    dat <- dat[!duplicated(dat$SAMPLE), ]
    dat<-dat[, !is.na(headerLine) & headerLine !="NA"]
    keptColNames <- c("SAMPLE","STATUS", "READS_RAW","READS_TRIMMED","READS_MAPPED","MAPPED","UNMAPPED",
                        "DUPLICATION","BASE_EXON", "BASE_INTRON", "BASE_INTERGENIC","BASE_CODING","BASE_UTR",
                        "BASE_RIBOSOMAL","FRAGMENT_COUNT")
    usedColNames <- intersect(keptColNames,colnames(dat))
    if (! "FRAGMENT_COUNT" %in% usedColNames){ # 3/28/2023
      dat$FRAGMENT_COUNT <- NA
    }
    dat$STATUS <- gsub(" ","",dat$STATUS)   ### to be check
    dat <- dat[, sort(colnames(dat))]  ### order columns by name
    return (dat)
}

######################################################
# MergeReports (from RNAseqERCC v1.1 stdRNAseq.Func.R)
MergeReports <- function(reportFiles, outFile=NULL){
    if (! is.vector(reportFiles) | length(reportFiles)<2 ){
      cat("\nError: Invalid report filelist.\n")
      return (FALSE)  # exit with error
    }
    if ( is.null(outFile)){
      cat("\n[MergeReports] Error: Undefined outFile.\n")
      return (FALSE)  # exit with error
    }

    locs<- NULL; dat <- NULL
    for(reportFile in reportFiles) {
        con <- file(reportFile,"r")
        items <- readLines(con,warn=FALSE) # to disable waring "incomplete final line"
        locs <- c(locs,items[1])
        close(con)
        tmp<- ReadReport(reportFile)
        if (is.null(dat)){
           dat <- tmp
        }else{

           if (!all(colnames(dat)==colnames(tmp))){
             cat("\n[MergeReports] Warn: report files have different format.\n")
             #return (FALSE)  # exit with error # 
           }
           dat<-bind_rows(dat,tmp) # [2024,July 24: merged df with column mismath which happens with different automapper report format]
         }
      }
     if (! is.null(outFile)){
        #-----------------------------mapping QC table--------------------------
        write.table(dat,outFile, sep="|",quote=F,row.names=F, col.names=T,append = F)
        con <- file(outFile,"r")
        items <- readLines(con,warn=FALSE) #
        close(con)
        #------------------------rewrite report like original format ------------
        cat(file=outFile,text=paste0(paste(locs,collapse=","),"\n\n"),append=FALSE) # lines#1-2 destination line followed by a blank line
        cat(file=outFile,text=paste0(items[1],"\n"),append=TRUE)                     # line#3  header
        cat(file=outFile,text=paste0(gsub("[A-Za-z_]","-",items[1]),"\n"),append=TRUE) # line#4: insert line |---|---
        cat(file=outFile,text=paste0(paste(items[2:length(items)],collapse="\n"),"\n"),append=TRUE) # real content of QC report
     }
     return (TRUE)
}

######################################################
# ParseMappingReport (from RNAseqERCC v1.1 stdRNAseq.Func.R)
ParseMappingReport <- function(reportFile, metaFile, outFile=NA)
{
    used<- ReadReport(reportFile)
    rownames(used) <- used$SAMPLE
    if ("FRAGMENT_COUNT" %in% colnames(used)){
        if (sum(is.na(used$FRAGMENT_COUNT)) == nrow(used)){
          used$FRAGMENT_COUNT <- NULL
        }
    }
    used$BASE_EXON <- as.numeric(used[,grep("CODING",colnames(used))]) + as.numeric(used[,grep("UTR",colnames(used))])
    meta<- ReadMeta(metaFile)$used  # SPADE ReadMeta() returns a list
    meta<- meta[,c("ID","NEWNAME","GROUP")]
    commonSamples <- intersect(used$SAMPLE,meta$ID)
    used<-used[commonSamples,]  #take subset according to meta$ID
    if(! all(meta$ID %in% used$SAMPLE)) {
        cat("\nERROR: meta$ID doesn\'t match report$SAMPLE !")
        cat("\n\n",setdiff(meta$ID,rownames(used)))
        stop("meta$ID does not match the mapping report SAMPLE column")
    }
    merged <- merge(x=meta,y=used,by.x="ID",by.y="SAMPLE")
    if (all(merged$ID == merged$NEWNAME)){
      merged$NEWNAME <- NULL
    }
    colnames(merged)[1] <- "SAMPLE"
    if(!is.na(outFile)){
      write.table(merged,outFile,quote=F,col.names=T,row.names=F,sep="\t")
    }
    return(merged)
}

######################################################
# MergeRSEMs (from RNAseqERCC v1.1 stdRNAseq.Func.R)
MergeRSEMs <- function(RSEMfiles, outFile=NULL, suffix="_RSEM_.*txt"){
    if (! is.vector(RSEMfiles) | length(RSEMfiles)<2 ){
      cat("\nError: Invalid RSEM RSEMfiles.\n")
      return (FALSE)  # exit with error
    }
    if ( is.null(outFile)){
      cat("\nError: Undefined outFile.\n")
      return (FALSE)  # exit with error
    }
    objs<- NULL
    keptColNames<-  c("geneID", "geneSymbol","bioType","annotationLevel")
    numsMatched <- NULL
    numRows <- NULL
    n<- 0
    for(inFile in RSEMfiles) {
        n <- n+1
        prefix<- paste0(gsub(suffix,"", basename(inFile)),"_",n)
        cat("\n",basename(inFile))
        dat <- read.table(inFile, sep="\t",header=TRUE,fill=TRUE,stringsAsFactors = FALSE, quote="",row.names=NULL ,check.names=F)
        assign (x=prefix, value=dat)
        objs <- c(objs, prefix)
        numsMatched <- c(numsMatched,  sum(keptColNames %in% colnames(dat)))
        numRows <- c(numRows, nrow(dat))
   }
    if (length(unique(numsMatched)) >1){
       cat("\nError: found different format of colnames in your RSEM files.\n")
       cat("\n\t Are you merging SRMs from different species/reference genomes?\n")
       return (FALSE)
    }
    if (length(unique(numRows)) >1){
       cat("\nError: found different number of rows in your RSEM files.\n")
       cat("\n\t Are you merging SRMs from different species/reference genomes?\n")
       return (FALSE)
    }
    matchedColNames <- keptColNames[keptColNames %in% colnames(dat)]
    if (!all(c("geneID", "geneSymbol") %in% matchedColNames)  ){
       cat("\nError: geneID, geneSymbol must be present as column names in your RSEM files.\n")
       return (FALSE)
    }
    finalDat <- NULL
     for(x in 2:length(objs)) {
        nextDat <- get(objs[x])
        if (! exists("curDat")){
          curDat <- get(objs[1])
        }
       if (exists("curDat")){
          nextDat<- nextDat[, ! colnames(nextDat) %in% c("geneSymbol","bioType","annotationLevel")]
          finalDat <- merge(x=curDat, y=nextDat, by.x="geneID", by.y="geneID")
          curDat <- finalDat
        }
     }
     write.table(finalDat,outFile, sep="\t",quote=F,row.names=F, col.names=T)
     return(file.exists(outFile))
}

######################################################
# PlotHeatmap (from RNAseqERCC v1.1 stdRNAseq.Func.R)
PlotHeatmap <- function(diff_table, meta,genes=NULL , outFile=NA, title=NA)
{
    if (length(genes)==0){
        return (NULL)
    }
     if (class(diff_table) != "data.frame"){
        cat("\n** Error: parameter [diff_table] must be a data.frame! **\n\n")
        stop("Exit...")
      }
      if (class(meta) != "data.frame"){
        cat("\n** Error: parameter [meta] must be a data.frame! **\n\n")
        stop("Exit...")
      }
      #validate input
      colnames(meta) <- toupper(colnames(meta))
      if (!all(c("ID","GROUP","COLOR") %in% colnames(meta))) {
        cat("\nInvalid meta file. ID,GROUP or COLOR column(s) not found.\n")
        stop("exit.")
      }
      isNewName <- length(intersect(colnames(diff_table), meta$NEWNAME )) > length(intersect(colnames(diff_table), meta$ID ))
      if (isNewName){
        cat("\n Use NEWNAME in diff_table table..\n")
        keptSamples <-  intersect(meta$NEWNAME,colnames(diff_table))
        rownames(meta) <- meta$NEWNAME
      }else{
        rownames(meta) <- meta$ID
        keptSamples <-  intersect(meta$ID,colnames(diff_table))
      }
      if (length(keptSamples) ==0){
        cat("\n ** Error: ID in metadata table doesn\'t match colnames in input matrix! **\n")
        stop ("Invalid metadata file!\n")
      }
      meta <- meta[keptSamples, ]
      dat <- diff_table[, keptSamples]
      colnames(dat) <- meta$NEWNAME #

      rm(list=c("diff_table"))
      #colnames(dat) <- keptSamples
     #-----------------------------------------------------
     usedGenes <- intersect(genes, rownames(dat))
     if (length(usedGenes)==0 ){
        cat("\n ** Error: no topn DEGs**\n")
         return(NA)
     }
     if (length(usedGenes)<2 ){
        cat("\n ** Warning: not enough topn DEGs [n=",length(usedGenes),"]**\n")
        cluster_rows <- FALSE
     }else{
        cluster_rows <- TRUE
    }
    used <- dat[usedGenes, ]
    if(length(unique(meta$COLOR)) == length(unique(meta$GROUP))){
        # automatic color. each group has different color
        annotation <- data.frame(GROUP = factor(meta$GROUP))
        rownames(annotation) <- meta$NEWNAME #meta$ID
        # change the color of annotation
        GROUP <- unique(meta$COLOR)
        names(GROUP) <- unique(meta$GROUP)
        anno_colors <- list(GROUP = GROUP)
    }else{
        # custom color. multiple groups could share same color
        annotation <- data.frame(GROUP = factor(paste(meta$COLOR,meta$GROUP,sep="_")))  #unique pseudo group
        uniqCombs <- meta[,c("COLOR","GROUP")]
        uniqCombs$pseudoGroup <- paste(meta$COLOR,meta$GROUP,sep="_")
        uniqCombs <- uniqCombs[!duplicated(uniqCombs$pseudoGroup), ]
        annotation <- data.frame(GROUP = factor(paste(meta$COLOR,meta$GROUP,sep="_")))  #unique pseudo group
        rownames(annotation) <- meta$NEWNAME #meta$ID
        GROUP <- uniqCombs$COLOR
        names(GROUP) <- uniqCombs$pseudoGroup
        anno_colors <- list(GROUP = GROUP)
    }
    hmCols <- c(rgb(seq(0.01, 0.99, 0.01), seq(0.01, 0.99, 0.01), 1),
                "white",
                rgb(1, 1 - seq(0.01, 0.99, 0.01), 1 - seq(0.01, 0.99, 0.01)))
    plotWidth <- 3+log10(ncol(used))*5
    plotHeight <- 3+ log10(nrow(used))*4
    fontsize <- (19-log2(nrow(used)))/2
    #cluster_rows <- dist(used, method = "euclidean")
    bk <- NA
    hm <- pheatmap(used,main=title,
      border_color = "grey", scale = "none", cluster_rows = cluster_rows, cluster_cols = T, fontsize=fontsize,
      cellheight = NA, cellwidth = NA,breaks = bk, show_rownames = T,
      display_numbers = F, color = hmCols, filename = outFile, clustering_distance_rows = "euclidean",
      annotation_col  = annotation, annotation_colors = anno_colors, clustering_method = "ward.D",
      width = plotWidth, height = plotHeight  )
    return(hm[[4]])
 }
