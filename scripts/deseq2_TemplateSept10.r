#!/usr/bin/Rscript

############ add in Ivan's DEG formatting package
#'** note the intended polarity of the comparison i.e Q175-Q20 would be positive if comparison is higher'

# On my mac I had to: install bimatRt, and openxlsx manually
# source('https://bioconductor.org/biocLite.R') biocLite('biomaRt')
# library('biomaRt') install.packages('openxlsx', dependencies=TRUE)

# open the tarball manually then install it install.packages(pkgs =
# '~/Documents/deg_reports/formatDESeq2', repos=NULL, type='source',
# dependencies=TRUE)

#'Ensembl', 'Gene', 'Description', 'Mean', 'Log2FC', 'Log2FC_SE','Stat','PValue','AdjPValue',\t'AbsLog2FC'
# “Ensembl”, “Gene”, “Log2FC”, “PValue”, and “AdjP”.
library(formatDESeq2)

merge_deseq2_results(file_dir="~/Documents/JOB_Current/code/DEG_format/", 
    deseq2_colnames=c(geneid="EnsemblID", gene="GeneName", lfc="log2FC", 
        pval="p", adjp="fdr"))

################################# 


# Load libraries
library(DESeq2)
library(biobroom)
library(patchwork)
library(plyr)
library(dplyr)
library(toolboxR)
library(tidyverse)
library(pca3d)

rm(list=ls())
options(stringsAsFactors=FALSE)
platform_name <- "RNASeq"
basedir <- "~/Documents/JOB_Current"
setwd(paste0(basedir, "/code/"))
outdir <- paste0(basedir, "/data/")
source("analysis/dan_deseq_functions.R")
source("analysis/util.r")  # using Jike's PCA utility


######## read TRAP DATA add new short name function to make the short name
ColShortName <- function(df) {
    nameLs <- names(df)  # get all of the colum names
    col_idx <- 0
    for (col in nameLs) {
        col_idx <- col_idx + 1
        colstr_pts <- strsplit(col, "[.]")[[1]]
        last_p <- length(colstr_pts)
        nameLs[col_idx] <- paste(colstr_pts[1], colstr_pts[last_p - 2], 
            colstr_pts[last_p - 1], colstr_pts[last_p], sep="_")
    }
    # overwrite the old names with the updated list
    names(df) <- nameLs
    return(df)
}

# add short name to metadata and additional slicers.
meta2 <- auto_read("../data/SampleManifest.xlsx") %>% 
  mutate(NewName=paste(prefix, Age, Mouse.Line, Replicate, sep="_")) %>% 
  mutate(Exp_type=factor(paste(prefix, Age, Mouse.Line, sep="_"))) %>% 
  mutate(Technique=factor(prefix))

data1_counts <- auto_read("../data/data1_counts.xlsx") %>% 
  column_to_rownames("geneid") %>% 
  ColShortName()

data_2counts <- auto_read("../data/data2_counts.xlsx") %>% 
  column_to_rownames("geneid") %>% 
  ColShortName()
# HSVcounts<-as.data.frame(HSVcounts)

Allcounts <- merge(data2_counts, data1_counts, by=0)
Allcounts <- Allcounts %>% 
  column_to_rownames("Row.names")

data1_normCounts <- auto_read("../data/data1_FPKM.xlsx") %>% 
  column_to_rownames("geneid") %>% 
  ColShortName()

data2_normCounts <- auto_read("../data/data2_FPKM.xlsx") %>% 
  column_to_rownames("geneid") %>% 
  ColShortName()

# excludedSamps=c('X1','X2')

#### get the gene symbol translation table
geneID2symbol <- read.delim("../data/ensmouse_table.txt")

# clean up duplicated indexes
dupInd <- which(duplicated(geneID2symbol[, 1]))
if (length(dupInd) > 0) {
    geneID2symbol <- geneID2symbol[-dupInd, ]
}
rownames(geneID2symbol) <- geneID2symbol[, 1]

####### compile a list of sample IDs and concatenate with control/ treatment
####### asssignment
analysisWrapper <- function(ctrllvl, trtlvl, File_title, trtCol, sampleCol, 
    sample_data, rawCounts2, gene_names_ref, base_name) {
    
    ### clean samples and set up 'Labels' list in order of columns
    excludedSamps <- c("data1_20_1", "data1_20_2", "data1_111_2", "data2_50_1", 
        "data2_50_4")
    
    levels_used <- c(ctrllvl, trtlvl)
    
    # Get lables list by stripping off the replicate number
    labels <- gsub("(.*)_.*", "\\1", c(colnames(rawCounts2)))
    
    # Cast this as a dataframe insted of tibble to suport the required row
    # names Change the names to the simplified names and output as counts
    # table > includes all the data
    tempDF <- as.data.frame(rawCounts2)
    write.table(as.data.frame(tempDF), file=paste("Reps_Allcounts.txt"), 
        sep="\t", quote=F, row.names=T)
    names(tempDF) <- labels
    write.table(as.data.frame(tempDF), file=paste("Norep_Allcounts.txt"), 
        sep="\t", quote=F, row.names=T)
    
    ## filter the levels and labels by the subset actually used
    samples_bool <- !(c(colnames(rawCounts2)) %in% excludedSamps)
    levels_bool <- (labels %in% levels_used)
    combined_bool <- as.logical(samples_bool * levels_bool)
    
    # filter out excluded samples
    labels_filt <- labels[combined_bool]
    rawCounts2_filt_cols <- colnames(rawCounts2)[combined_bool]
    
    # select the columns being used in the current expereiment cast this as
    # a dataframe insted of tibble to suport the required row names
    rawCounts2 <- as.data.frame(rawCounts2)
    rawCounts2 <- rawCounts2[, rawCounts2_filt_cols]
    
    ###### Run It ###
    dsRe <- deseq2v2(rawCounts2, labels_filt, ctrl.label=ctrllvl, case.label=trtlvl, 
        FCs=1.5, p.fdr="fdr", ofBaseName=base_name)
    
    dsRe <- dsRe %>% 
      mutate(EnsemblID=GeneID) %>% 
      join(geneID2symbol, by="EnsemblID", match="all") %>% 
      select("EnsemblID", "GeneName", "Description", everything(), -"GeneID")
    
    ## Ensembl”, “Gene”, “Log2FC”, “PValue”, and “AdjP”.  get the gene
    ## symbol translation table
    ## geneID2symbol=read.delim('ensebl_geneid_to_symbol.tsv')
    ## dupInd=which(duplicated(geneID2symbol[,1])) if(length(dupInd)>0){
    ## geneID2symbol=geneID2symbol[-dupInd,]}
    
    # rownames(geneID2symbol)=geneID2symbol[,1]
    
    # for the heatmap plots, read the normalized counts
    # normCounts=readRDS(paste0(BaseName,'.normalized.counts.rds'))
    
    ############################ analyze the results and make figs bind the gene names and gene
    ############################ descriptions to the output
    
    
    ### OUTPUT REPORTS UNDER CUTOFFS full data set
    write.table(dsRe, file=paste(File_title, "Full_deseq2.txt"), sep="\t", 
        quote=F, row.names=F)
    
    # FDR
    dsReFDR05 <- dim(dsRe[which(dsRe$fdr <= 0.05), ])
    print(paste("after filter with FDR<=.05, the count is: ", dim(dsRe[which(dsRe$fdr <= 
        0.05), ])[1]))
    write.table(dsReFDR05, file=paste(File_title, "FDR<=0.05.deseq2.txt"), 
        sep="\t", quote=F, row.names=F)
    
    # padjBH
    dsReBH <- (dsRe[which(dsRe$padj_BH <= 0.05), ])
    print(paste("after filter with padj_BH<=.05, the count is: ", dim(dsReBH)[1]))
    write.table(dsReBH, file=paste(File_title, "_adj_BH<=0.05.deseq2.txt"), 
        sep="\t", quote=F, row.names=F)
    
    dsRe1 <- dsRe[which(dsRe$p <= 0.05 & abs(dsRe$FC) >= 1.5), ]
    print(paste("after filter with p05 and FC, the count is: ", dim(dsRe1)[1]))
    replace(dsRe1$geneSymbol, is.na(dsRe1$geneSymbol), "")  #need to remove NAs to make this work
    # Use gene symbol unless it is missing, then fill in with the geneID
    # (Ensembl id) Not sure that you need this >->->->-> >>>#
    # dsRe1[which(dsRe1$geneSymbol==''),'geneSymbol']=dsRe1[which(dsRe1$geneSymbol==''),'GeneID']
    write.table(dsRe1, file=paste(File_title, "_p05fc1.5_deseq2.txt"), 
        sep="\t", quote=F, row.names=F)
    
    ############################## Heatmap ############ -takes in norm counts from before
    normCounts <- readRDS(paste0(base_name, case.label, "vs", ctrl.label, 
        ".normalized.counts.rds"))
    
    normCounts <- normCounts[dsRe$GeneID[1:100], ]
    # > rownames(normCounts)=dsRe$geneSymbol[1:100] ####NEED TO ADD GENE
    # SYMBOL + DESCR
    
    loadPackage(c("pheatmap", "gplots", "RColorBrewer"))
    # library(pheatmap)
    
    # Match the filtered column labels to the experiment type.
    # colnames(normCounts) <-labels_filt # strip the replicate information
    # from the columns
    
    # mat_col <- data.frame(group =
    # sample_data[colnames(normCounts),trtCol])
    mat_col <- data.frame(group=colnames(normCounts), labels_filt)
    
    # Data frame with column annotations.
    mat_col <- data.frame(group=labels_filt)
    rownames(mat_col) <- colnames(normCounts)
    
    # List with colors for each annotation.
    nUniq <- length(unique(mat_col$group))
    mat_colors <- list(group=brewer.pal(nUniq, "Set1")[1:nUniq])
    names(mat_colors$group) <- unique(mat_col$group)
    
    png(file=paste(File_title, ".heatmap.png"), width=1000, height=800)
    pheatmap(normCounts, col=greenred(75), scale="row", annotation_col=mat_col, 
        annotation_colors=mat_colors, border_color=NA, fontsize=7, 
        main=paste(File_title, "RNA-seq-", trtlvl, "_vs.", ctrllvl, "_top_100_DE_genes"))
    dev.off()
    return(dsRe)
}

# requres the initial file to be there
complieTargetGene <- function(data, outdir, filename, targetGene, currExpt) {
    Mtx <- auto_read(paste0(outdir, filename))
    Ln <- data[targetGene, ]
    # add name of the experiment
    Ln["ExperimentName"] <- currExpt
    print(Ln)
    # add it to the auto matirix
    HttMtx <- rbind(Mtx, Ln)
    write.table(Mtx, file=paste0(outdir, filename), sep="\t", quote=F, row.names=F)
}

################################################# 
targetGene1 <- "ENSMUSG000000XXXX"
targetGeneNm1 <- "HttMatrix.txt"
rawCounts2 <- Allcounts
trtCol <- "Exp_type"
sampleCol <- "NewName"
sample_data <- meta2
gene_names_ref <- geneID2symbol
base_name <- "Base_tech"

# ctrllvl='X20' trtlvl='X111' trtlvl_50='X50'

File_title <- "Comparison1-X50vsX20"
deseq2_out <- analysisWrapper("X20", "X50", File_title, trtCol, sampleCol, 
    sample_data, rawCounts2, gene_names_ref, base_name)
complieTargetGene(deseq2_out, outdir, targetGeneNm1, targetGene1, File_title)

###################### EXECUTE THE RUNS ####################### rawCounts2=HSVcounts

File_title <- "Comparison2-X50vsX20"
Basename <- "BaseExpt"
analysisWrapper("X20", "X50", File_title, trtCol, sampleCol, sample_data, 
    rawCounts2, gene_names_ref, base_name)
complieTargetGene(deseq2_out, outdir, targetGeneNm1, targetGene1, File_title)

########### post processing you will have to arrange files for the DEG report
########### step DEG_format > comp1.txt comp2.txt comp3.txt ...
########### Reps_Allcounts.txt

library(formatDESeq2)
merge_deseq2_results(file_dir="~/Documents/JOB_Current/code/DEG_format/", 
    deseq2_colnames=c(geneid="EnsemblID", gene="GeneName", lfc="log2FC", 
        pval="p", adjp="fdr"))

########## comparing last time and this time note there are 2M and 6m levels for
########## some of these

function() 
AssayPairs <- list(c("X50_sig_results.txt", "Y50_6M_results.txt"), c("X111_sig_results.txt", 
    "Y111_2M_results.txt"), c("X111_sig_results.txt", "Y111_6M_results.txt"))

datafolder <- "~/Documents/JOB_Current/code/Result_comparisons/"
i <- 1
CurrPr <- AssayPairs[[i]]
first_assay <- CurrPr[2]
second_assay <- CurrPr[1]

second_data <- auto_read(paste0(datafolder, second_assay)) %>% 
filter()

# dplyr::select('Ensembl','Gene','AdjP')

first_data <- auto_read(paste0(datafolder, first_assay)) %>% 
dplyr::select("Ensembl", "Gene", "AdjP")

merged_data <- first_data %>% 
  full_join(second_data, by="Ensembl", suffix=c(".first", ".second")) %>%
    mutate(NewName=paste(prefix, Age, Mouse.Line, Replicate, sep="_")) %>% 
    mutate(Exp_type=factor(paste(prefix, Age, Mouse.Line, sep="_"))) %>% 
    mutate(Technique=factor(prefix))

# combine