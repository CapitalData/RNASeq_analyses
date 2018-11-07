loadPackage <- function(packageNames=NULL, update.these=F, update.all=F, 
    update.bioc=F, install.bioc=F, lib.loc=NULL, ask=F, ...) {
    
    source("http://www.bioconductor.org/biocLite.R")
    tryCatch({
        if (install.bioc) 
            biocLite(ask=ask, lib.loc=lib.loc)
        if (update.all) 
            biocLite(ask=ask, lib.loc=lib.loc)  #update current packages
        if (update.bioc) 
            biocLite("BiocUpgrade", ask=ask, lib.loc=lib.loc)  #update bioconductor
    }, error=function(cond) {
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)
    }, warning=function(cond) {
        message("Here's the original warning message:")
        message(cond)
        return(NULL)
    })
    
    for (pn in packageNames) {
        if (!require(pn, character.only=TRUE) || update.these) {
            warning("package ", pn, " is not installed, installing now!")
            biocLite(pn, ask=ask, lib.loc=lib.loc, ...)
            library(pn, character.only=TRUE)
        }
    }
}

# yx.ratio is the fig ratio, can be 'yx.ratio=2/3' hjust=0.5: center
# the plot title
myggTheme <- function(fs=15, hideLegend=F, hideLegendTitle=F, legendPosition="top", 
    hideXlabel=F, hideYlabel=F, rotateYlabel=F, rotateXtickLabel=F, 
    rotateAngle=90, yx.ratio=NULL) {
    library(ggplot2)
    if (hideLegend) 
        legendPosition <- "none"
    myT <- theme(axis.text=element_text(size=fs), axis.title=element_text(size=fs + 
        3), plot.title=element_text(size=fs + 5, hjust=0.5), legend.title=element_text(size=fs + 
        3), legend.text=element_text(size=fs), legend.position=legendPosition)
    if (hideLegendTitle) 
        myT <- myT + theme(legend.title=element_blank())
    if (hideXlabel) 
        myT <- myT + theme(axis.title.x=element_blank())
    if (hideYlabel) 
        myT <- myT + theme(axis.title.y=element_blank())
    if (rotateYlabel) 
        myT <- myT + theme(axis.title.y=element_text(angle=rotateAngle, 
            hjust=1))
    if (!is.null(yx.ratio)) 
        myT <- myT + theme(aspect.ratio=yx.ratio)
    if (rotateXtickLabel) 
        myT <- myT + theme(axis.text.x=element_text(angle=rotateAngle, 
            hjust=1))
    
    return(myT)
}

# Remove rows or columns of 0 variance gep: a numberical matrix or
# data.frame
rmVar0 <- function(gep, margin=1) {
    if (!is.matrix(gep)) {
        gep <- as.matrix(gep)
    }
    
    vars.genes <- apply(gep, margin, var, na.rm=T)
    ids.var0 <- which(vars.genes == 0)
    if (length(ids.var0) > 0) {
        if (margin == 1) {
            gep <- gep[-ids.var0, ]
        } else {
            gep <- gep[, -ids.var0]
        }
    }
    return(gep)
}

# gep: a matrix or data frame defined by gepRowIsSample samps: a data
# frame: row: samples; col: features; rownames: sampleID doScale: Scale
# variance or not revPCdir: reverse the direction of PCs
do.PCA <- function(gep, samps, gepRowIsSample=F, revPCdir=F, doScale=T, 
    PCs=1:3, ...) {
    if (!gepRowIsSample) {
        gep <- t(gep)
    }
    comSamps <- intersect(rownames(gep), rownames(samps))
    if (length(comSamps) < 3) 
        stop("Less than 3 shared samps, something wrong!")
    gep <- gep[comSamps, ]
    samps <- samps[comSamps, , drop=F]
    
    gep <- rmVar0(gep, 2)
    pcaOut <- prcomp(gep, scale=doScale, ...)
    pcScores <- predict(pcaOut)
    if (revPCdir) 
        pcScores <- -pcScores
    pcaSmry <- summary(pcaOut)
    var.perc <- round(100 * pcaSmry$importance, digits=1)
    
    if (all(rownames(pcScores) == rownames(samps))) {
        pcScores <- cbind(pcScores[, PCs], samps)
    } else {
        stop("pcaScore rownames not equal samps rownames!")
    }
    return(list(score=pcScores, var=var.perc))
}



# pcaReturn: returned list from do.PCA() figMain: the title of the plot
# ofBaseName: the base name for output files classCol: the column name
# in sample table used for point color shapeCol: the column name in
# sample table used for point shpae
draw.PCA <- function(pcaReturn, PC=1:3, figMain=NULL, ofBaseName=NULL, 
    classCol=NULL, shapeCol=NULL, pointSize=3, textSize=20) {
    
    pcScores <- pcaReturn[[1]]
    var.perc <- pcaReturn[[2]]
    
    loadPackage("ggplot2")
    pcInd <- combn(1:3, 2)
    for (i in PC) {
        pcx <- paste("PC", pcInd[1, i], sep="")
        pcy <- paste("PC", pcInd[2, i], sep="")
        
        if (is.null(figMain)) {
            figTitle <- paste(pcx, "vs", pcy, sep=" ")
        } else {
            figTitle <- figMain
        }
        p <- ggplot(data=pcScores, mapping=aes_string(pcx, pcy, color=classCol, 
            shape=shapeCol)) + geom_point(size=pointSize, alpha=8/10) + 
            ggtitle(figTitle) + myggTheme(textSize) + xlab(paste(pcx, " : ", 
            var.perc[2, pcInd[1, i]], "%", sep="")) + ylab(paste(pcy, 
            " : ", var.perc[2, pcInd[2, i]], "%", sep=""))
        
        ggsave(paste(ofBaseName, ".PC", pcInd[1, i], pcInd[2, i], ".png", 
            sep=""), p, dpi=600)
    }
    
    return(p)
}

# countsTable: dataframe or matrix of gene counts. cols: sample ids
# ONLY; rows: gene ids. samCondition: a vector of sample group labels
# corresponding to columns of countsTable ctrl.label: sample group
# label of control; default: the 1st level in factor(samCondition)
# p.fdr: use p or fdr for DEG useCooksFilter: if one or more samples
# for a row have a distance higher than Cook's distance, the p-value
# for the row is set to NA; see 'cooksCutoff' for detail useShrunkenFC:
# if T, the FC isn't the ratio of case/ctrl, but the shrunken fold
# change for the goal of reducing the FC of low expression genes, see
# http://genomebiology.com/content/pdf/s13059-014-0550-8.pdf
# ofBaseName: output file base name; if set, output files are saved.
# Five output files are two volcano plot, normalized gene count,
# differential results for all genes and DEGs only.

# Deseq2 set the FDR of some low expression genes to NA in order to
# increase the number of genes with significant FDR
deseq2 <- function(countsTable, samCondition, ctrl.label=NULL, case.label=NULL, 
    FCs=c(1.2, 1.5), p.fdr="p", ofBaseName=NULL, useCooksFilter=T, 
    useShrunkenFC=T) {
    if (is.null(ctrl.label)) {
        labels.uni <- unique(samCondition)
        ctrl.label <- labels.uni[1]
        case.label <- labels.uni[2]
    }
    
    loadPackage("DESeq2")
    countsTable <- as.matrix(countsTable)
    countsTable <- round(countsTable)
    samCondition <- as.data.frame(samCondition, stringsAsFactors=T)
    dds <- DESeqDataSetFromMatrix(countData=countsTable, colData=samCondition, 
        design=~samCondition)
    # this is necessary because the comparison will be made to be the 2nd
    # level over the 1st level; so this is to make sure the comparison is
    # set correctly.
    dds$samCondition <- factor(dds$samCondition, levels=c(ctrl.label, 
        case.label))
    
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds, betaPrior=useShrunkenFC)
    # or use a higher cutoff: p <- 2 m <- 10 dds <- nbinomWaldTest(dds,
    # cooksCutoff=qf(.95,p,m-p))
    res <- results(dds, cooksCutoff=useCooksFilter)
    # mcols(res)$description [1] 'the base mean over all rows' [2] 'log2
    # fold change (MAP): condition treated vs untreated' [3] 'standard
    # error: condition treated vs untreated' [4] 'Wald statistic: condition
    # treated vs untreated' [5] 'Wald test p-value: condition treated vs
    # untreated' [6] 'BH adjusted p-values'
    res <- as.data.frame(res)
    colnames(res)[c(2, 5, 6)] <- c("log2FC", "p", "fdr")
    isCtrl <- (colData(dds)$samCondition == ctrl.label)
    counts.norm <- counts(dds, normalized=TRUE)
    res[, ctrl.label] <- rowMeans(counts.norm[, isCtrl])
    isCase <- (colData(dds)$samCondition == case.label)
    res[, case.label] <- rowMeans(counts.norm[, isCase])
    res[, "FC"] <- 2^res[, "log2FC"]
    inds <- which(res[, "FC"] < 1)
    res[inds, "FC"] <- -1/(res[inds, "FC"])
    res$GeneID <- rownames(res)
    res <- res[, c("GeneID", "baseMean", case.label, ctrl.label, "FC", 
        "log2FC", "p", "fdr")]
    if (!is.null(ofBaseName)) {
        write.table(res, file=paste(ofBaseName, case.label, "vs", ctrl.label, 
            "deseq2.txt", sep="."), sep="\t", quote=F, row.names=F)
        saveRDS(counts.norm, paste0(ofBaseName, ".normalized.counts.rds"))
        volcano.plot(res$p, res$log2FC, figFileName=paste(ofBaseName, 
            case.label, "vs", ctrl.label, "volcano.png", sep="."))
        volcano.plot(res$p, res$FC, xlabStr="Fold Change", figFileName=paste(ofBaseName, 
            case.label, "vs", ctrl.label, "volcano.FC.png", sep="."))
    }
    
    for (FC in FCs) {
        degs <- res[which(abs(res[, "FC"]) > FC & res[, p.fdr] < 0.05), 
            ]
        if (!is.null(ofBaseName)) 
            write.table(degs, file=paste(ofBaseName, case.label, "vs", 
                ctrl.label, "DEG.FC", FC, p.fdr, "0.05.txt", sep="."), 
                sep="\t", quote=F, row.names=F)
    }
    
    return(res)
}


volcano.plot <- function(pvals, FC, p.cutoff=0.05, FC.cutoff=1.5, figFileName=NULL, 
    titleStr="volcano plot", ylabStr="-log10(pvals)", xlabStr="log2FoldChange") {
    inds1 <- which(!is.na(pvals) & !is.nan(pvals) & !is.null(pvals) & is.finite(pvals) & 
        pvals < 1 & pvals > 0)
    inds2 <- which(!is.na(FC) & !is.nan(FC) & !is.null(FC) & is.finite(FC))
    
    inds <- intersect(inds1, inds2)
    pvals <- pvals[inds]
    FC <- FC[inds]
    
    pointsDF <- data.frame(x=FC, y=-log10(pvals))
    pointsDF$sig <- (pvals < p.cutoff & abs(FC) > log2(FC.cutoff))
    
    loadPackage("ggplot2")
    p <- ggplot(pointsDF, aes(x=x, y=y, colour=as.factor(sig))) + 
        geom_point(alpha=1/3, show.legend=FALSE) + scale_color_manual(values=c("black", 
        "red")) + geom_hline(yintercept=c(-log10(c(0.05, 0.01, 0.001))), 
        linetype=2, lwd=1) + geom_vline(xintercept=c(FC.cutoff, -FC.cutoff), 
        linetype=2, lwd=1) + myggTheme(20) + ggtitle(titleStr) + xlab(xlabStr) + 
        ylab(ylabStr) + scale_y_continuous(breaks=scales::pretty_breaks(n=6)) + 
        scale_x_continuous(breaks=scales::pretty_breaks(n=6))
    
    if (!is.null(figFileName)) {
        ggsave(figFileName, p, dpi=600, width=7, height=10)
    } else {
        p
    }
}



##########################################################
# for GSEA we are sorting by a fold change * log10Pvalue combination as suggested in 
# https://www.ncbi.nlm.nih.gov/pubmed/20660011 (plaiser et al. 2010)
# seee https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/#appendix
# the sorting column is called FC_PvalCombo
#########################################################
# GSEA Functions                                             #
#########################################################

preFmt4fgsea <- function(dataset, format_tag, gene_name_map) {
  # takes in full output from nanostringDIFF or deseq2 
  # threshods by fold change, p value and prepares data for ingestion by fgsea
  # requres:
  # input format argument ("Deseq2format", "nanostringformat")
  # p value and fold change thresholds
  # expects 
  # gene name map with hgnc name 'Approved_symbol' 
  # and/or 'Ensembl_gene_ID' as well as Entrz gene number 'NCBI_gene_ID'
  
  # do the input specific clean up
  if (format_tag == "Deseq2format") { 
    print(paste("analyzing", format_tag, 'data', sep=" "))
    dataset<-dataset %>%
      #mutate(EntzNum=as.numeric(substring(GeneID , 5, length(GeneID)))) %>%
      mutate(ensembl_gene_id=GeneID) %>%
      left_join(gene_name_map, by='ensembl_gene_id') %>%
      filter(!is.na(FC)) %>% # remove genes with no data
      mutate(geneSymbol=hgnc_symbol)  ### using this becasue the number is not the entrezID and fails
      #mutate(geneSymbol=NCBI_gene_ID)  # fixed
  }
  
  else if (format_tag =="nanostringformat") {
    print(paste("analyzing", format_tag, 'data', sep=" "))
    dataset<-dataset %>%
      rename(hgnc_symbol=X) %>%
      rename(p=pvalue) %>%
      left_join(gene_name_map, by='hgnc_symbol') %>% # add gene names to dataset
      rowwise() %>%
      mutate(FC=BacktoFC(logFC,2)) %>%
      filter(!is.na(FC)) %>% # remove genes with no data
      mutate(geneSymbol=hgnc_symbol)  ### using this because the number is not the entrezID and fails
      #mutate(geneSymbol=NCBI_gene_ID) 
  }
}

prep4fgsea <- function(dataset){    
  # general reformatting for FGSEA ingestion
  # requires FC, p, GeneSymbol (either ENTREZ integer part or gene name)
  dataset<-dataset %>%
    mutate(rank_metric=log10(as.numeric(p))*(FC)) %>%
    select(geneSymbol, FC, p, rank_metric) %>%
    filter(!is.na(p)) %>%
    filter(geneSymbol != "") %>%
    filter(!duplicated(geneSymbol))%>% # remove any duplicated gene suymbols
    arrange(rank_metric)
}

RunGSEA <- function(filename, genes, gen_sym_type, pThresh, adjpThresh, MinSz, nperms){
  # make a named gene list of the ranking metric
  rank_list <- genes$rank_metric
  names(rank_list) <- genes$geneSymbol
  
  if (gen_sym_type == "entrez_gene") {
    # get the msigdb lists after filtering for applicable sets
    msigdb_pathways <- msigdbr(species="Homo sapiens") %>%
      filter(entrez_gene %in% genes$geneSymbol) %>%
      split(x=.$entrez_gene, f=.$gs_name)
  }
  
  if (gen_sym_type == "gene_symbol") {
    msigdb_pathways <- msigdbr(species="Homo sapiens") %>%
      filter(gene_symbol %in% genes$geneSymbol) %>%
      split(x=.$gene_symbol, f=.$gs_name)
  }
  
  fgseaRes <-  fgsea(pathways=msigdb_pathways,
                     stats=rank_list,
                     minSize=MinSz,
                     maxSize=5000,
                     nperm=nperms) %>%
    arrange(pval)
  
  #report signicant count
  print(paste("the number of postive adjusted pvalue pathways under threshold:", adjpThresh, sep=" "))
  print(sum(fgseaRes$padj < adjpThresh))
  print(paste("the number of postive pvalue pathways under threshold:", pThresh, sep=" "))
  print(sum(fgseaRes$pval < pThresh))
  
  y<-fgseaRes[fgseaRes$ES > 0,]
  topPathwaysUp <-head(y[order(y$pval),]$pathway)
  
  z<-fgseaRes[fgseaRes$ES < 0,]
  topPathwaysDown <-head(z[order(z$pval),]$pathway)
  
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(msigdb_pathways[topPathways], rank_list, fgseaRes, gseaParam=0.5)
  fwrite(fgseaRes, file=paste(filename, "fgseaRes.txt", sep=""), sep="\t", sep2=c("", " ", ""))
}

BacktoFC <- function(logval, baseval) {
  if (abs(logval) == 0){ 
    FCout=(baseval**abs(logval)) # for no FC the value is 1 (the same)
  } else {
    FCout=sign(logval)*(baseval**abs(logval))
  }
}

dropNA_avg <- function (val1, val2) {
  #take the average but if there are NAs drop whichver is an NA 
  #and divide by the number of terms remaining
  #works for two values currently
  a=c(val1,val2) 
  a=a[!is.na(a)] 
  return(sum(a)/length(a))
}

Merge_Avg <- function(L_data, R_data, merge_col) {
  L_data<-L_data %>%
    # Merged datasets - left (outer) join RNAseq to Nanostring data
    left_join(R_data, by=merge_col) %>%
    # where there are FC and p values for both, take the average and make new column (._combo)
    rowwise() %>%
    mutate(FC_combo=dropNA_avg(FC.x, FC.y)) %>%
    mutate(p_combo=dropNA_avg(p.x, p.y)) %>%
    select(-p.x, -p.y, -FC.x, -FC.y) %>%
    rename(p=p_combo) %>% 
    rename(FC=FC_combo) %>% 
    rename(geneSymbol=hgnc_symbol)
  # where there is only one value take the sigle value as the combo value
  # pass the ._combo values into the ranking function for GSEA
}


#########################################################
# Gene list functions                                        #
#########################################################

filt_list <- function(dataset, pval_TH, FC_TH) {
  dataset1<-dataset %>%
    filter(abs(FC) >= FC_TH) %>%
    filter(p<pval_TH)
  return(dataset1)
}

TwolvVenn <- function(Figname, xlist) {
  # xlist has the form   x=list(X=1:150, Y=121:180)
  venn.diagram(xlist,
               filename=paste(Figname,".tiff", sep=""),
               lwd=4,
               fill=c("cornflowerblue", "darkorchid1"),
               alpha=0.75,
               label.col="black",
               cex=3,
               fontfamily="serif",
               fontface="bold",
               cat.col=c("cornflowerblue", "darkorchid1"),
               cat.cex=3,
               cat.fontfamily="serif",
               cat.fontface="bold",
               cat.dist=c(0.03, 0.03),
               cat.pos=c(-20, 14)
  )
}

######################## ######################## ########################
################# deseq functions ######################## ###############
######################## ######################## ########################