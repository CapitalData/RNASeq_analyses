# Load libraries
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library(DESeq2)
library(biobroom)

# custom tools from Dan
# devtools::install_github("dkrozelle/toolboxR")
library(toolboxR)
library(tidyverse)

# Current requirements
# counts table with first column of identifiers named `geneid`
#

pretty_fraction <- function(v){
  total <- sum(v)
  paste0(v," (",sprintf("%.1f", v/total*100),"%)")
}

check_dir <- function(dir){
  if(!dir.exists(dir)){dir.create(dir)}
  return(dir)
}

make_rownames <- function(tbl, column_name = "geneid"){
  tbl %>%
    as.data.frame() %>%
    column_to_rownames(column_name)
}

tidy_rownames <- function(tbl, column_name = "geneid"){
  tbl %>%
    as.data.frame() %>%
    rownames_to_column(column_name)
}

# Takes a long form table and makes a data matrix with defined rownames.
# This will only work if nrow = unique(row_id)*unique(col_id)
bio_matrix <- function(df){
  row_id  <- names(df)[1]
  col_id  <- names(df)[2]
  val     <- names(df)[3]

  if(nrow(df) == length(unique(df[[col_id]])) * length(unique(df[[row_id]]))){
    df %>%
      arrange_(col_id) %>%
      spread_(col_id, val) %>%
      make_rownames(column_name = row_id)
  }else{
    message("spread operation is not possible with non-unique row identifiers")
  }


}



# result is a data frame with the columns:
#
# gene        gene ID
# baseMean    mean abundance level
# estimate    estimated log2 fold change
# stderror    standard error in log2 fold change estimate
# statistic   test statistic
# p.value     p-value
# p.adjusted  adjusted p-value
tidy_results <- function(control, test, dset){
  results(dset, contrast=c("group", test, control),
          alpha = p.adj.cutoff,
          pAdjustMethod = p.adj.method) %>%
    biobroom::tidy.DESeqResults() %>%
    mutate(sig = case_when(
      (p.adjusted < p.adj.cutoff) & (estimate >= fc_cutoff) ~ "up",
      (p.adjusted < p.adj.cutoff) & (estimate <= -fc_cutoff) ~ "down",
      TRUE ~ "")  ) %>%
    arrange(p.adjusted)
}

# provide results tables for set A and set B, return overlap gene list
overlap <- function(resultsA, resultsB){

  cnames <- c("gene", "p.adjusted", "p.value")

  if( !all(cnames %in% names(resultsA)) ){
    stop(paste("overlap requires columns:", paste(cnames, collapse = ", ")))
  }
    if( !all(cnames %in% names(resultsB)) ){
    stop(paste("overlap requires columns:", paste(cnames, collapse = ", ")))
    }

  one <- intersect(
    filter(resultsA, (p.adjusted < p.adj.cutoff) & (abs(estimate) >= fc_cutoff) )$gene ,
    filter(resultsB, (p.value < p.adj.cutoff) & (abs(estimate) >= fc_cutoff) )$gene  )

   two <- intersect(
    filter(resultsB, (p.adjusted < p.adj.cutoff) & (abs(estimate) >= fc_cutoff) )$gene ,
    filter(resultsA, (p.value < p.adj.cutoff) & (abs(estimate) >= fc_cutoff) )$gene  )

   union(one, two)
}

to_file <- function(df, name, output_dir){
  p <- file.path(output_dir, paste0(name, ".tsv"))

  df %>%
    dplyr::rename(ensembl_gene_id = gene) %>%
    left_join(ensg_symbol_map, by = "ensembl_gene_id") %>%
    transmute(
      ensembl_gene_id = ensembl_gene_id,
      hgnc_symbol     = hgnc_symbol,
      mean_counts     = baseMean,
      fold_change     = 2^estimate,
      log2_fold_change= estimate,
      lfc_stderror    = stderror,
      statistic       = statistic,
      pvalue	        = p.value,
      padj            = p.adjusted ) %>%
    auto_write(file = p)
  df
}



# requires a tidy results tablewith gene, estimate, p.adjusted columns
# returns a character vector of gene identifier strings
top_genes <- function(df, n = 1, direction = "up"){
  if(direction == "up"){
    df <- filter(df, estimate > 0)
  }else{
    df <- filter(df, estimate < 0)}
  df %>%
    arrange(p.adjusted) %>%
    dplyr::slice(1:n) %>%
    .[["gene"]]
}




ma_plot <- function(df, title){
  ggplot(df, aes(x=baseMean, y = estimate)) +
    geom_point(aes(color = sig), alpha = 0.2) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
    scale_x_log10() +
    labs(x = "mean counts",
         y = "estimated log2 fold change" ,
         title = title) +
    theme( panel.background = element_blank(),
           panel.border = element_rect(fill = NA, color = "black", size = 2),
           legend.position = "none",
           axis.text = element_text(size = 9),
           axis.title = element_text(size = 9)
    )
}


pc_pct_label <- function(var, n){
  paste0("PC",n," (",sprintf("%.1f", var[2,n]),"%)")
}


plot_fpkm <- function(data = data, ensg){
  foo <- filter(data, geneid %in% ensg) %>%
    mutate(trend = paste(platform, repeats))

  grp_foo <- foo %>%
    group_by(platform, repeats, stage, hgnc_symbol) %>%
    summarise(fpkm = mean(fpkm, na.rm=T)) %>%
    ungroup()

  ggplot(foo, aes(x=stage, y = fpkm))+
    geom_point(aes(shape = platform, color = repeats), size = 3) +
    geom_line(data = grp_foo, aes(color = repeats, group = paste0(platform, repeats))) +
    facet_grid(~ hgnc_symbol)  +
    theme_bw() +
    theme(
      axis.text   = element_text(size = 13),
      axis.title  = element_text(size = 15),
      legend.text = element_text(size = 13),
      legend.position = "bottom" )
}
# plot_fpkm(data = data, ensg = hd_specific_genes$both_gene_list[[1]][1:5] )
# plot_fpkm(data = data, ensg = hd_specific_genes$both_gene_list[[2]][1:5] )
# plot_fpkm(data, c( "ENSG00000258017", "ENSG00000171316","ENSG00000134323"))
# plot_fpkm(data, "ENSG00000258017")
# plot_fpkm(data, "ENSG00000100365")
# plot_fpkm(data, "ENSG00000259661")


# this is optimized to be performed as part of a contrast table using pmap
write_heatmap <- function(gene_list, control_group, test_group, data, meta,
                          n=1000, output_dir){

  test_map <- tibble(comparison = c("control group", "test group"),
                     group = c(control_group, test_group))

  if(length(gene_list) < n ){
    n <- length(gene_list) }

  if(n<100){
    show_rownames <- TRUE
  }else{
    show_rownames <- FALSE
  }

  filtered_data <- data %>%
    filter(dataset == "task107") %>%
    filter( geneid %in% gene_list[1:n] ) %>%
    select(geneid, sample, fpkm) %>%
    bio_matrix()

  # remove genes that have no variance, not sure why you'd want to plot anyway ;)
  ind<- apply(filtered_data, 1, sd) == 0
  filtered_data <- filtered_data[!ind,]
  n <- dim(filtered_data)[1]

  gene_ids <- tibble(geneid = row.names(filtered_data)) %>%
  left_join(ensg_symbol_map, by = c("geneid" = "ensembl_gene_id"))

  pheatmap(mat = filtered_data,
           scale = "row",
           color = colorRamps::blue2yellow(10),

           annotation_col = meta %>%
             mutate(group = as.character(group)) %>%
             left_join(test_map, by = "group") %>%
             select(sample, stage, comparison) %>%
             make_rownames("sample"),

           labels_col = task107_meta$fig_label,
           labels_row = gene_ids$hgnc_symbol,

           show_rownames = show_rownames,
           fontsize_row = 9,

           main = paste0("Top ",n," genes from comparison of ",control_group," vs. ",test_group),
           width = 8,
           length = 11,
           filename = file.path(output_dir,  paste0("top_",n,"_genes_from_",control_group,"_vs_",test_group,".pdf"))
           )
  return("")
}

