# Clear workspace --------------------------------------------------------------
#rm(list = ls())
dir="~/Desktop/Work/inflammatome/inflammatome_resource"
setwd(dir)

# Install packages (if not already present) and Load libraries ---------------------------------------------------------------
list.of.packages <- c("ggplot2","ggrepel","msigdbr","tidyverse","GO.db","org.Hs.eg.db","readxl","biomaRt","clusterProfiler","enrichplot","tidytext","dplyr","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, library, character.only=TRUE)

BiocManager::install("biomaRt")
remove(biomaRt)

# Define functions -------------------------------------------------------------
source(file = "scripts/99_project_functions.R")

# Read inflammatome list  --------------------------------------------------------------------
#all.genes <- read_tsv("data/04_rank_agg_list.tsv") 
all.genes <- read_tsv("data/05_agg_list_unionmarkers_latest.tsv",show_col_types = FALSE) 

all.sets <- read_tsv("data/inflammationGeneSets.tsv",show_col_types = FALSE)

################################################################################
######### Function for GSEA, GSEA dotplot and volcano plot #####################
################################################################################

process_data <- function(data, id){
  #  Put here code for pre-processing - mapping ids, etc.  
  #data[[id]] <- as.character(data[[id]])
  return(data)
}

run_gsea <- function(results_list, gene_sets, id_col_name, sorting_value_col_name){
  results_list <- results_list[order(results_list[[sorting_value_col_name]], decreasing = TRUE), ]
  ranked_gene_list <- setNames(results_list[[sorting_value_col_name]], results_list[[id_col_name]])
  # Run GSEA on the combined gene sets
  set.seed(100)
  gsea_results <- GSEA(
    ranked_gene_list, 
    TERM2GENE = gene_sets, 
    minGSSize = 0, 
    maxGSSize = 2000, 
    pvalueCutoff = 10,
    eps = 0
  )
  return(gsea_results)
}

plot_gsea <- function(gsea_results, name){
  
  gsea_results$Description <- as.character(gsea_results$Description)
  gsea_results$Description <- gsub("_"," ",gsea_results$Description)
  gsea_results$Description <- tolower(gsea_results$Description)
  gsea_results$Description <- sub("^(.)", "\\U\\1", gsea_results$Description, perl = TRUE)  # Capitalize first letter
  gsea_results$Description <- gsub("Gobp","GOBP",gsea_results$Description)
  gsea_results$Description <- gsub("Wp","WP",gsea_results$Description)
  
  # Wrap text at 32 characters
  gsea_results$Description <- str_wrap(gsea_results$Description, width = 35)
  
  p <- ggplot(gsea_results, aes( x=NES, y = reorder(Description,NES))) + 
    geom_point(aes(fill = NES, size = -log10(p.adjust)), shape = 21, colour = "black",alpha = 0.8) +
    scale_fill_gradient2(midpoint = 0,low = "blue4",mid = "white",high = "red4",space = "Lab") +
    theme_light() + 
    ylab(NULL) + xlab("NES") 

  # Add black dots for p.adjust < 0.05
  p <- p + geom_point(data = gsea_results[gsea_results$p.adjust < 0.05,], 
                      aes(x = NES, y = Description), 
                      color = "black", size = 1)
  print(p)
  ggsave(paste0("figures/gsea_",name,".png"), p, h = 4.5, w = 5.8)
}


data=igan

plot_volcano <- function(data, logFC_col_name="log2FoldChange", pval_col_name="pvalue", name="data") {
    
  theme_custom <- theme_classic() + theme(
      plot.title = element_text(size = 10))
  
  data <- data %>%
    arrange(desc(.data[[pval_col_name]])) %>%
    mutate(idx = row_number()) %>%
    mutate(
      top_2000 = if_else(ENSG.ID %in% filter(all.genes, Position <= 2000)$ENSG.ID, "yes", "no"),
      top_100 = if_else(ENSG.ID %in% filter(all.genes, Position <= 100)$ENSG.ID, "yes", "no")
    ) %>%
    mutate(
      top_2000 = factor(top_2000, levels = c("no", "yes")),  # Ensure "no" is plotted first
      top_100 = factor(top_100, levels = c("no", "yes"))
    )  %>%
    mutate(
      # **Cap -log10(p-value) at 20**
      adj_pval = pmin(-log10(.data[[pval_col_name]]), 20),
      # **Cap log fold change between -10 and 10**
      adj_logFC = pmax(pmin(.data[[logFC_col_name]], 10), -10)
    )
  
  p1 <- ggplot(data %>% arrange(top_2000), aes(x = adj_logFC, y = adj_pval, color = top_2000, alpha = top_2000)) +
    geom_point(size = .4, alpha=0.7) +
    scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
    scale_alpha_manual(values = c("no" = 0.5, "yes" = 0.5)) +
    labs(
      x = "Log Fold Change",
      y = "-log10(P-Value)",
      title = paste0(name,", inflammatome (top 2000)")
    ) +
    theme_custom
    #coord_cartesian(ylim = c(0, 20))
  
  print(p1)
  ggsave(paste0("figures/volcano_2000_",name,".png"), p1, height = 3, width = 3)
  
  p2 <- ggplot(data %>% arrange(top_100), aes(x = adj_logFC, y = adj_pval, color = top_100, alpha = top_100)) +
    geom_point(size = .5) +
    scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
    scale_alpha_manual(values = c("no" = 0.3, "yes" = 1)) +
    labs(
      x = "Log Fold Change",
      y = "-log10(P-Value)",
      title = paste0(name,", inflammation signature (top 100)")
      ) +
    theme_custom
    #coord_cartesian(ylim = c(0, 20))
  
  print(p2)
  ggsave(paste0("figures/volcano_100_",name,".png"), p2, height = 3, width = 3)
  
  combined_plot <- p1 + p2 + plot_layout(ncol = 2,guides = "collect")
  
  # Print and save the combined plot
  print(combined_plot)
  ggsave(paste0("figures/volcano_",name,".png"), combined_plot, height = 3, width = 6)  # Wider panel
  
}



########## User input data #########################
####################################################


## Run for different input -------------------------------------------------------------

igan <- read_tsv("data/test_datasets/02_GSE175759_IgAN_ctl.tsv",show_col_types = FALSE) 

data <- process_data(igan,id="ENSG.ID")
gsea <- run_gsea(results_list=data, gene_sets=all.sets, id_col_name="ENSG.ID", sorting_value_col_name = "stat")
plot_gsea(gsea@result, "igan")
plot_volcano(data, logFC_col_name = "log2FoldChange", pval_col_name = "pvalue", name="igan")


ms <- read_tsv("data/test_datasets/02_GSE138614_MS_CTL.tsv", show_col_types = FALSE)

data <- process_data(ms,id="ENSG.ID")
gsea <- run_gsea(results_list=data, gene_sets=all.sets, id_col_name="ENSG.ID", sorting_value_col_name = "stat")
plot_gsea(gsea@result, "ms")
plot_volcano(data, logFC_col_name = "log2FoldChange", pval_col_name = "pvalue", name="ms")

UC.Andersen <- read.delim("data/test_datasets/DE.res.UC.Andersen.raw.tsv", row.names = 1, header = TRUE, sep = "\t")
UC.Andersen$id = rownames(UC.Andersen)



