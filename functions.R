# Remove setwd line! 
dir="~/Desktop/Work/inflammatome/inflammatome_resource"
setwd(dir)


#figures=paste0(dir,"/figures/")
#resultdir=paste0(dir,"results/")
#ifelse(!dir.exists(figures), dir.create(figures), FALSE)
#ifelse(!dir.exists(resultdir), dir.create(resultdir), FALSE)

# Install packages (if not already present) and Load libraries ---------------------------------------------------------------
list.of.packages <- c("here","fs","ggplot2","ggrepel","patchwork","msigdbr","tidyverse","GO.db","org.Hs.eg.db","readxl",
                      "biomaRt","clusterProfiler","enrichplot","tidytext","dplyr","stringr", "DOSE")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, library, character.only=TRUE)

here::i_am("functions.R")
figures_dir <- here("figures")
fs::dir_create(figures_dir)

# Read inflammatome list -------------------------------
data_file <- here("data", "ranked_list_inflammatome.tsv")
if (file_exists(data_file)) {
  ranked.list <- read_tsv(data_file, show_col_types = FALSE)
} else {
  stop("Error: Data file not found. Ensure 'data/ranked_list_inflammatome.tsv' exists.")
}


# Convert gene names to Entrez IDs ------
symbol.entrez = bitr(ranked.list$Gene.name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ranked.list.entrez.symbol = merge(ranked.list, symbol.entrez, by.x = "Gene.name", by.y = "SYMBOL", all.x = TRUE)

# Get top 2000 genes (inflammatome) and top 100 (inflammation signature) -------
top2000 <- ranked.list.entrez.symbol %>% 
  filter(Position <= 2000)
top100 <- ranked.list.entrez.symbol %>% 
  filter(Position <= 100)



################################################################################
######### Function for GSEA, GSEA dotplot and volcano plot #####################
################################################################################


# Preprocess input file and prepare gene sets for GSEA - not used yet
prepare_data <- function(data, id_col_name, keytype, mean_expr_col_name=NULL, pval_col_name=NULL){
  
  data <- process_input_data(data, id_col_name, keytype, mean_expr_col_name, pval_col_name)
  gene_sets <- get_gene_sets(keytype)
  
  return(list(data, gene_sets))
}


# Preprocess input files  
process_input_data <- function(data, id_col_name, keytype, mean_expr_col_name=NULL, pval_col_name=NULL){
  
  keytype <- tolower(keytype)
  
  if (keytype == "uniprot") {
    
    # TO DO: Check if duplicates and give warning that there are duplicated ids
  
    id.mapping <- bitr(data[[id_col_name]], 
                       fromType = "UNIPROT", 
                       toType = "ENTREZID", 
                       OrgDb = org.Hs.eg.db)
    
    # Merge the mapped IDs with the original data
    data <- merge(data, id.mapping, by.x=id_col_name, by.y = "UNIPROT")
    
    # Ensure ENTREZID is unique by keeping the row with the smallest UNIPROT alphabetically
    data <- data %>%
      group_by(ENTREZID) %>%
      arrange(!!sym(id_col_name)) %>%  # Sort by UniProt ID alphabetically
      slice(1) %>%  # Keep only the first (smallest UniProt)
      ungroup()
    
    data$id.mapped = data$ENTREZID
    
  }
  
  if (keytype == "symbol") {
    
    # TO DO: Check if duplicates and give warning that there are duplicated ids
    
    id.mapping <- bitr(data[[id_col_name]], 
                       fromType = "SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = org.Hs.eg.db)
    
    # Merge the mapped IDs with the original data
    data <- merge(data, id.mapping, by.x=id_col_name, by.y = "UNIPROT")
    
    # Ensure ENTREZID is unique by keeping the row with the smallest UNIPROT alphabetically
    data <- data %>%
      group_by(ENTREZID) %>%
      arrange(!!sym(id_col_name)) %>%  # Sort by UniProt ID alphabetically
      slice(1) %>%  # Keep only the first (smallest UniProt)
      ungroup()
    
    data$id.mapped = data$ENTREZID
    
  }
  
  if (keytype == "symbol") {
    
    # TO DO: Check if duplicates and give warning that there are duplicated ids
    
    id.mapping <- bitr(data[[id_col_name]], 
                       fromType = "SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = org.Hs.eg.db)
    
    # Merge the mapped IDs with the original data
    data <- merge(data, id.mapping, by.x=id_col_name, by.y = "UNIPROT")
    
    # Ensure ENTREZID is unique by keeping the row with the smallest UNIPROT alphabetically
    data <- data %>%
      group_by(ENTREZID) %>%
      arrange(!!sym(id_col_name)) %>%  # Sort by UniProt ID alphabetically
      slice(1) %>%  # Keep only the first (smallest UniProt)
      ungroup()
    
    data$id.mapped = data$ENTREZID
    
  }
  
  if (keytype == "refseq") {
    
    # TO DO: Check if duplicates and give warning that there are duplicated ids
    
    id.mapping <- bitr(data[[id_col_name]], 
                       fromType = "REFSEQ", 
                       toType = "ENTREZID", 
                       OrgDb = org.Hs.eg.db)
    
    # Merge the mapped IDs with the original data
    data <- merge(data, id.mapping, by.x=id_col_name, by.y = "UNIPROT")
    
    # Ensure ENTREZID is unique by keeping the row with the smallest id alphabetically
    data <- data %>%
      group_by(ENTREZID) %>%
      arrange(!!sym(id_col_name)) %>%  # Sort by UniProt ID alphabetically
      slice(1) %>%  # Keep only the first (smallest UniProt)
      ungroup()
    
    data$id.mapped = data$ENTREZID
    
  }
  
  
  if (keytype == "entrez") {
    
    # TO DO: Check if duplicates and give warning that there are duplicated ids
    
    data$id.mapped = data[[id_col_name]]
    
  }
  
  if (keytype == "ensembl") {
    #check and give warning if duplicate gene names
    data$id.mapped = data[[id_col_name]]
  }
  
  return(data)
  
}

# Extract gene sets for GSEA
get_gene_sets <- function(keytype){ 
  
  keytype <- tolower(keytype)
  
  # if the keytype is ensembl, extract gene sets using ensembl ids
  if (keytype == "ensembl") {
    
  go <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
    dplyr::filter(gs_exact_source %in% c("GO:0002534", "GO:0006954", "GO:0002526", "GO:0002544")) %>%
    dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_exact_source) %>%
    dplyr::rename(gene = ensembl_gene)
  
  wp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%
    dplyr::filter(gs_exact_source %in% c("WP4493", "WP530", "WP5198", "WP453")) %>%
    dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_exact_source) %>%
    dplyr::rename(gene = ensembl_gene)
  
  reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
    dplyr::filter(gs_exact_source %in% c("R-HSA-622312", "R-HSA-913531", "R-HSA-9020702", "R-HSA-75893")) %>%
    dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_exact_source) %>%
    dplyr::rename(gene = ensembl_gene)
  
  h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, ensembl_gene, entrez_gene, gene_symbol)
  
  msigdb.markers <- h[h$gs_name == "HALLMARK_INFLAMMATORY_RESPONSE",]
  #msigdb.markers = msigdb.markers[c(2,3,4)]
  
  top.100 = filter(ranked.list, Position <= 100)$ENSG.ID
  go = dplyr::select(go, gs_name, gene)
  wp = dplyr::select(wp, gs_name, gene)
  reactome = dplyr::select(reactome, gs_name, gene)
  msigdb = msigdb.markers$ensembl_gene
  
  all.sets <- rbind(
    data.frame(gs_name = "Inflammation signature (top100)", gene = top.100) ,
    data.frame(gs_name = "MSigDB hallmark inflammatory response", gene = msigdb),
    go,
    wp,
    reactome
  )
  }
  
  # Use Entrez if the original identifier is not Ensembl
  if (keytype != "ensembl") {
    go <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
      dplyr::filter(gs_exact_source %in% c("GO:0002534", "GO:0006954", "GO:0002526", "GO:0002544")) %>%
      dplyr::select(gs_name, entrez_gene, gene_symbol, gs_exact_source) %>%
      dplyr::rename(gene = entrez_gene)
    
    wp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%
      dplyr::filter(gs_exact_source %in% c("WP4493", "WP530", "WP5198", "WP453")) %>%
      dplyr::select(gs_name, entrez_gene, gene_symbol, gs_exact_source) %>%
      dplyr::rename(gene = entrez_gene)
    
    reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
      dplyr::filter(gs_exact_source %in% c("R-HSA-622312", "R-HSA-913531", "R-HSA-9020702", "R-HSA-75893")) %>%
      dplyr::select(gs_name, entrez_gene, gene_symbol, gs_exact_source) %>%
      dplyr::rename(gene = entrez_gene)
    
    h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
      dplyr::select(gs_name, entrez_gene, entrez_gene, gene_symbol)
    msigdb.markers <- h[h$gs_name == "HALLMARK_INFLAMMATORY_RESPONSE",]
    #msigdb.markers = msigdb.markers[c(2,3,4)]
    
    top.100 = filter(ranked.list.entrez.symbol, Position <= 100)$ENTREZID
    go = dplyr::select(go, gs_name, gene)
    wp = dplyr::select(wp, gs_name, gene)
    reactome = dplyr::select(reactome, gs_name, gene)
    msigdb = msigdb.markers$entrez_gene
    
    all.sets <- rbind(
      data.frame(gs_name = "Inflammation signature (top100)", gene = top.100) ,
      data.frame(gs_name = "MSigDB hallmark inflammatory response", gene = msigdb),
      go,
      wp,
      reactome
    )
  }
  
  return(all.sets)
}


# Run GSEA and create a dotplot 
gsea_analysis <- function(data, gene_sets, sorting_value_col_name, name="data"){
  gsea_results <- run_gsea(data, gene_sets, sorting_value_col_name)
  gsea_results_df <- gsea_results@result
  plot_gsea(gsea_results_df, name)
  return(gsea_results_df)
}


# Run GSEA
run_gsea <- function(data, gene_sets, sorting_value_col_name){
  data <- data[order(data[[sorting_value_col_name]], decreasing = TRUE), ]
  ranked_gene_list <- setNames(data[[sorting_value_col_name]], data$id.mapped)
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


# Dotplot for GSEA results
plot_gsea <- function(gsea_result_df, name){
  
  gsea_result_df$Description <- as.character(gsea_result_df$Description)
  gsea_result_df$Description <- gsub("_"," ",gsea_result_df$Description)
  gsea_result_df$Description <- tolower(gsea_result_df$Description)
  gsea_result_df$Description <- sub("^(.)", "\\U\\1", gsea_result_df$Description, perl = TRUE)  # Capitalize first letter
  gsea_result_df$Description <- gsub("Gobp","GOBP",gsea_result_df$Description)
  gsea_result_df$Description <- gsub("Wp","WP",gsea_result_df$Description)
  
  # Wrap text at 32 characters
  gsea_result_df$Description <- str_wrap(gsea_result_df$Description, width = 35)
  
  p <- ggplot(gsea_result_df, aes( x=NES, y = reorder(Description,NES))) + 
    geom_point(aes(fill = NES, size = -log10(p.adjust)), shape = 21, colour = "black",alpha = 0.8) +
    scale_fill_gradient2(midpoint = 0,low = "blue4",mid = "white",high = "red4",space = "Lab") +
    theme_light() + 
    ylab(NULL) + xlab("NES") 

  # Add black dots for p.adjust < 0.05
  p <- p + geom_point(data = gsea_result_df[gsea_result_df$p.adjust < 0.05,], 
                      aes(x = NES, y = Description), 
                      color = "black", size = 1)
  print(p)
  ggsave(paste0("figures/gsea_",name,".png"), p, h = 4.5, w = 5.8)
}


# Volcano plot
plot_volcano <- function(data, keytype="Ensembl", logFC_col_name="log2FoldChange", pval_col_name="pvalue", name="data"){ 
  
  theme_custom <- theme_classic() + theme(
      plot.title = element_text(size = 8))
  
  keytype <- tolower(keytype)
  
  if (keytype == "ensembl") {
  data <- data %>%
    arrange(desc(.data[[pval_col_name]])) %>%
    mutate(idx = row_number()) %>%
    mutate(
      top_2000 = if_else(id.mapped %in% filter(top2000, Position <= 2000)$ENSG.ID, "yes", "no"),
      top_100 = if_else(id.mapped %in% filter(top2000, Position <= 100)$ENSG.ID, "yes", "no")
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
    )}
  #pval_col_name="P.Value"; logFC_col_name="logFC"; name="UC.proteomics" 
  else{
      data <- data %>%
      arrange(desc(.data[[pval_col_name]])) %>%
      mutate(idx = row_number()) %>%
      mutate(
        top_2000 = if_else(id.mapped %in% filter(top2000, Position <= 2000)$ENTREZID, "yes", "no"),
        top_100 = if_else(id.mapped %in% filter(top2000, Position <= 100)$ENTREZID, "yes", "no")
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
      )  }
  
  p1 <- ggplot(data %>% arrange(top_2000), aes(x = adj_logFC, y = adj_pval, color = top_2000 ))+ # , alpha = top_2000)) +
    geom_point(size = .4, alpha=0.6) +
    scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
    #scale_alpha_manual(values = c("no" = 0.5, "yes" = 0.5)) +
    labs(
      x = "Log Fold Change",
      y = "-log10(P-Value)",
      #title = paste0(name,", inflammatome (top 2000)")
      title = paste0("inflammatome (top 2000)")
    ) +
    theme_custom
    #coord_cartesian(ylim = c(0, 20))
  
  print(p1)
  ggsave(paste0("figures/volcano_2000_",name,".png"), p1, height = 3, width = 3)
  
  p2 <- ggplot(data %>% arrange(top_100), aes(x = adj_logFC, y = adj_pval, color = top_100))+ #, alpha = top_100)) +
    geom_point(size = .4, alpha=0.6) +
    scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
    labs(
      x = "Log Fold Change",
      y = "-log10(P-Value)",
      #title = paste0(name,", inflammation signature (top 100)")
      title = paste0("inflammation signature (top 100)")
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

plot_volcano_and_stripplot <- function(data,  keytype="Ensembl", logFC_col_name="log2FoldChange", pval_col_name="pvalue", stat_col_name="stat"){ 
  
  theme_custom <- theme_classic() + theme(
    plot.title = element_text(size = 8))
  
  theme_custom_strip_plot <- theme_classic() + theme(
    plot.title = element_blank(),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.line.y = element_blank()  # Add this line to hide the y-axis line
    
  )
  
  keytype <- tolower(keytype)
  
  if (keytype == "ensembl") {
    data <- data %>%
      arrange(desc(.data[[pval_col_name]])) %>%
      mutate(idx = row_number()) %>%
      mutate(
        top_2000 = if_else(id.mapped %in% filter(top2000, Position <= 2000)$ENSG.ID, "yes", "no"),
        top_100 = if_else(id.mapped %in% filter(top2000, Position <= 100)$ENSG.ID, "yes", "no")
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
      )}
  #pval_col_name="P.Value"; logFC_col_name="logFC"; name="UC.proteomics" 
  else{
    data <- data %>%
      arrange(desc(.data[[pval_col_name]])) %>%
      mutate(idx = row_number()) %>%
      mutate(
        top_2000 = if_else(id.mapped %in% filter(top2000, Position <= 2000)$ENTREZID, "yes", "no"),
        top_100 = if_else(id.mapped %in% filter(top2000, Position <= 100)$ENTREZID, "yes", "no")
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
      )  }
  
  p1 <- ggplot(data %>% arrange(top_2000), aes(x = adj_logFC, y = adj_pval, color = top_2000 ))+ # , alpha = top_2000)) +
    geom_point(size = .4, alpha=0.6) +
    scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
    #scale_alpha_manual(values = c("no" = 0.5, "yes" = 0.5)) +
    labs(
      x = logFC_col_name,
      y = "-log10(p-value)",
      #title = paste0(name,", inflammatome (top 2000)")
      title = paste0("inflammatome (top 2000)")
    ) +
    theme_custom

  #print(p1)
  ggsave(paste0("figures/volcano_2000",".png"), p1, height = 3, width = 3)
  
  p2 <- ggplot(data %>% arrange(top_100), aes(x = adj_logFC, y = adj_pval, color = top_100))+ #, alpha = top_100)) +
    geom_point(size = .4, alpha=0.6) +
    scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
    labs(
      x = logFC_col_name,
      y = "-log10(p-value)",
      #title = paste0(name,", inflammation signature (top 100)")
      title = paste0("inflammation signature (top 100)")
    ) +
    theme_custom

  #print(p2)
  ggsave(paste0("figures/volcano_100",".png"), p2, height = 3, width = 3)
  
  max_abs_stat <- max(abs(data[[stat_col_name]]), na.rm = TRUE)
  
  p3 <- ggplot(data, aes(y = factor(1), x = .data[[stat_col_name]], color = top_2000)) +
    geom_jitter(width = 0, height = 0.2, size = 0.4) +
    scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
    geom_vline(xintercept = 0, color = "darkgrey") +
    scale_x_continuous(limits = c(-max_abs_stat, max_abs_stat)) +
    labs(
      x = stat_col_name,
      y = "",
      title = ""
    ) +
    theme_custom_strip_plot
  
  p4 <- ggplot(data, aes(y = factor(1), x = .data[[stat_col_name]], color = top_100)) +
    geom_jitter(width = 0, height = 0.2, size = 0.4) +
    scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
    geom_vline(xintercept = 0, color = "darkgrey") +
    scale_x_continuous(limits = c(-max_abs_stat, max_abs_stat)) +
    labs(
      x = stat_col_name,
      y = "",
      title = ""
    ) +
    theme_custom_strip_plot
  
  #print(p3)
  ggsave(paste0("figures/strip_plot_2000",".png"), p3, height = 3, width = 3)
  
  #print(p4)
  ggsave(paste0("figures/strip_plot_100",".png"), p4, height = 3, width = 3)
  
  combined_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2, guides = "collect", heights = c(2, 1))
  
  # Print and save the combined plot
  print(combined_plot)
  ggsave(paste0("figures/volcano",".png"), combined_plot, height = 3, width = 6)  # Wider panel
  
}





