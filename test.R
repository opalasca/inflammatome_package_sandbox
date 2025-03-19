
#ranked.list <- read_tsv("data/05_agg_list_unionmarkers_latest.tsv",show_col_types = FALSE) 

#all.sets.ens <- read_tsv("data/inflammationGeneSetsEnsembl.tsv",show_col_types = FALSE)
#all.sets.entrez <- read_tsv("data/inflammationGeneSetsEntrez.tsv",show_col_types = FALSE)


########## User input data #########################
####################################################

## Run for different input -------------------------------------------------------------

igan <- read_tsv("data/test_datasets/02_GSE175759_IgAN_ctl.tsv",show_col_types = FALSE) 

data <- process_input_data(igan,id="ENSG.ID", keytype="ensembl")
gene_sets <- get_gene_sets("Ensembl")
gsea_results <- gsea_analysis(data, gene_sets, sorting_value_col_name = "stat", name="igan")
plot_volcano(data, keytype="Ensembl", logFC_col_name = "log2FoldChange", pval_col_name = "pvalue")
plot_volcano_and_stripplot(data, keytype="Ensembl", logFC_col_name = "log2FoldChange", pval_col_name = "pvalue", stat_col_name = "stat")


ms <- read_tsv("data/test_datasets/02_GSE138614_MS_CTL.tsv", show_col_types = FALSE)

data <- process_input_data(ms,id="ENSG.ID",keytype="ensembl")
gene_sets <- get_gene_sets("Ensembl")
gsea_results <- gsea_analysis(data, gene_sets, sorting_value_col_name = "stat", name="ms")
plot_volcano(data, logFC_col_name = "log2FoldChange", pval_col_name = "pvalue", name="ms",keytype="Ensembl")

UC.Andersen <- read.delim("data/test_datasets/DE.res.UC.Andersen.raw.tsv", header = TRUE, sep = "\t")

data <- process_input_data(UC.Andersen,id="id",keytype="Uniprot")
gene_sets <- get_gene_sets("Uniprot")
gsea_results <- gsea_analysis(data, gene_sets, sorting_value_col_name = "t", name="data")
plot_volcano(data, logFC_col_name = "logFC", pval_col_name = "P.Value", name="UC.proteomics",keytype="Entrez")



install.packages("shiny")
library(shiny)
runGitHub("inflammatome_package_sandbox", "opalasca")
runApp("app")


## More datasets to test ------------------------------------------------------------------

UC.hansen.pre <- read_excel("data/proteomics/UC_Hansen/supp.table.processed.xlsx", sheet = "All proteins")

### UC Hansen data processing --------------------------------------------------
# filter proteins in over 70% samples (they do this in the paper)
UC.hansen.pre <- filter(UC.hansen.pre, Number.of.samples.quantified.in >= 23)

# doing this to avoid weird duplicates
multiple.genes <- filter(UC.hansen.pre, str_detect(Gene.names, ";"))
single.gene <- filter(UC.hansen.pre, !(Gene.names %in% multiple.genes$Gene.names))

multiple.genes <- multiple.genes %>%  
  mutate(Gene.name = Gene.names) %>%
  separate_rows(Gene.name, sep = ";") %>%
  filter(!(Gene.name %in% single.gene$Gene.names))

single.gene <- single.gene %>%
  mutate(Gene.name = Gene.names)

UC.hansen <- rbind(multiple.genes, single.gene) %>%
  filter(Gene.name %in% all.genes$Gene.name)  %>% # after this there are still some duplicate Gene.name, then just keep the one with lowest pval like Oana did
  mutate(pvalue = 10 ^ (-minus.log.pvalue),
         adj_pvalue = p.adjust(pvalue, method = "BH"), 
         logFC = log2(Ratio.UC.H),
         sorting.value = logFC*-log(adj_pvalue)) %>%
  drop_na(Gene.name) %>% 
  group_by(Gene.name) %>% 
  filter(pvalue == min(pvalue)) %>%
  distinct(Gene.name, .keep_all = TRUE) %>% # Keep first occurrence bc there is one weird dup with same pval 
  arrange(desc(sorting.value)) %>%
  ungroup()
# create sorting value:  logFC*-log(pval) 

UC.hansen <- UC.hansen %>%
  rename(stat = sorting.value) %>%
  inner_join(Final_Annotation_List)




