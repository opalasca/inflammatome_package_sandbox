
dir="~/Desktop/Work/inflammatome/inflammatome_resource"
setwd(dir)

# Install packages (if not already present) and Load libraries ---------------------------------------------------------------
list.of.packages <- c("ggplot2","ggrepel","msigdbr","tidyverse","GO.db","org.Hs.eg.db","readxl","biomaRt","clusterProfiler","enrichplot","tidytext","dplyr","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, library, character.only=TRUE)



# Define functions -------------------------------------------------------------
#source(file = "scripts/99_project_functions.R")

# Read inflammatome list  --------------------------------------------------------------------
all.genes <- read_tsv("data/05_agg_list_unionmarkers_latest.tsv",show_col_types = FALSE) 


########## Prepare inflammation gene sets ############
######################################################

## GO sets ---------------------------------------------------------------------

go <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_exact_source) %>%
  filter(gs_exact_source %in% c("GO:0002534", "GO:0006954", "GO:0002526", "GO:0002544")) %>%
  dplyr::rename(gene = ensembl_gene,
         gs = gs_name)
length(intersect(go$gene[go$gs == "GOBP_INFLAMMATORY_RESPONSE"], filter(all.genes, Position <= 2000)$ENSG.ID))
length(unique(filter(go, gs == "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE")$gene_symbol))

setdiff(filter(go, gs == "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE")$gene, Final_Annotation_List$ENSG.ID)
# some duplicate gene symbols with different ensembl id because they are not in chromosome
setdiff(filter(go, gs == "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE")$gene_symbol, Final_Annotation_List$Gene.name)
# there are also some non-coding RNAs

a<-go %>% 
  group_by(gs) %>%
  dplyr::count(gene)


## Wikipathways ----------------------------------------------------------------
# WP4493 WP_CELLS_AND_MOLECULES_INVOLVED_IN_LOCAL_ACUTE_INFLAMMATORY_RESPONSE
# WP530 WP_CYTOKINES_AND_INFLAMMATORY_RESPONSE
# WP5198 WP_INFLAMMATORY_BOWEL_DISEASE_SIGNALING ?? not found
# WP453 WP_INFLAMMATORY_RESPONSE_PATHWAY

wp <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_exact_source) %>%
  filter(gs_exact_source %in% c("WP4493", "WP530", "WP5198", "WP453")) %>%
  dplyr::rename(gene = ensembl_gene,
         gs = gs_name)

wp %>% 
  dplyr::count(gs)

a<-wp %>% 
  group_by(gs) %>%
  dplyr::count(gene) 


## Reactome --------------------------------------------------------------------
# R-HSA-622312 REACTOME_INFLAMMASOMES
# R-HSA-913531 REACTOME_INTERFERON_SIGNALING
# R-HSA-9020702 REACTOME_INTERLEUKIN_1_SIGNALING
# R-HSA-75893 REACTOME_TNF_SIGNALING

reactome <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_exact_source) %>%
  filter(gs_exact_source %in% c("R-HSA-622312", "R-HSA-913531", "R-HSA-9020702", "R-HSA-75893")) %>%
  dplyr::rename(gene = ensembl_gene,
         gs = gs_name)

reactome %>% 
  dplyr::count(gs)

a<-reactome %>% 
  group_by(gs) %>%
  dplyr::count(gene) 

setdiff(unique(reactome$gene), Final_Annotation_List$ENSG.ID) 



## Bind all sets ---------------------------------------------------------------
top.100 = filter(all.genes, Position <= 100)$ENSG.ID
go = dplyr::select(go, gs, gene)
wp = dplyr::select(wp, gs, gene)
reactome = dplyr::select(reactome, gs, gene)
msigdb = msigdb.markers$ENSG.ID

all.sets <- rbind(
  data.frame(gs = "Inflammation signature (top100)", gene = top.100) ,
  data.frame(gs = "MSigDB hallmark inflammatory response", gene = msigdb),
  go,
  wp,
  reactome
)

all.sets %>% dplyr::count(gs)
all.sets %>% dplyr::count(gs) %>% 
  ggplot(aes(n)) + 
  geom_histogram() +
  theme_bw() +
  labs(title = "Distribution of gene set sizes")

## Save file -------- 
write_tsv(all.sets, "data/inflammationGeneSets.tsv") 

#all.sets <- read_tsv("data/inflammationGeneSets.tsv")




#

########## Extract biomart  ##########################
######################################################

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
a=ensembl@attributes

genes.uniprot.biomart <- getBM(attributes = c("uniprotswissprot", "ensembl_gene_id", "hgnc_symbol","gene_biotype"), mart = ensembl)
genes.refseq.biomart <- getBM(attributes = c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol","gene_biotype"), mart = ensembl)
genes.refseq.ncrna.biomart <- getBM(attributes = c("refseq_ncrna", "ensembl_gene_id", "hgnc_symbol","gene_biotype"), mart = ensembl)
genes.entrez.biomart <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "hgnc_symbol","gene_biotype"), mart = ensembl)
genes.uniprot.refseq.biomart <- getBM(attributes = c("uniprotswissprot", "ensembl_gene_id", "refseq_mrna"), mart = ensembl)


write_tsv(genes.uniprot, "data/genes.uniprot.tsv") 
write_tsv(genes.refseq, "data/genes.refseq.tsv") 
write_tsv(genes.entrez, "data/genes.entrez.tsv") 

genes.uniprot.biomart <- read_tsv("data/genes.uniprot.tsv",show_col_types = FALSE) 
genes.refseq.biomart <- read_tsv("data/genes.refseq.tsv",show_col_types = FALSE)
genes.entrez.biomart <- read_tsv("data/genes.entrez.tsv",show_col_types = FALSE)

genes.uniprot = genes.uniprot.biomart[!is.na(genes.uniprot.biomart$uniprotswissprot),]
genes.refseq = genes.refseq.biomart[!is.na(genes.refseq.biomart$refseq_mrna),]


genes.uniprot <- read_tsv("data/genes.uniprot.tsv",show_col_types = FALSE) 
genes.refseq <- read_tsv("data/genes.refseq.tsv",show_col_types = FALSE)
genes.entrez <- read_tsv("data/genes.entrez.tsv",show_col_types = FALSE)



length(unique(genes.refseq$refseq_mrna))
length(unique(genes.refseq$ensembl_gene_id))
length(unique(genes.uniprot$ensembl_gene_id))
length(unique(genes.uniprot$uniprotswissprot))

#write.table(genes.ensembl, paste(datadir,"genes.ensembl.UC.Andersen.tsv",sep=""), sep = "\t", row.names = F, quote=F)
genes.ensembl.in.list <- genes.ensembl[genes.ensembl$ensembl_gene_id %in% ranked.list$ENSG.ID, ]

# Remove duplicates; for each uniprot ID keep the first ensembl ID, alphabetically
genes.ensembl.in.list <- genes.ensembl.in.list[order(genes.ensembl.in.list$uniprotswissprot, genes.ensembl.in.list$ensembl_gene_id), ]
genes.ensembl.in.list.dedup <- genes.ensembl.in.list[!duplicated(genes.ensembl.in.list$uniprotswissprot), ]

df.prot.dedup=merge(df.prot, genes.ensembl.in.list.dedup[c(1,2)],by.x="UniprotID",by.y="uniprotswissprot")
length(unique(df.prot.dedup$UniprotID))










## UC Andersen -----------------------------------------------------------------
UC.andersen <- UC.andersen %>% 
  group_by(ENSG.ID) %>% 
  filter (P.Value == min(P.Value))

gsea_UC.andersen <- run_gsea(UC.andersen, all.sets, "stat")
p_andersen <- dotplot(gsea_UC.andersen, x = "NES", size= "GeneRatio",showCategory = 13) + ggtitle("UC Andersen stat value") 
p_andersen

gsea_UC.andersen.lfc <- run_gsea(UC.andersen, all.sets, "logFC")
dotplot(gsea_UC.andersen.lfc, x = "NES", size= "GeneRatio",showCategory = 13) + ggtitle("UC Andersen LFC") 

ggsave("figures/05_UC_andersen_dotplot.png", 
       p_andersen,
       h = 9,
       w = 10)






# Volcano plots ----------------------------------------------------------------
## IgA neph --------------------------------------------------------------------
igan <- igan %>%
  arrange(desc(stat)) %>%
  mutate(idx = 1:length(stat)) %>%
  mutate(top_2000 = if_else(ENSG.ID %in% filter(all.genes, Position <= 2000)$ENSG.ID, "yes", "no"),
         top_100 = if_else(ENSG.ID %in% filter(all.genes, Position <= 100)$ENSG.ID, "yes", "no"))

ggplot(igan, aes(x = log2FoldChange, y = -log10(pvalue), color = top_2000, alpha = top_2000)) +
#ggplot(igan, aes(x = log2FoldChange, y = -log10(padj), color = top_2000, alpha = top_2000)) +
  geom_point(size = .3) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.3, "yes" = 1)) +
  labs(
    x = "Log Fold Change",
    y = "-log10(P-Value)",
    title = "IgA nephropathy"
  ) +
  theme_minimal() +
 coord_cartesian(ylim = c(0,20))

## UC Andersen ------------------------------------------------------------------
UC.andersen <- UC.andersen %>%
  arrange(desc(stat)) %>%
  mutate(idx = 1:length(stat)) %>%
  mutate(top_2000 = if_else(ENSG.ID %in% filter(all.genes, Position <= 2000)$ENSG.ID, "yes", "no"),
         top_100 = if_else(ENSG.ID %in% filter(all.genes, Position <= 100)$ENSG.ID, "yes", "no"))

#ggplot(UC.andersen, aes(x = logFC, y = -log10(adj.P.Val), color = top_2000, alpha = top_2000)) +
  ggplot(UC.andersen, aes(x = logFC, y = -log10(P.Value), color = top_2000, alpha = top_2000)) +
  geom_point(size = .3) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.2, "yes" = 1)) +
  labs(
    x = "Log Fold Change",
    y = "-log10(P-Value)",
    title = "UC Andersen"
  ) +
  theme_minimal()

## Multiple sclerosis ----------------------------------------------------------
MS <- MS %>%
  arrange(desc(stat)) %>%
  mutate(idx = 1:length(stat)) %>%
  mutate(top_2000 = if_else(ENSG.ID %in% filter(all.genes, Position <= 2000)$ENSG.ID, "yes", "no"),
         top_100 = if_else(ENSG.ID %in% filter(all.genes, Position <= 100)$ENSG.ID, "yes", "no"))

ggplot(MS, aes(x = log2FoldChange, y = -log10(pvalue), color = top_2000, alpha = top_2000)) +
  #ggplot(MS, aes(x = log2FoldChange, y = -log10(padj), color = top_2000, alpha = top_2000)) +
  geom_point(size = .3) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.3, "yes" = 1)) +
  labs(
    x = "Log Fold Change",
    y = "-log10(P-Value)",
    title = "Multiple sclerosis"
  ) +
  theme_minimal()

## UC Hansen -------------------------------------------------------------------
UC.hansen<- UC.hansen %>%
  arrange(desc(stat)) %>%
  mutate(idx = 1:length(stat)) %>%
  mutate(top_2000 = if_else(ENSG.ID %in% filter(all.genes, Position <= 2000)$ENSG.ID, "yes", "no"),
         top_100 = if_else(ENSG.ID %in% filter(all.genes, Position <= 100)$ENSG.ID, "yes", "no"))

# no adj pval
# The “Significant” column contains a “+” if a protein met the selected 
# significance threshold (usually q-value). Additionally, p-values 
# (probability of type I error) and the corresponding q-values (corrected p-value) 
# are provided in the output table. https://link.springer.com/protocol/10.1007/978-1-4939-7493-1_7#Sec16 
# from perseus protocol, I assume p values here are not corrected
UC.hansen$adj_pvalue <- p.adjust(UC.hansen$pvalue, method = "BH")

UC.hansen$significance <- ifelse(
  UC.hansen$adj_pvalue < 0.01, 
  ifelse(UC.hansen$top_2000 == "yes", "significant_inflammatory", "significant"), 
  "not_significant"
)

write_tsv(UC.hansen, "data/05_UC_hansen_processed.tsv")

ggplot(UC.hansen, aes(x = logFC, y = -log10(adj_pvalue), color = significance, alpha = significance)) +
  geom_point(size = .3) +
  scale_color_manual(values = c(
    "significant" = "black", 
    "significant_inflammatory" = "red", 
    "not_significant" = "grey"
  )) +
  scale_alpha_manual(values = c("significant" = .2, 
                                "significant_inflammatory" = 1,
                                "not_significant" = .2)) +
  labs(
    x =  expression(log[2]("fold change")),
    y = expression(-log[10]("adj. p-value")),
    title = "Ulcerative colitis vs. healthy controls (Schniers et al.)"
  ) +
  geom_hline(yintercept = 2, linewidth = .3, linetype = "dashed") +
  coord_cartesian(xlim=c(-2, 2)) +
  plotTheme +
  theme(legend.position = "none")
# weird that some pvals are == 1

ggsave("figures/05_volcano_UCprot_usecase.png", h = 8, w = 10)

## Plot all --------------------------------------------------------------------
igan <- igan %>%
  mutate(case = "igan")

UC.andersen <- UC.andersen %>%
  mutate(case = "UC.andersen")

AH <- AH %>%
  mutate(case = "AH")

MS <- MS %>%
  mutate(case = "MS")

UC.hansen <- UC.hansen %>%
  mutate(case = "UC.hansen")

test <- rbind(dplyr::select(igan, ENSG.ID, stat, top_2000, top_100, case),
              dplyr::select(UC.andersen, ENSG.ID, stat,top_2000, top_100, case),
              dplyr::select(AH, ENSG.ID, stat, top_2000, top_100, case),
              dplyr::select(MS, ENSG.ID, stat, top_2000, top_100, case),
              dplyr::select(UC.hansen, ENSG.ID, stat, top_2000, top_100, case))

test %>%
  filter(case != "UC.hansen") %>%
  ggplot(aes(x=case, y = stat, color = top_2000, alpha = top_2000, size = top_2000)) +
  geom_jitter() +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.3, "yes" = 1)) +
  scale_size_manual(values = c("no" = 0.1, "yes" = .3)) +
  labs(
    x = "Dataset",
    y = "stat"
  ) +
  theme_minimal() 

test %>%
  filter(case != "UC.hansen") %>%
  ggplot(aes(x=case, y = stat, color = top_2000, alpha = top_2000, size = top_2000)) +
  geom_jitter() +
  geom_hline(yintercept = 0, linewidth = .3) +
  scale_color_manual(values = c("no" = "black", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.1, "yes" = 1)) +
  scale_size_manual(values = c("no" = 0.1, "yes" = .3)) +
  labs(
    x = "Dataset",
    y = "stat"
  ) +
  theme_minimal() 

test %>%
  filter(case == "UC.hansen") %>%
  ggplot(aes(x=case, y = stat, color = top_2000, alpha = top_2000, size = top_2000)) +
  geom_jitter() +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("no" = "black", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.2, "yes" = 1)) +
  scale_size_manual(values = c("no" = 0.1, "yes" = .3)) +
  labs(
    x = "Dataset",
    y = expression(logFC %*% -log(pval))
  ) +
  theme_minimal() 

for(case2 in distinct(test, case)$case){
  print(case2)
  p <- test %>%
    filter(case == case2) %>%
    ggplot(aes(x=case, y = stat, color = top_2000, alpha = top_2000, size = top_2000)) +
    geom_jitter() +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = c("no" = "black", "yes" = "red")) +
    scale_alpha_manual(values = c("no" = 0.1, "yes" = 1)) +
    scale_size_manual(values = c("no" = 0.1, "yes" = .4)) +
    labs(
      x = "Dataset",
      y = "stat"
    ) +
    theme_minimal() 
  
  print(p)
}



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

# Prepare gene sets ------------------------------------------------------------
# top 100, msigdb hallmark inflammatory response, GO:BP inflammatory response
### UC Andersen old results format ------------------
#UC.andersen <- read_tsv("data/raw/DE.res.UC.Andersen.tsv")  # DE by Oana, t stat
#UC.andersen <- UC.andersen %>%
#  rename(ENSG.ID = ensembl_gene_id,
#         stat = t) %>%
#  filter(ENSG.ID %in% Final_Annotation_List$ENSG.ID)



