
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


########## Extract biomart  ##########################
######################################################

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
a=ensembl@attributes
genes.uniprot <- getBM(attributes = c("uniprotswissprot", "ensembl_gene_id", "hgnc_symbol","gene_biotype"), mart = ensembl)
genes.refseq <- getBM(attributes = c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol","gene_biotype"), mart = ensembl)
genes.refseq.ncrna <- getBM(attributes = c("refseq_ncrna", "ensembl_gene_id", "hgnc_symbol","gene_biotype"), mart = ensembl)
genes.entrez <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "hgnc_symbol","gene_biotype"), mart = ensembl)

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



