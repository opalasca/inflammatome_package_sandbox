# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load packages ----------------------------------------------------------------
library(tidyverse)
library(biomaRt)
library(httr)
library(jsonlite)

# Gene annotation list  ----------------------------------------------------

#Created 02/06/2023

#Annotation database: BioMart
mart.human <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                               dataset = "hsapiens_gene_ensembl",
                               host = "https://www.ensembl.org")

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",version="112")

#look for individual attributes
#Attributes <- mart.human@attributes

genes.only.human.symbol <- getBM(
  attributes = c("ensembl_gene_id" ,"gene_biotype", "external_gene_name", "chromosome_name"),
  mart = ensembl)

colnames(genes.only.human.symbol)=c("ENSG.ID", "Gene.type", "Gene.name", "Chromosome")

#genes.protein.human.symbol <- getBM(
#attributes = c("ensembl_gene_id" ,"gene_biotype", "external_gene_name", "ensembl_peptide_id", "uniprotswissprot", "chromosome_name"),
#mart = mart.human)

#colnames(genes.protein.human.symbol)=c("ENSG.ID", "Gene.type", "Gene.name", "ENSP.ID", "UniProt.ID", "Chromosome")


#Creating final annotation list by...
## - converting missing Gene.names to NA (ENSG.ID available)
## - filtering for main ENSEMBL entries of a gene (by chromosome positions)
## - keeping NA values separate until adding them in the end
## - removing pseudogenes for all Gene.name duplicates (= Gene.name groups with more than one member)
## - choosing lowest ENSG.ID number entry for genes with multiple chromosomal locations

Chrom <- str_extract(genes.only.human.symbol$Chromosome, "^.{1,5}$") %>% 
  unique()
ChromOI <- Chrom[!is.na(Chrom)]
ChromOI

Final_Annotation_List <- genes.only.human.symbol %>% 
  mutate(Gene.name = na_if(Gene.name, "")) %>% 
  filter(Chromosome %in% ChromOI)

Annotation_NAs <- Final_Annotation_List %>% 
  filter(is.na(Gene.name))

Final_Annotation_List <- Final_Annotation_List %>%
  filter(!is.na(Gene.name)) %>%
  group_by(Gene.name) %>%
  filter(case_when(
    n() > 1 ~ !str_detect(Gene.type, "pseudogene"),
    n() == 1 ~ str_detect(Gene.type, ".")
  )) %>%
  ungroup() %>%
  mutate(ENSG.number = str_extract(ENSG.ID, "\\d+")) %>%
  group_by(Gene.name) %>%
  filter(ENSG.number == min(ENSG.number)) %>%
  ungroup() %>% 
  dplyr::select(-ENSG.number) %>%
  base::rbind(Annotation_NAs)

# Load annotation list from 02/06/2023
Final_Annotation_List <- read_tsv(paste0("data/", "Final_Annotation_List_BioMart.tsv"))