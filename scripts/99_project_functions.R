library(tidyverse)

# Load annotations and top 100 markers ------------------------------------------

Final_Annotation_List <- read_tsv("data/Final_Annotation_List_Biomart.tsv",show_col_types = FALSE) %>%
  filter(Gene.type == "protein_coding")

## text-mined markers ----------------------------------------------------------
Inflamm_Top100_marker <- read_tsv("data/Inflammation_markers_Top100.tsv",show_col_types = FALSE)

length(intersect(Inflamm_Top100_marker$ENSG.ID, Final_Annotation_List$ENSG.ID))

STRING_ENSP_to_ENSG <- read_tsv( "data/human.aliases.filtered.txt",
                                 col_names = c("SPLIT", "ENSG.ID", "data.type"), show_col_types = FALSE) %>% 
  separate_wider_delim(SPLIT, ".", names = c("tax", "ENSP.ID")) %>%
  dplyr::select(-c(tax, data.type)) 

## msigdb inflammation markers -------------------------------------------------
library(msigdbr)

h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, ensembl_gene, gene_symbol)

infl.markers <- h[h$gs_name == "HALLMARK_INFLAMMATORY_RESPONSE",]
msigdb.markers = infl.markers[c(2,3)]
colnames(msigdb.markers)=c("ENSG.ID","Gene.name")
# check if any are non-coding
msigdb.markers %>% 
  left_join(Final_Annotation_List) %>%
  dplyr::count(Gene.type)
# there are 22 with NA in gene.type
msigdb.markers %>% 
  left_join(Final_Annotation_List) %>%
  filter(is.na(Gene.type)) %>% print(n=22)
# if you look at gene.name they are all in annotation as protein coding,
# but the ENSG.IDs don't match

# keep only markers present in our list of protein coding genes
msigdb.markers <- msigdb.markers %>% 
  filter(ENSG.ID %in% Final_Annotation_List$ENSG.ID)

## Text-mined and msigdb union -------------------------------------------------
union <- full_join(msigdb.markers, Inflamm_Top100_marker)

## GO terms --------------------------------------------------------------------
go.markers <- read_tsv("data/markers.go.tsv", show_col_types = FALSE)
go.markers.pos <- go.markers %>%
  filter(gene_biotype =="protein_coding", 
         GO.name != "negative regulation of inflammatory response", 
         ENSG.ID %in% Final_Annotation_List$ENSG.ID)

go.markers.neg <- go.markers %>%
  filter(gene_biotype =="protein_coding", 
         GO.name == "negative regulation of inflammatory response", 
         ENSG.ID %in% Final_Annotation_List$ENSG.ID)

go.unique <- data.frame(ENSG.ID = unique(go.markers$ENSG.ID)) %>%
  filter(ENSG.ID %in% Final_Annotation_List$ENSG.ID)

# Functions --------------------------------------------------------------------

## Expression data preprocessing -----------------------------------------------
geo2r_preprocess <- function(expset){
  # make proper column names to match toptable 
  fvarLabels(expset) <- make.names(fvarLabels(expset)) 
  
  # log2 transformation
  ex <- exprs(expset)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { ex[which(ex <= 0)] <- NaN 
  #exprs(expset) <- log2(ex) 
  ex <- log2(ex)
  }
  
  # check if transformation was applied
  if (identical(ex, exprs(expset))){print("expression data not transformed")
  }else{print("log2 transformation applied")}
  
  return(ex)
}

## merge annotation when no annotationdb avail ---------------------------------
ps_to_ensg_fromfile <- function(expdata, mapping){
  # cases where one ensg matches multiple probes take the one with hightest ave_exp
  # cases when one probe matches multiple ensgs remove
  # first keep unique probe set - ensg pairs
  mapping <- mapping %>%
    distinct(entrezg, ensg) %>%
    drop_na()
  
  mapping <- mapping %>%
    right_join(rownames_to_column(as.data.frame(expdata), var="entrezg")) %>%
    drop_na()
  
  # remove probes matching multiples ensg ids
  dup <- mapping %>%
    dplyr::count(entrezg) %>%
    filter(n>1)
  
  mapping <- mapping %>% filter(!(entrezg %in% dup$entrezg)) 
  
  # choose highest average expression probe set 
  mapping <- mapping %>% 
    mutate(ave_exp = apply(mapping[,pd$geo_accession], 1, mean)) %>%
    group_by(ensg) %>%
    arrange(desc(ave_exp), .by_group=TRUE) %>% 
    filter(row_number()==1) %>%
    ungroup() %>%
    dplyr::select(-c(entrezg, ave_exp)) %>%
    mutate(ensg = gsub("_at", "", ensg)) %>%
    column_to_rownames(var = "ensg")
  
  return(mapping)
}

## Individual benchmark plot ---------------------------------------------------
benchmark_plot <- function(data, dataset_id){
  gold_set <- union
  sorted <- data[order(-data$t),] # sort by decreasing lfc
  n <- dim(sorted)[1] # number of rows
  sorted$index = c(1:n)
  sorted$enr = 0; prev = 0
  rank.pval.thr = 0; rank.pval.thr.y = 0
  thr = 0.05
  for (i in 1:n){
    if (sorted[i,"ENSG.ID"] %in% gold_set$ENSG.ID){
      sorted[i,"enr"] = prev+1; prev = prev+1 } else {sorted[i,"enr"] = prev} # counting gold set occurences 
    if (!is.na(sorted[i,"adj.P.Val"]) & sorted[i,"adj.P.Val"]>thr & rank.pval.thr==0) {
      rank.pval.thr=i
      rank.pval.thr.y=sorted$enr[i]
    }}
  segment.benchmark <- data.frame(x1 = 0,
                                  y1 = 0,
                                  x2 = dim(filter(Final_Annotation_List, Gene.type == "protein_coding"))[1],
                                  y2 = dim(union)[1])
  #Create benchmarking plots 
  ggplot(sorted, aes(x = index, y = enr)) +
    geom_point(alpha = 0.3) + 
    geom_point(aes(x = rank.pval.thr, y = rank.pval.thr.y, size = 2,colour="red"))+
    geom_segment(aes(x = x1,
                     y = y1,
                     xend = x2,
                     yend = y2), data = segment.benchmark, colour = "red", linetype="dotted")+
    geom_vline(aes(xintercept = rank.pval.thr), linetype = "dotted")+
    labs(title = paste(dataset_id, " benchmark"),
          y = "Occurence of inflammation markers",
         x = "DEGs sorted by stat value") +
    theme_classic() +
    guides(size = "none", linetype = "none", colour = "none")
}

## Individual benchmark plot for RNAseq ----------------------------------------
benchmark_plot_RNA <- function(data, dataset_id){
  gold_set <- union
  sorted <- data[order(-data$stat),] # sort by decreasing stat
  n <- dim(sorted)[1] # number of rows
  sorted$index = c(1:n)
  sorted$enr = 0; prev = 0
  rank.pval.thr = 0; rank.pval.thr.y = 0
  thr = 0.05
  segment.benchmark <- data.frame(x1 = 0,
                                  y1 = 0,
                                  x2 = dim(filter(Final_Annotation_List, Gene.type == "protein_coding"))[1],
                                  y2 = dim(union)[1])
  for (i in 1:n){
    if (sorted[i,"ENSG.ID"] %in% gold_set$ENSG.ID){
      sorted[i,"enr"] = prev+1; prev = prev+1 } else {sorted[i,"enr"] = prev} # counting gold set occurences 
    if (!is.na(sorted[i,"padj"]) & sorted[i,"padj"]>thr & rank.pval.thr==0) {
      rank.pval.thr=i
      rank.pval.thr.y=sorted$enr[i]
    }}
  #Create benchmarking plots 
  ggplot(sorted, aes(x = index, y = enr)) +
    geom_point(alpha = 0.3) + 
    geom_point(aes(x = rank.pval.thr, y = rank.pval.thr.y, size = 2,colour="red"))+
    geom_segment(aes(x = x1,
                    y = y1,
                    xend = x2,
                    yend = y2), data = segment.benchmark, colour = "red", linetype="dotted")+
    geom_vline(aes(xintercept = rank.pval.thr), linetype = "dotted")+
    labs(title = paste(dataset_id, " benchmark"),
         y = "Occurence of inflammation markers",
         x = "DEGs sorted by stat value") +
    theme_classic() +
    guides(size = "none", linetype = "none", colour = "none")
}

## probe sets to ENSG ids for HG-U133_Plus_2 annotation ------------------------
#
#ps_to_ensg <- function(toptable){
#
#  res <- toptable %>%
#    rownames_to_column(var = "ps") %>%
#    left_join(distinct(group_by(HGU133Plus2, ENSG.ID, ps))) %>% 
#    drop_na() %>%
#    group_by(ENSG.ID) %>%
#    arrange(desc(abs(logFC)), .by_group=TRUE) %>% # cases where one ensg matches multiple probes take the one with hightest LFC
#    filter(row_number()==1) %>%
#    ungroup() %>%
#    group_by(ps) %>%
#    mutate(ps_count = n()) %>%
#    ungroup() %>%
#    filter(ps_count == 1) %>%
#    mutate(ENSG.ID = gsub("_at", "", ENSG.ID))
#  
##    # multiple ensg match same probe --> discard
#    # res <- res[!duplicated(res$ps), ] %>%
#    # mutate(ENSG.ID = gsub("_at", "", ENSG.ID))
#    
#    return(res)
#}
#
## Calculate correlation of gene exp levels with inflammation score ------------

get_corr <- function(pdata, exdata, score){
  # subset and transpose expression data
  t_ex <- t(exdata[,pdata$geo_accession])
  # calculate correlation
  res <- as.data.frame(cor(t_ex, score)) %>% 
    dplyr::rename(corr = V1) %>% 
    drop_na() %>%
    rownames_to_column(var = "ps")
  # for chosing ps that match same ensg.id
  res$av_exp <- apply(t_ex[,res$ps], 2, mean)
  
  return(res)
}

## benchmark plot for correlation approach -------------------------------------

benchmark_corr <- function(df, gold_set, set_str, group_str){
  # initialize 
  count <- 0 
  df <- df %>% 
    arrange(desc(corr)) %>%
    add_column(marker_count = NA,
               idx = 1:length(df$ENSG.ID)) 
  
  for (i in 1:length(df$ENSG.ID)){
    if (df$ENSG.ID[i] %in% gold_set$ENSG.ID){
      count <- count+1
    }
    df$marker_count[i] <- count
  }
  
  #Create benchmarking plots 
  ggplot(df, aes(x = idx, y = marker_count)) +
    geom_step() + 
    geom_segment(aes(x = 0, y = 0,
                     xend = max(idx),
                     yend = max(marker_count),
                     colour = "red", linetype="dotted"))+
    labs(title = paste(set_str, group_str, "benchmark", sep = " "),
         y = "Occurence of inflammation markers",
         x = "Genes sorted by correlation to severity score") +
    theme_classic() +
    guides(size = "none", linetype = "none", colour = "none")
}

## Function to count marker occurences in a data set ---------------------------
# updated to include blood rna seq datasets
count_markers <- function(filepath, gold_set, coding_only, intersect_only, intersect_ids){
  print(filepath)
  df <- read_tsv(filepath) %>%
    mutate(Chromosome = as.character(Chromosome),
           Contrast = gsub("\\w+\\/([\\w.]+)\\.tsv",
                          "\\1",
                          filepath, perl = T))
  if(coding_only == T){
    df <- df %>% filter(Gene.type == "protein_coding")
  }
  if(intersect_only == T){
    df <- df %>% 
      filter(ENSG.ID %in% intersect_ids)
  }
  try({
    df <- df %>% dplyr::rename(t = stat, 
                 adj.P.Val = padj)
  })
 
  df <- df %>% 
    arrange(desc(t)) %>%
    add_column(marker_count = NA) 
  df$idx <- 1:length(df$ENSG.ID)
  
  # ROC-like version
  # initialize 
  pos_count <- 0 
  neg_count <- 0
  marker_count <-c()
  idx_ROC <- c()
  
  for (i in 1:length(df$ENSG.ID)){
    if (df$ENSG.ID[i] %in% gold_set$ENSG.ID){
      pos_count <- pos_count+1
    } else {
      neg_count <- neg_count+1
    }
    df$marker_count[i] <- pos_count
    df$idx_ROC[i] <- neg_count
  }
  
  df <- df %>%
    mutate(non_sig_y = ifelse((idx == max(filter(df, adj.P.Val > 0.05)$idx) |
                                 idx == min(filter(df, adj.P.Val > 0.05)$idx)),
                              marker_count,
                              NA
    ))
}

## BIRRA function --------------------------------------------------------------
#BIRRA aggregation method (performs well for large set number and heterogeneity)
##https://sites.pitt.edu/~mchikina/BIRRA/runBIRRA.R

BIRRA=function(data, prior=0.05, num.bin=50, num.iter=15, return.all=F, plot.legend=F, grp=NULL, cor.stop=1, ...){
  nr=nrow(data)
  nrp=floor(nrow(data)*prior)
  #data=apply(-data,2,rank)/nr #Replaced by the following line of code 
  data = data #Input file already consists of normalized ranks per dataset
  nc=ncol(data)
  TPR=FPR=Bayes.factors=matrix(ncol=nc, nrow=num.bin)
  binned.data=ceiling(data*num.bin)
  bayes.data=matrix(nrow=nrow(data), ncol=ncol(data))
  guess=apply(data,1,mean)
  cprev=0
  #par(mfrow=c(floor(sqrt(num.iter)), ceiling(sqrt(num.iter))), mai=rep(0.7,4))
  for ( iter in 1:num.iter){
    if((cor.stop-cprev)>1e-15){
      guesslast=guess
      oo=order(guess)
      guess[oo[1:nrp]]=1
      guess[oo[(nrp+1):nr]]=0
      for (i in 1:nc){  
        for (bin in 1:num.bin){
          frac=bin/num.bin
          TPR=sum(guess[binned.data[,i]<=bin])
          FPR=sum((!guess)[binned.data[,i]<=bin])
          Bayes.factors[bin,i]=log((TPR+1)/(FPR+1)/(prior/(1-prior)))
        }
      }
      Bayes.factors=apply(Bayes.factors,2,smooth)
      Bayes.factors=apply(Bayes.factors,2,function(x){rev(cummax(rev(x)))})
      # Plot TPR vs bin for each data set
      if(is.null(grp)){
        matplot(1:num.bin, Bayes.factors, type="l", lwd=2, ...)
      }
      else{
        matplot(1:num.bin, Bayes.factors, type="l", lwd=2, lty=grp, col=grp)
      }
      title(paste("Iteration", iter))
      if (iter==1&plot.legend){
        legend("topright", col=1:5, lty=1:4, legend=colnames(data), lwd=2, ncol=2)
      }
      for (bin in 1:num.bin){
        oo=order(Bayes.factors[bin,], decreasing=T)
        Bayes.factors[bin, oo[1]]=Bayes.factors[bin, oo[2]]
      }
      for (i in 1:nc){
        bayes.data[,i]=Bayes.factors[binned.data[,i],i]
      }
      bb=exp(apply(bayes.data,1, sum))
      f=prior/(1-prior)
      prob=bb*f/(1+bb*f)
      exp=sort(prob, decreasing=F)[nrp]
      guess=rank(-apply(bayes.data,1, sum))
      cprev=cor(guess, guesslast)
      message("correlation with pervious iteration=",cprev)
    }
    else{
      message("Converged");
      break
    }
  }
  if(return.all){
    return(list(result=guess, data=bayes.data, BF=Bayes.factors))
  }
  else{
    guess
  }
}

BIRRA_return_all=function(data, prior=0.05, num.bin=50, num.iter=15, return.all=T, plot.legend=F, grp=NULL, cor.stop=1, ...){
  nr=nrow(data)
  nrp=floor(nrow(data)*prior)
  #data=apply(-data,2,rank)/nr #Replaced by the following line of code 
  data = data #Input file already consists of normalized ranks per dataset
  nc=ncol(data)
  TPR=FPR=Bayes.factors=matrix(ncol=nc, nrow=num.bin)
  binned.data=ceiling(data*num.bin)
  bayes.data=matrix(nrow=nrow(data), ncol=ncol(data))
  guess=apply(data,1,mean)
  cprev=0
  #par(mfrow=c(floor(sqrt(num.iter)), ceiling(sqrt(num.iter))), mai=rep(0.7,4))
  for ( iter in 1:num.iter){
    if((cor.stop-cprev)>1e-15){
      guesslast=guess
      oo=order(guess)
      guess[oo[1:nrp]]=1
      guess[oo[(nrp+1):nr]]=0
      for (i in 1:nc){  
        for (bin in 1:num.bin){
          frac=bin/num.bin
          TPR=sum(guess[binned.data[,i]<=bin])
          FPR=sum((!guess)[binned.data[,i]<=bin])
          Bayes.factors[bin,i]=log((TPR+1)/(FPR+1)/(prior/(1-prior)))
        }
      }
      Bayes.factors=apply(Bayes.factors,2,smooth)
      Bayes.factors=apply(Bayes.factors,2,function(x){rev(cummax(rev(x)))})
      # Plot TPR vs bin for each data set
      if(is.null(grp)){
        matplot(1:num.bin, Bayes.factors, type="l", lwd=2, ...)
      }
      else{
        matplot(1:num.bin, Bayes.factors, type="l", lwd=2, lty=grp, col=grp)
      }
      title(paste("Iteration", iter))
      if (iter==1&plot.legend){
        legend("topright", col=1:5, lty=1:4, legend=colnames(data), lwd=2, ncol=2)
      }
      for (bin in 1:num.bin){
        oo=order(Bayes.factors[bin,], decreasing=T)
        Bayes.factors[bin, oo[1]]=Bayes.factors[bin, oo[2]]
      }
      for (i in 1:nc){
        bayes.data[,i]=Bayes.factors[binned.data[,i],i]
      }
      bb=exp(apply(bayes.data,1, sum))
      f=prior/(1-prior)
      prob=bb*f/(1+bb*f)
      exp=sort(prob, decreasing=F)[nrp]
      guess=rank(-apply(bayes.data,1, sum))
      cprev=cor(guess, guesslast)
      message("correlation with pervious iteration=",cprev)
    }
    else{
      message("Converged");
      break
    }
  }
  if(return.all){
    return(list(result=guess, data=bayes.data, BF=Bayes.factors))
  }
  else{
    guess
  }
}

## Contrast selection ----------------------------------------------------------
# function to calculate auc for each contrast
calculate_auc <- function(contrast){
  curve <- approxfun(contrast$idx, contrast$marker_count, ties = "ordered")
  auc <- integrate(curve, min(contrast$idx), max(contrast$idx), subdivisions = 10000)
  
  ROC_curve <- approxfun(contrast$idx_ROC, contrast$marker_count, ties = "ordered")
  ROC_auc <- integrate(ROC_curve, min(contrast$idx_ROC), max(contrast$idx_ROC), subdivisions = 10000)
  
  return(data.frame(AUC = auc$value, AUC_ROC = ROC_auc$value))
}

calculate_auc_p <- function(contrast){
  curve <- approxfun(contrast$p_annotation, contrast$p_markers, ties = "ordered")
  auc <- integrate(curve, min(contrast$p_annotation), max(contrast$p_annotation), subdivisions = 10000)
  ROC_curve <- approxfun(contrast$p_ROC, contrast$p_markers, ties = "ordered")
  ROC_auc <- integrate(ROC_curve, min(contrast$p_ROC), max(contrast$p_ROC), subdivisions = 10000)
  return(data.frame(AUC_p = auc$value, AUC_ROC_p = ROC_auc$value))
}

# function to select contrast
# argument is auc cutoff based on best auc --> change it to be above expected by chance
# when more than one contrast from same dataset and disease the one with highest AUC is chosen
select_contrast <- function(df, auc.cutoff){
  df_select <- df %>%
    group_by(Contrast) %>%
    nest() %>%
    mutate(AUC = map(data, calculate_auc)) %>%
    unnest(cols = c(data, AUC)) %>%
    group_by(Contrast) %>%
    nest() %>%
    mutate(AUC_p = map(data, calculate_auc_p)) %>%
    unnest(cols = c(data, AUC_p)) %>%
    group_by(Contrast) %>%
    drop_na(non_sig_y) %>%
    mutate(sig_thres = min(non_sig_y)) %>%
    dplyr::select(Contrast, sig_thres, AUC, AUC_ROC, AUC_p, AUC_ROC_p) %>%
    distinct() %>%
    ungroup() %>%
    mutate(include = if_else(AUC_ROC_p >= auc.cutoff, 
                             "yes",
                             "no")) %>% 
    separate_wider_delim(Contrast, "_", names = c(NA, "dataset", "case", "control"),
                         too_many = "merge",
                         cols_remove = F) %>%
    mutate(include = if_else(AUC_p == max(AUC_p) & include == "yes", "yes", "no"),
           .by = c("dataset", "case")) %>%
    mutate(tissue = case_when(case %in% c("PSO", "AD", "Acontact", "mixed", "Contact") ~ "skin",
                              case %in% c("MS", "Alz") ~ "CNS",
                              case %in% c("OA", "RA") ~ "synovium",
                              case %in% c("CD", "UC") ~ "intestine",
                              case %in% c("SD", "COPD", "TB", "AC", "CHP", "IPF", "SSc") ~ "lung",
                              case %in% c("IgAN", "DiabNeph", "FSGS", "LupusNeph", "MCD", "MN") ~ "kidney",
                              case %in% c("HCV", "TB", "NAFLD", "NASH", "AH", "Cirr") | (case == "OB" & dataset == "GSE126848") ~ "liver",
                              (case == "OB" & dataset == "GSE244120") ~ "muscle",
                              (case == "OB" & dataset == "GSE244118") ~ "fat",
                              case == "SS" ~ "salivary_glands"))
  return(df_select)
}

# function to calculate auc for each rank agg list (used for choosing AUC threshold for contrast inclusion)
calculate_auc_rank <- function(rank){
  curve <- approxfun(rank$Position, rank$enr, ties = "ordered")
  auc <- integrate(curve, min(rank$Position), max(rank$Position), subdivisions = 1000)
  ROC_curve <- approxfun(rank$idx_ROC, rank$enr, ties = "ordered")
  ROC_auc <- integrate(ROC_curve, min(rank$idx_ROC), max(rank$idx_ROC), subdivisions = 1000)
  return(data.frame(AUC = auc$value, AUC_ROC = ROC_auc$value))
}

# function to create a color palette asigning grey to excluded contrasts
assign_color <- function(x, df, pal) {
  ifelse(x %in% df, pal[match(x, df)],
         "gray")
}

## Produce rank aggregated list from contrast list -----------------------------

rank_aggregation <- function(contrast_list, gene_type){
  set_pool <- list()
  
  for(d in contrast_list){
    print(d)
    set_df <- read_tsv(paste("data/", d, ".tsv", sep = "")) %>% 
      dplyr::filter(Gene.type == gene_type)
    try({set <- set_df %>% 
      arrange(desc(stat))})
    try({set <- set_df %>% 
      arrange(desc(t))})
    set_pool <- c(set_pool, list(set$ENSG.ID))
  }
  
  ## Assign normalized ranks (between 0 and 1); all NAs receive rank "1" 
  Rank_sets <- RobustRankAggreg::rankMatrix(set_pool) 
  
  Rank_results <- BIRRA(Rank_sets) # increased number of iterations
  entities <- rownames(Rank_sets)
  rankedEntities <- entities[order(Rank_results)]
  
  head(entities, 10)
  head(rankedEntities, 10)
  
  Aggreg_list <- as.data.frame(rankedEntities)
  colnames(Aggreg_list) <- "ENSG.ID"
  Aggreg_list <- rownames_to_column(Aggreg_list, var = "Position")
  
  Aggreg_list_final_pos <- left_join(Aggreg_list, Final_Annotation_List, by = "ENSG.ID")
  Aggreg_list_final_pos <- left_join(Aggreg_list_final_pos, STRING_ENSP_to_ENSG, by = "ENSG.ID")
  Aggreg_list_final_pos$Position <- as.numeric(Aggreg_list_final_pos$Position)
  
  #artificially add "stat" column such as to put together with the other plots
  Aggreg_list_final_pos$stat=max(Aggreg_list_final_pos$Position)-Aggreg_list_final_pos$Position +1
  
  return(Aggreg_list_final_pos)
  
}

## Count marker occurences in rank aggregated list -----------------------------
count_markers_rank <- function(list, gold_set){
  n <- dim(list)[1]
  # list$enr = 0; prev = 0
  # for (i in 1:n){
  #   if (list[i,"ENSG.ID"] %in% gold_set$ENSG.ID){
  #     list[i,"enr"] = prev+1; prev = prev+1 } else {list[i,"enr"] = prev}
  # }
  # return(list)
  
  pos_count <- 0 
  neg_count <- 0
  enr <-c()
  idx_ROC <- c()
  
  for (i in 1:n){
    if (list[i,"ENSG.ID"] %in% gold_set$ENSG.ID){
      pos_count <- pos_count+1
    } else {
      neg_count <- neg_count+1
    }
    list$enr[i] <- pos_count
    list$idx_ROC[i] <- neg_count
  }
  
  return(list)
}

## Dataset QC functions --------------------------------------------------------
### PCA ------------------------------------------------------------------------
# performs vst, saves pca plot to figures/ and returns sample scores(pc1 pc2 coords)
# takes deseq2 object, grouping variable and geo accession as input
PCA_norm <- function(dds, condition, filename){
  
  ## Count transformation
  vsd <- vst(dds, blind = FALSE)
  
  ## Plot
  PCA <- plotPCA(vsd, intgroup = condition, ntop = nrow(dds))
  ggsave(paste("figures/", filename, ".png", sep=""), PCA)
  
  ## Sample info/coordinates
  return(plotPCA(vsd, intgroup = condition, ntop = nrow(dds), returnData = TRUE))
}

### Boxplots -------------------------------------------------------------------
# Boxplot of input counts (before DESeq2, cts: raw counts)

Boxplot_prior <- function(cts, filename){
  
  set.seed(1234)
  png(paste("figures/", filename, ".png", sep = ""))
  if(dim(cts)[2] > 10){
    boxplot(log2(cts[, sample(1:dim(cts)[2], 10)]+1), 
                               ylab = "RNAseq counts [log2(n+1)]",
                               xaxt = "n",
                               main = "Count distribution of input data (10 random samples)")
  } else {
    boxplot(log2(cts+1), 
            ylab = "RNAseq counts [log2(n+1)]",
            xaxt = "n",
            main = "Count distribution of input data (10 random samples)")
    }
  
  dev.off()
}

# Boxplot of normalized counts (after DESeq2, cts.norm: normalized for library size, NOT transformed [log2, vst, etc])

Boxplot_post <- function(cts.norm, filename){
  
  set.seed(1234)
  png(paste("figures/", filename, ".png", sep = ""))
  if(dim(cts.norm)[2] > 10){boxplot(log2(cts.norm[, sample(1:dim(cts.norm)[2], 10)]+1), 
                               ylab = "RNAseq counts [log2(n+1)]",
                               xaxt = "n",
                               main = "Count distribution of DESeq2 output (10 random samples)")
  }else{
    boxplot(log2(cts.norm+1), 
                ylab = "RNAseq counts [log2(n+1)]",
                xaxt = "n",
                main = "Count distribution of DESeq2 output (10 random samples)")}
  
  dev.off()
}

## Correlation of topx genes to severity scores -------------------------------
score_corr_rna <- function(x, score, inf, metadata, vst, vst.center, random, list){
  if(random == "yes"){
    topx <- row.names(vst)[runif(x, min=1, max=length(row.names(vst)))]
  } else {
    topx <- intersect(row.names(vst), inf$ENSG.ID[1:x])
  }
  if(length(topx) == 1){
    med.corr <- cor(meta[score], vst[topx,])
    med.center.corr <- cor(meta[score], vst.center[topx,])
    mean.corr <- med.corr
    mean.center.corr <- med.center.corr
  } else {
    med.corr <- cor(meta[score], colMedians(vst[topx,], useNames =T))
    med.center.corr <- cor(meta[score], colMedians(vst.center[topx,], useNames =T))
    mean.corr <- cor(meta[score], colMeans(vst[topx,]))
    mean.center.corr <- cor(meta[score], colMeans(vst.center[topx,]))
  }
  res <- data.frame("corr.med"= med.corr, "corr.med.center" = med.center.corr,
                    "corr.mean" = mean.corr, "corr.mean.center" = mean.center.corr,
                    "top_x"=x, "intersect" = length(topx), inf=list)
  return(res)
}

## Lmer coefficient topx ~ severity -------------------------------------------
# In dataset with paired samples, patient is modeled as random effects
score_lmer_rna <- function(x, score, inf, metadata, vst, vst.center, random, list){
  if(random == "yes"){
    topx <- row.names(vst)[runif(x, min=1, max=length(row.names(vst)))]
  } else {
    topx <- intersect(row.names(vst), inf$ENSG.ID[1:x])
  }
  if(length(topx) == 1){
    df <- data.frame(center = vst.center[topx,],
                     no.center = vst[topx,],
                     tlss = metadata[score],
                     patient = as.factor(metadata$patient))
    mod.center <- lmer(center ~ total_tlss + (1|patient), data = df)
    mod.no.center <- lmer(no.center ~ total_tlss + (1|patient), data = df)
    res <- data.frame("beta.med"= mod.no.center@beta[2], "beta.med.center" = mod.center@beta[2],
                      "beta.mean" = mod.no.center@beta[2], "beta.mean.center" = mod.center@beta[2],
                      "top_x"=x, "intersect" = length(topx), inf=list) 
  } else {
    df <- data.frame(med.center = colMedians(vst.center[topx,], useNames = T),
                     med = colMedians(vst[topx,], useNames = T),
                     mean.center = colMeans(vst.center[topx,]),
                     mean = colMeans(vst[topx,]),
                     tlss = metadata[score],
                     patient = as.factor(metadata$patient))
    mod.med.center <- lmer(med.center ~ total_tlss + (1|patient), data = df)
    mod.med <- lmer(med ~ total_tlss + (1|patient), data = df)
    mod.mean.center <- lmer(mean.center ~ total_tlss + (1|patient), data = df)
    mod.mean <- lmer(mean ~ total_tlss + (1|patient), data = df)
    res <- data.frame("beta.med"= mod.med@beta[2], "beta.med.center" = mod.med.center@beta[2],
                      "beta.mean" = mod.mean@beta[2], "beta.mean.center" = mod.mean.center@beta[2],
                      "top_x"=x, "intersect" = length(topx), inf=list) 
  }
  return(res)
}

## Sliding window --------------------------------------------------------------
# Define your function to process each window with a dynamic meta column
process_window <- function(i, k, inf, vst, meta, score_column) {
  window <- inf[i:(k+(i-1))]
  score <- meta[[score_column]]
  corr.med <- cor(colMedians(vst[window,], useNames = TRUE), score)
  corr.mean <- cor(colMeans(vst[window,]), score)
  return(data.frame(corr.med = corr.med, corr.mean = corr.mean, window_length = length(window), start_pos = i, end_pos = k + (i - 1)))
}

# General function to apply the process
apply_window_process <- function(k, inf, vst, meta, score_column) {
  results <- lapply(1:(length(inf) - (k - 1)), function(i) process_window(i, k, inf, vst, meta, score_column))
  df <- do.call(rbind, results)
  return(df)
}