

# LANDSCAPES' CORRELATIONS

library(dplyr) ; library(parallel) ; library(png) ; library(readr) ; library(ggplot2) ; library(plotly)

source("lib/lib_lands.R")

# Load landscapes data
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
L_MR_M <- readRDS("out/landscapes_generation/M/L_MR_M.rds")
L_MR_H <- readRDS("out/landscapes_generation/H/L_MR_H.rds")
L_MR_O <- readRDS("out/landscapes_generation/O/L_MR_O.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# NON-PAIRWWISE CORRELATIONS
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

# load orthologs table
ortholog <- as.data.frame(read.csv("data/MH_orthologs.csv"))
# check which genes have a exact same across species (besides capitalization)
ortholog$samename <- ifelse(tolower(ortholog$mmu)==tolower(ortholog$hsa),"yes","no")

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 



# CALCULATE NON-PAIRWISE CORRELATIONS - SAME GENE ACROSS SPECIES CORRELATIONS 
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

# create empty dataframe with the needed shape
df_corr_nopairwise <- data.frame(genename_h=c(ortholog$hsa[ortholog$samename=="yes"]),genename_m=c(ortholog$mmu[ortholog$samename=="yes"]),
                                 corr_HvsM=NA,corr_HvsO=NA,corr_MvsO=NA,stringsAsFactors=FALSE)

# calculate non-pairwise correlations and populate df_corr_nopairwise
for (i in 1:nrow(df_corr_nopairwise)){
  
  print(df_corr_nopairwise$genename_h[i])
  
  if( any(df_corr_nopairwise$genename_h[i] %in% rownames(L_MR_H)) & any(df_corr_nopairwise$genename_m[i] %in% rownames(L_MR_M)) ){
    df_corr_nopairwise$corr_HvsM[i] <- L_TH_BoosCorr(gene1=df_corr_nopairwise$genename_h[i],gene2=df_corr_nopairwise$genename_m[i],
                                                     threshold=0.3,n_boosts=10,boostsize=2000,L1=L_MR_H,L2=L_MR_M)$corr
  }else{"do nothing"}
  
  if( any(df_corr_nopairwise$genename_h[i] %in% rownames(L_MR_H)) & any(df_corr_nopairwise$genename_h[i] %in% rownames(L_MR_O)) ){
    df_corr_nopairwise$corr_HvsO[i] <- L_TH_BoosCorr(gene1=df_corr_nopairwise$genename_h[i],gene2=df_corr_nopairwise$genename_h[i],
                                                     threshold=0.3,n_boosts=10,boostsize=2000,L1=L_MR_H,L2=L_MR_O)$corr
  }else{"do nothing"}
  
  if( any(df_corr_nopairwise$genename_m[i] %in% rownames(L_MR_M)) & any(df_corr_nopairwise$genename_h[i] %in% rownames(L_MR_O)) ){
    df_corr_nopairwise$corr_MvsO[i] <- L_TH_BoosCorr(gene1=df_corr_nopairwise$genename_m[i],gene2=df_corr_nopairwise$genename_h[i],
                                                     threshold=0.3,n_boosts=10,boostsize=2000,L1=L_MR_M,L2=L_MR_O)$corr
  }else{"do nothing"}

}

saveRDS(df_corr_nopairwise,"out/correlations/df_corr_nopairwise.rds")

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# PAIRWWISE CORRELATIONS ARE ONLY TO DISPLAY IN WEBSITE

