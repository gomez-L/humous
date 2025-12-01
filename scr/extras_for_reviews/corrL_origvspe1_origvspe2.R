

# CORRELATION ANALYSIS ON RPH analysis - original landscapes vs pe1, original landscapes vs pe1, only patterned genes

library(dplyr) ; library(parallel) ; library(png) ; library(readr) ; library(ggplot2) ; library(plotly)
source("humous_v3/lib/lib_lands.R")

# load data
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
L_MR_H <- readRDS("~/humous/humous_v4/out/landsH/L_MR_H.rds")
L_MR_H_pe1 <- readRDS("humous_v4/out/raw_permuted_humous/RPH_pe1/LandsH_pe1/L_MR_H_pe1.rds")
L_MR_H_pe2 <- readRDS("humous_v4/out/raw_permuted_humous/RPH_pe2/LandsH_pe2/L_MR_H_pe2.rds")
# subset original landscapes array to keep only patterned genes, as in permutation arrays
L_MR_H_patterned <- L_MR_H[rownames(L_MR_H_pe1),,]
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# correlations original lands vs pe1 and orig lands vs pe2
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

# create empty dataframe with the needed shape
df_corr_nopairwise_origvspe1 <- data.frame(genename=rownames(L_MR_H_patterned),
                                 corrH_origvspe1=NA,corrH_origvspe2=NA,stringsAsFactors=FALSE)

# calculate non-pairwise correlations and populate df_corr_nopairwise_origvspe1
for (i in 1:nrow(df_corr_nopairwise_origvspe1)){
  
  print(i)

  df_corr_nopairwise_origvspe1$corrH_origvspe1[i] <- L_TH_BoosCorr(gene1=df_corr_nopairwise_origvspe1$genename[i],
                                                                   gene2=df_corr_nopairwise_origvspe1$genename[i],
                                                   threshold=0.3,n_boosts=10,boostsize=2000,L1=L_MR_H_patterned,L2=L_MR_H_pe1)$corr
  
  df_corr_nopairwise_origvspe1$corrH_origvspe2[i] <- L_TH_BoosCorr(gene1=df_corr_nopairwise_origvspe1$genename[i],
                                                                   gene2=df_corr_nopairwise_origvspe1$genename[i],
                                                                   threshold=0.3,n_boosts=10,boostsize=2000,L1=L_MR_H_patterned,L2=L_MR_H_pe2)$corr
}

#write.csv(df_corr_nopairwise_origvspe1,"humous_v4/out/raw_permuted_humous/df_corr_nopairwise_origvspe1_origvspe2.csv")

# density plots to show the distribution of correlations
ggplot(df_corr_nopairwise_origvspe1) + geom_density(aes(corrH_origvspe1)) + theme_minimal()
ggplot(df_corr_nopairwise_origvspe1) + geom_density(aes(corrH_origvspe2)) + theme_minimal()

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 



# explore results - which genes have negative, zero, or positive correlations? any pattern?
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

df_corr_nopairwise_origvspe1_origvspe2 <- read_csv("humous_v4/out/raw_permuted_humous/df_corr_nopairwise_origvspe1_origvspe2.csv")

# do correlation values change a lot between pe1 and pe2?
ggplot(df_corr_nopairwise_origvspe1_origvspe2) + geom_point(aes(corrH_origvspe1,corrH_origvspe2)) + theme_minimal()


gridExtra::grid.arrange(as.array(L_MR_H_patterned["CSMD3",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster() + ggtitle("CSMD3 H Original"),
                        as.array(L_MR_H_pe1["CSMD3",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("CSMD3 H pe1") ,
                        as.array(L_MR_H_pe2["CSMD3",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("CSMD3 H pe2"),ncol=3 )

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 






