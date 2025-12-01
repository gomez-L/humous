
# For human, patterned genes only, change the number of neighbors for landscape averaging
  # does this result in differences in correlation results?

source("humous_v3/lib/lib_misc.R")
source("humous_v3/lib/lib_lands.R")

library(Seurat) ; library(SeuratObject) ; library(dplyr) ; library(ggplot2) ; library(png) ; library(readr)

# load data, keep only patterned genes and normalize-scale
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
IntH_ordi <- readRDS("~/humous/humous_v3/out/ordiH/IntH_ordi_new.rds")
landinfo <- read_csv("humous_v4/data/landinfo.csv")
DefaultAssay(IntH_ordi) <- "RNA" ; IntH_ordi[['integrated']] <- NULL 
IntH_ordi_patterned <- subset(IntH_ordi, features= landinfo$gene[landinfo$cond=="H" & !landinfo$AnnotTypes %in% c("Non_patterned","NonSpe")])
IntH_ordi_patterned <- IntH_ordi_patterned %>% NormalizeData() %>% 
  ScaleData(vars.to.regress="nFeature_RNA",features=rownames(IntH_ordi_patterned)) # rang_expr <- range(IntH_ordi@assays$RNA@scale.data)
IntH_ordi_patterned@assays$RNA@scale.data <- IntH_ordi_patterned@assays$RNA@scale.data + abs(min(IntH_ordi_patterned@assays$RNA@scale.data))
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# generate grids and landscapes using 50 and 200 neighbors (normal version has 100)
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

gridlistH_test <- list(grid_medR_50=knn_array_medres(IntH_ordi_patterned$ordi_age_norm,IntH_ordi_patterned$ordi_diff_norm,k=50L),
                  grid_medR_200=knn_array_medres(IntH_ordi_patterned$ordi_age_norm,IntH_ordi_patterned$ordi_diff_norm,k=200L))
#saveRDS(gridlistH_test,"humous_v4/out/knn_averaging_test/gridlistH_test.rds")

L_MR_H_50 <- local({
  m_50 <- IntH_ordi_patterned@assays$RNA@scale.data
  landscapes_MR_H_50 <- knn_rowMeans(m_50,gridlistH_test$grid_medR_50) ; rownames(landscapes_MR_H_50) <- rownames(m_50) ; landscapes_MR_H_50
})
#saveRDS(L_MR_H_50,"humous_v4/out/knn_averaging_test/L_MR_H_50.rds")


L_MR_H_200 <- local({
  m_200 <- IntH_ordi_patterned@assays$RNA@scale.data
  landscapes_MR_H_200 <- knn_rowMeans(m_200,gridlistH_test$grid_medR_200) ; rownames(landscapes_MR_H_200) <- rownames(m_200) ; landscapes_MR_H_200
})
#saveRDS(L_MR_H_200,"humous_v4/out/knn_averaging_test/L_MR_H_200.rds")

# ---# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 



# correlate 50 and 200 results with 100 result and visualize
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

L_MR_H <- readRDS("~/humous/humous_v4/out/landsH/L_MR_H.rds")
L_MR_H_patterned <- L_MR_H[rownames(L_MR_H_200),,]

gridExtra::grid.arrange(as.array(L_MR_H_patterned["SOX2",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster() + ggtitle("SOX2 H knn100"),
                        as.array(L_MR_H_50["SOX2",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("SOX2 H knn50") ,
                        as.array(L_MR_H_200["SOX2",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("SOX2 H knn200"),ncol=3 )

# create empty dataframe with the needed shape
df_corr_nopairwise_50vs100_200vs100 <- data.frame(genename=rownames(L_MR_H_patterned),
                                           corrH_origvs50=NA,corrH_origvs200=NA,stringsAsFactors=FALSE)

# calculate non-pairwise correlations and populate df_corr_nopairwise_origvspe1
for (i in 1:nrow(df_corr_nopairwise_50vs100_200vs100)){
  
  print(i)
  
  df_corr_nopairwise_50vs100_200vs100$corrH_origvs50[i] <- L_TH_BoosCorr(gene1=df_corr_nopairwise_50vs100_200vs100$genename[i],
                                                                   gene2=df_corr_nopairwise_50vs100_200vs100$genename[i],
                                                                   threshold=0.3,n_boosts=10,boostsize=2000,L1=L_MR_H_patterned,L2=L_MR_H_50)$corr
  
  df_corr_nopairwise_50vs100_200vs100$corrH_origvs200[i] <- L_TH_BoosCorr(gene1=df_corr_nopairwise_50vs100_200vs100$genename[i],
                                                                   gene2=df_corr_nopairwise_50vs100_200vs100$genename[i],
                                                                   threshold=0.3,n_boosts=10,boostsize=2000,L1=L_MR_H_patterned,L2=L_MR_H_200)$corr
}

#write.csv(df_corr_nopairwise_50vs100_200vs100,"humous_v4/out/knn_averaging_test/df_corr_nopairwise_50vs100_200vs100.csv")

# density plots to show the distribution of correlations
ggplot(df_corr_nopairwise_50vs100_200vs100) + geom_density(aes(corrH_origvs50)) + theme_minimal()
ggplot(df_corr_nopairwise_50vs100_200vs100) + geom_density(aes(corrH_origvs200)) + theme_minimal()

# see them all
ggplot(df_corr_nopairwise_50vs100_200vs100) + geom_point(aes(corrH_origvs50,corrH_origvs200)) + theme_minimal()

# explore genes with low correlation
gridExtra::grid.arrange(as.array(L_MR_H_patterned["TPPP3",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster() + ggtitle("CSMD3 H Original"),
                        as.array(L_MR_H_50["TPPP3",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("CSMD3 H pe1") ,
                        as.array(L_MR_H_200["TPPP3",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("CSMD3 H pe2"),ncol=3 )

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

