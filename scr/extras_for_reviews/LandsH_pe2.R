


# LANDSCAPES CALCULATION HUMAN - raw permutation humous - permutation 2 human - only patterned genes

source("humous_v3/lib/lib_misc.R")
source("humous_v3/lib/lib_lands.R")

library(Seurat) ; library(SeuratObject) ; library(dplyr) ; library(ggplot2) ; library(png)

# Load data (result from integration and ordis) - 
# contains cells for landscapes and ordis normalized
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
IntH_pe2_ordi <- readRDS("~/humous/humous_v4/out/raw_permuted_humous/RPH_pe2/ordiH_pe2/IntH_pe2_ordi.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# Preprocess datasets for landscapes
# Keep only patterned genes, use RNA assay
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

landinfo <- read_csv("humous_v4/data/landinfo.csv")
DefaultAssay(IntH_pe2_ordi) <- "RNA" ; IntH_pe2_ordi[['integrated']] <- NULL 
table(landinfo$AnnotTypes,landinfo$cond)
IntH_pe2_ordi_patterned <- subset(IntH_pe2_ordi, features= landinfo$gene[landinfo$cond=="H" & !landinfo$AnnotTypes %in% c("Non_patterned","NonSpe")])
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 



# DATA SCALING AND SAVING
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
# scale data slot for nfeatures on all genes and correct it so it doesnt have negative values
LandsS_H_pe2 <- IntH_pe2_ordi_patterned %>% NormalizeData() %>% 
  ScaleData(vars.to.regress="nFeature_RNA",features=rownames(IntH_pe2_ordi_patterned)) 

LandsS_H_pe2@assays$RNA@scale.data <- LandsS_H_pe2@assays$RNA@scale.data + abs(min(LandsS_H_pe2@assays$RNA@scale.data))
#saveRDS(LandsS_H_pe2,"humous_v4/out/raw_permuted_humous/RPH_pe2/LandsH_pe2/LandsS_H_pe2.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# GRIDS generation
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
gridlistH_pe2 <- list(grid_medR=knn_array_medres(LandsS_H_pe2$ordi_age_norm,LandsS_H_pe2$ordi_diff_norm,k=100L))
#saveRDS(gridlistH_pe2,"humous_v4/out/raw_permuted_humous/RPH_pe2/LandsH_pe2/gridlistH_pe2.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# COMPUTE LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

# trial for key genes
L_MR_H_pe2 <- local({
  m <- LandsS_H_pe2@assays$RNA@scale.data[c("SOX2","EOMES","NEUROD6","JUNB","APBB2","MEF2C"),]
  landscapes_MR_H <- knn_rowMeans(m,gridlistH_pe2$grid_medR) ; rownames(landscapes_MR_H) <- rownames(m) ; landscapes_MR_H
})

gridExtra::grid.arrange(as.array(L_MR_H_pe2["SOX2",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster() + ggtitle("SOX2"),
                        as.array(L_MR_H_pe2["EOMES",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("EOMES") ,
                        as.array(L_MR_H_pe2["NEUROD6",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("NEUROD6"), 
                        as.array(L_MR_H_pe2["JUNB",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster() + ggtitle("JUNB"),
                        as.array(L_MR_H_pe2["APBB2",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("APBB2") ,
                        as.array(L_MR_H_pe2["MEF2C",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("MEF2C"),ncol=3 )


# MEDIUM RESOLUTION, ALL PATTERNED GENES
L_MR_H_pe2 <- local({
  m <- LandsS_H_pe2@assays$RNA@scale.data
  landscapes_MR_H <- knn_rowMeans(m,gridlistH_pe2$grid_medR) ; rownames(landscapes_MR_H) <- rownames(m) ; landscapes_MR_H
})
#saveRDS(L_MR_H_pe2,"humous_v4/out/raw_permuted_humous/RPH_pe2/LandsH_pe2/L_MR_H_pe2.rds")

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 