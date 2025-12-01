


# LANDSCAPES CALCULATION HUMAN - raw permutation humous - permutation 1 human - only patterned genes

source("humous_v3/lib/lib_misc.R")
source("humous_v3/lib/lib_lands.R")

library(Seurat) ; library(SeuratObject) ; library(dplyr) ; library(ggplot2) ; library(png)

# Load data (result from integration and ordis) - 
# contains cells for landscapes and ordis normalized
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
IntH_pe1_ordi <- readRDS("~/humous/humous_v4/out/raw_permuted_humous/RPH_pe1/ordiH_pe1/IntH_pe1_ordi.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# Preprocess datasets for landscapes
# Keep only patterned genes, use RNA assay
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

landinfo <- read_csv("humous_v4/data/landinfo.csv")
DefaultAssay(IntH_pe1_ordi) <- "RNA" ; IntH_pe1_ordi[['integrated']] <- NULL 
table(landinfo$AnnotTypes,landinfo$cond)
IntH_pe1_ordi_patterned <- subset(IntH_pe1_ordi, features= landinfo$gene[landinfo$cond=="H" & !landinfo$AnnotTypes %in% c("Non_patterned","NonSpe")])
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 



# DATA SCALING AND SAVING
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
# scale data slot for nfeatures on all genes and correct it so it doesnt have negative values
LandsS_H_pe1 <- IntH_pe1_ordi_patterned %>% NormalizeData() %>% 
  ScaleData(vars.to.regress="nFeature_RNA",features=rownames(IntH_pe1_ordi_patterned)) 

LandsS_H_pe1@assays$RNA@scale.data <- LandsS_H_pe1@assays$RNA@scale.data + abs(min(LandsS_H_pe1@assays$RNA@scale.data))
#saveRDS(LandsS_H_pe1,"humous_v4/out/raw_permuted_humous/RPH_pe1/LandsH_pe1/LandsS_H_pe1.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# GRIDS generation
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
gridlistH_pe1 <- list(grid_medR=knn_array_medres(LandsS_H_pe1$ordi_age_norm,LandsS_H_pe1$ordi_diff_norm,k=100L))
#saveRDS(gridlistH_pe1,"humous_v4/out/raw_permuted_humous/RPH_pe1/LandsH_pe1/gridlistH_pe1.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# COMPUTE LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

# trial for key genes
L_MR_H_pe1 <- local({
  m <- LandsS_H_pe1@assays$RNA@scale.data[c("SOX2","EOMES","NEUROD6","JUNB","APBB2","MEF2C"),]
  landscapes_MR_H <- knn_rowMeans(m,gridlistH_pe1$grid_medR) ; rownames(landscapes_MR_H) <- rownames(m) ; landscapes_MR_H
})

gridExtra::grid.arrange(as.array(L_MR_H_pe1["SOX2",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster() + ggtitle("SOX2"),
                        as.array(L_MR_H_pe1["EOMES",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("EOMES") ,
                        as.array(L_MR_H_pe1["NEUROD6",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("NEUROD6"), 
                        as.array(L_MR_H_pe1["JUNB",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster() + ggtitle("JUNB"),
                        as.array(L_MR_H_pe1["APBB2",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("APBB2") ,
                        as.array(L_MR_H_pe1["MEF2C",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster()+ ggtitle("MEF2C"),ncol=3 )


# MEDIUM RESOLUTION, ALL PATTERNED GENES
L_MR_H_pe1 <- local({
  m <- LandsS_H_pe1@assays$RNA@scale.data
  landscapes_MR_H <- knn_rowMeans(m,gridlistH_pe1$grid_medR) ; rownames(landscapes_MR_H) <- rownames(m) ; landscapes_MR_H
})
#saveRDS(L_MR_H_pe1,"humous_v4/out/raw_permuted_humous/RPH_pe1/LandsH_pe1/L_MR_H_pe1.rds")

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
