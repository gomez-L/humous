

# LANDSCAPES GENERATION HUMAN

source("lib/lib_misc.R")
source("lib/lib_lands.R")

library(Seurat) ; library(SeuratObject) ; library(dplyr) ; library(ggplot2) ; library(png)

# Load data (result from integration and ordis) - contains cells for landscapes and ordis normalized
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
IntH_ordi <- readRDS("out/ordinals_age_diff/H/IntH_ordi.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# Preprocess datasets for landscapes
# Using RNA assay, keep genes expressed in all datasets and scale data slot of RNA for nfeatures
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

DefaultAssay(IntH_ordi) <- "RNA" ; IntH_ordi[['integrated']] <- NULL 

# GENE SELECTION
###################

# subset to keep just genes expressed in more than 1% of the cells across datasets 
num_positive <- apply(IntH_ordi@assays$RNA@counts > 0, 1, sum)
genestokeep_h <- rownames(IntH_ordi)[num_positive > round(1/100 * ncol(IntH_ordi)) ] # 1%

LandsS_H <- subset(IntH_ordi,features = genestokeep_h)
###################


# DATA SCALING AND SAVING
###################
# scale data slot for nfeatures on all genes and correct it so it doesnt have negative values
LandsS_H <- LandsS_H %>% NormalizeData() %>% ScaleData(vars.to.regress="nFeature_RNA",features=rownames(LandsS_H)) # rang_expr <- range(IntH_ordi@assays$RNA@scale.data)
LandsS_H@assays$RNA@scale.data <- LandsS_H@assays$RNA@scale.data + abs(min(LandsS_H@assays$RNA@scale.data))
#saveRDS(LandsS_H,"out/landscapes_generation/H/LandsS_H.rds")
###################


# VISUALIZATION
LandsS_H$SOX2 <- LandsS_H@assays$RNA@scale.data["SOX2",] ; LandsS_H$EOMES <- LandsS_H@assays$RNA@scale.data["EOMES",] ; LandsS_H$NEUROD6 <- LandsS_H@assays$RNA@scale.data["NEUROD6",]
LandsS_H$JUNB <- LandsS_H@assays$RNA@scale.data["JUNB",] ; LandsS_H$APBB2 <- LandsS_H@assays$RNA@scale.data["APBB2",] ; LandsS_H$MEF2C <- LandsS_H@assays$RNA@scale.data["MEF2C",]
ggplot(LandsS_H@meta.data) + geom_point(aes(ordi_age_norm,ordi_diff_norm,color=EOMES)) + scale_color_gradient2(low="white",mid="grey",high="red",midpoint = 6)

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# GRIDS generation
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
gridlistH <- list(grid_medR=knn_array_medres(LandsS_H$ordi_age_norm,LandsS_H$ordi_diff_norm,k=100L),
                  grid_lowR=knn_array_lowres(LandsS_H$ordi_age_norm,LandsS_H$ordi_diff_norm,k=100L),
                  grid_highR=knn_array_highres(LandsS_H$ordi_age_norm,LandsS_H$ordi_diff_norm,k=100L))
#saveRDS(gridlistH,"out/landscapes_generation/H/gridlistH.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# COMPUTE LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

# MEDIUM RESOLUTION, ALL GENES
L_MR_H <- local({
  m <- LandsS_H@assays$RNA@scale.data
  landscapes_MR_H <- knn_rowMeans(m,gridlistH$grid_medR) ; rownames(landscapes_MR_H) <- rownames(m) ; landscapes_MR_H
})
#saveRDS(L_MR_H,"out/landscapes_generation/H/L_MR_H.rds")

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# STANDARD DEVIATION LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
sd_L_H_df <- data.frame(gene=rownames(L_MR_H),sd_L_H=NA)
sd_L_H <- function(gene){sd(as.vector(L_MR_H[gene,,] ) ) }
a <- mapply(sd_L_H,sd_L_H_df$gene)
sd_L_H_df$sd_L_H <- a[match(sd_L_H_df$gene,names(a))]
#saveRDS(sd_L_H_df,"out/landscapes_generation/H/sd_L_H_df.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# PNG LANDSCAPES FOR WEBSITE
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

parallel::detectCores()
mclapply(mc.cores=64, rownames(L_MR_H),function(i){
  png(paste0("out/landscapes_generation/pngLs/web_pngLs/H/",i,"_H.png"),width=500,height=500) 
  print(as.array(L_MR_H[i,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% website_ggraster() + ggtitle(paste(i," (Hm)")) + 
          xlab("Age") + ylab("Differentiation") + theme(axis.title = element_text(size=40)))
  dev.off()
})

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# PNG LANDSCAPES FOR DEEPLEARNING
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

parallel::detectCores()
mclapply(mc.cores=16, rownames(L_MR_H),function(i){
  png(paste0("out/landscapes_generation/pngLs/ml_pngLs/H/",i,"_H.png"),width=224,height=224) 
  print(as.array(L_MR_H[i,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% ML_ggraster())
  dev.off()
})

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# ENTROPY LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
pngs_H <- list.files("out/landscapes_generation/pngLs/ml_pngLs/H/",full.names = TRUE,recursive = TRUE,pattern = ".png")
ent_L_H_df <- data.frame(gene=mgsub::mgsub(pattern=c("out/landscapes_generation/pngLs/ml_pngLs/H//",".png"),replacement=c("",""),string=pngs_H),entropy=NA)
for (i in levels(as.factor(pngs_H))){
  print(i)
  ent_L_H_df$entropy[ent_L_H_df$gene==mgsub::mgsub(pattern=c("out/landscapes_generation/pngLs/ml_pngLs/H//",".png"),replacement=c("",""),string=i)] <- entropy::entropy(as.matrix(readPNG(i)),method="CS")
}
#saveRDS(ent_L_H_df,"out/landscapes_generation/H/ent_L_H_df.rds")

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 











