

# LANDSCAPES GENERATION MOUSE

source("lib/lib_misc.R")
source("lib/lib_lands.R")

library(Seurat) ; library(SeuratObject) ; library(dplyr) ; library(ggplot2) ; library(stringr) ; library(png)

# Load data (result from integration and ordis) - contains cells for landscapes and ordis normalized
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
IntM_ordi <- readRDS("out/ordinals_age_diff/M/IntM_ordi.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# Preprocess dataset for landscapes
DefaultAssay(IntM_ordi) <- "RNA" ; IntM_ordi[['integrated']] <- NULL 

# GENE SELECTION
###################
# subset to keep just genes expressed in more than 1% of the cells across datasets 
num_positive <- apply(IntM_ordi@assays$RNA@counts > 0, 1, sum)
genestokeep_m <- rownames(IntM_ordi)[num_positive > round(1/100 * ncol(IntM_ordi)) ] # 1%
LandsS_M <- subset(IntM_ordi,features = genestokeep_m)
###################


# DATA SCALING AND SAVING
###################
# scale data slot for nfeatures on all genes and correct it so it doesnt have negative values
LandsS_M <- LandsS_M %>% NormalizeData() %>% ScaleData(vars.to.regress="nFeature_RNA",features=rownames(LandsS_M)) # rang_expr <- range(IntM_ordi@assays$RNA@scale.data)
LandsS_M@assays$RNA@scale.data <- LandsS_M@assays$RNA@scale.data + abs(min(LandsS_M@assays$RNA@scale.data))
#saveRDS(LandsS_M,"out/landscapes_generation/M/LandsS_M.rds")
###################


# VISUALIZATION
LandsS_M$Sox2 <- LandsS_M@assays$RNA@scale.data["Sox2",] ; LandsS_M$Eomes <- LandsS_M@assays$RNA@scale.data["Eomes",] ; LandsS_M$Neurod6 <- LandsS_M@assays$RNA@scale.data["Neurod6",]
LandsS_M$Junb <- LandsS_M@assays$RNA@scale.data["Junb",] ; LandsS_M$Apbb2 <- LandsS_M@assays$RNA@scale.data["Apbb2",] ; LandsS_M$Mef2c <- LandsS_M@assays$RNA@scale.data["Mef2c",]
ggplot(LandsS_M@meta.data) + geom_point(aes(ordi_age_norm,ordi_diff_norm,color=Eomes)) + scale_color_gradient2(low="white",mid="grey",high="red",midpoint = 6)


# GRIDS generation
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
gridlistM <- list(grid_medR=knn_array_medres(LandsS_M$ordi_age_norm,LandsS_M$ordi_diff_norm,k=100L),
                 grid_lowR=knn_array_lowres(LandsS_M$ordi_age_norm,LandsS_M$ordi_diff_norm,k=100L),
                 grid_highR=knn_array_highres(LandsS_M$ordi_age_norm,LandsS_M$ordi_diff_norm,k=100L))
#saveRDS(gridlistM,"out/landscapes_generation/M/gridlistM.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# COMPUTE LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

# MEDIUM RESOLUTION, ALL GENES
L_MR_M <- local({
  m <- LandsS_M@assays$RNA@scale.data
  landscapes_MR_M <- knn_rowMeans(m,gridlistM$grid_medR) ; rownames(landscapes_MR_M) <- rownames(m) ; landscapes_MR_M
})
#saveRDS(L_MR_M,"out/landscapes_generation/M/L_MR_M.rds")

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# STANDARD DEVIATION LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
sd_L_M_df <- data.frame(gene=rownames(L_MR_M),sd_L_M=NA)
sd_L_M <- function(gene){sd(as.vector(L_MR_M[gene,,] ) ) }
a <- mapply(sd_L_M,sd_L_M_df$gene)
sd_L_M_df$sd_L_M <- a[match(sd_L_M_df$gene,names(a))]
#saveRDS(sd_L_M_df,"out/landscapes_generation/M/sd_L_M_df.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# PNG LANDSCAPES FOR WEBSITE
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

parallel::detectCores()
mclapply(mc.cores=64, rownames(L_MR_M),function(i){
  png(paste0("out/landscapes_generation/pngLs/web_pngLs/M/",i,"_M.png"),width=500,height=500) 
  print(as.array(L_MR_M[i,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% website_ggraster() + ggtitle(paste(str_to_title(i)," (Ms)")) + 
          xlab("Age") + ylab("Differentiation") + theme(axis.title = element_text(size=40)))
  dev.off()
})
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 



# PNG LANDSCAPES FOR DEEPLEARNING
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

parallel::detectCores()
mclapply(mc.cores=16, rownames(L_MR_M),function(i){
  png(paste0("out/landscapes_generation/pngLs/ml_pngLs/M/",i,"_M.png"),width=224,height=224) 
  print(as.array(L_MR_M[i,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% ML_ggraster())
  dev.off()
})

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# ENTROPY LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

pngs_M <- list.files("humous_v4/out/pngLs/ml_pngLs/M/",full.names = TRUE,recursive = TRUE,pattern = ".png")
ent_L_M_df <- data.frame(gene=mgsub::mgsub(pattern=c("humous_v4/out/pngLs/ml_pngLs/M//",".png"),replacement=c("",""),string=pngs_M),entropy=NA)
for (i in levels(as.factor(pngs_M))){
  ent_L_M_df$entropy[ent_L_M_df$gene==mgsub::mgsub(pattern=c("humous_v4/out/pngLs/ml_pngLs/M//",".png"),replacement=c("",""),string=i)] <- entropy::entropy(as.matrix(readPNG(i)),method="CS")
}
#saveRDS(ent_L_M_df,"out/landscapes_generation/M/ent_L_M_df.rds")

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 










