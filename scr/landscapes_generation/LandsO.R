

# LANDSCAPES GENRATION HUMAN ORGS

source("lib/lib_misc.R")
source("lib/lib_lands.R")

library(Seurat) ; library(SeuratObject) ; library(dplyr) ; library(ggplot2) ; library(png)

# Load data (result from integration and ordis) - 
# contains cells for landscapes and ordis normalized
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
IntO_ordi <- readRDS("out/ordinals_age_diff/O/IntO_ordi.rds")


# Preprocess datasets for landscapes
# assemble and annotate (go back to raw data to make sure cells selected are correct)
DefaultAssay(IntO_ordi) <- "RNA" ; IntO_ordi[['integrated']] <- NULL
IntO_ordi_f <- subset(IntO_ordi,dataset %in% c("Ck19","Kp19"))
Pt20 <- readRDS("data/raw/Pt20/Pt20.rds") ; Pt20 <- subset(Pt20,cells=intersect(colnames(Pt20),colnames(IntO_ordi))) 
Av19 <- readRDS("data/raw/Av19/Av19.rds") ; Av19 <- subset(Av19,cells=intersect(colnames(Av19),colnames(IntO_ordi))) 
IntO_ordi_2 <- merge(IntO_ordi_f,list(Pt20,Av19))
IntO_ordi_2$dataset <- IntO_ordi$dataset[match(colnames(IntO_ordi_2),colnames(IntO_ordi))] 
IntO_ordi_2$dataset_proto <- IntO_ordi$dataset_proto[match(colnames(IntO_ordi_2),colnames(IntO_ordi))] 
IntO_ordi_2$dataset_proto2 <- IntO_ordi$dataset_proto2[match(colnames(IntO_ordi_2),colnames(IntO_ordi))] ; table(IntO_ordi_2$dataset_proto2, useNA="always" )
IntO_ordi_2$ordi_age <- IntO_ordi$ordi_age[match(colnames(IntO_ordi_2),colnames(IntO_ordi))] 
IntO_ordi_2$ordi_age_norm <- IntO_ordi$ordi_age_norm[match(colnames(IntO_ordi_2),colnames(IntO_ordi))] 
IntO_ordi_2$ordi_diff <- IntO_ordi$ordi_diff[match(colnames(IntO_ordi_2),colnames(IntO_ordi))] 
IntO_ordi_2$ordi_diff_norm <- IntO_ordi$ordi_diff_norm[match(colnames(IntO_ordi_2),colnames(IntO_ordi))] 
IntO_ordi_2$diff_ek <- IntO_ordi$diff_ek[match(colnames(IntO_ordi_2),colnames(IntO_ordi))] 
IntO_ordi_2$age_ek <- IntO_ordi$age_ek[match(colnames(IntO_ordi_2),colnames(IntO_ordi))] 
IntO_ordi_2$y_age <- IntO_ordi$y_age[match(colnames(IntO_ordi_2),colnames(IntO_ordi))]  
IntO_ordi_2$y_diff <- IntO_ordi$y_diff[match(colnames(IntO_ordi_2),colnames(IntO_ordi))] 
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# GENE SELECTION
###################
# subset to keep just genes expressed in more than 1% of the cells across datasets 
num_positive <- apply(IntO_ordi_2@assays$RNA@counts > 0, 1, sum)
genestokeep_o <- rownames(IntO_ordi_2)[num_positive > round(1/100 * ncol(IntO_ordi_2)) ] # 1%
LandsS_O <- subset(IntO_ordi_2,features = genestokeep_o)
###################


# DATA SCALING AND SAVING
###################
# scale data slot for nfeatures on all genes and correct it so it doesnt have negative values
LandsS_O <- LandsS_O %>% NormalizeData() %>% ScaleData(vars.to.regress="nFeature_RNA",features=rownames(LandsS_O)) # rang_expr <- range(IntH_ordi@assays$RNA@scale.data)
LandsS_O@assays$RNA@scale.data <- LandsS_O@assays$RNA@scale.data + abs(min(LandsS_O@assays$RNA@scale.data))
#saveRDS(LandsS_O,"out/landscapes_generation/O/LandsS_O.rds")
###################


# VISUALIZATION
LandsS_O$SOX2 <- LandsS_O@assays$RNA@scale.data["SOX2",] ; LandsS_O$EOMES <- LandsS_O@assays$RNA@scale.data["EOMES",] ; LandsS_O$NEUROD6 <- LandsS_O@assays$RNA@scale.data["NEUROD6",]
LandsS_O$JUNB <- LandsS_O@assays$RNA@scale.data["JUNB",] ; LandsS_O$APBB2 <- LandsS_O@assays$RNA@scale.data["APBB2",] ; LandsS_O$MEF2C <- LandsS_O@assays$RNA@scale.data["MEF2C",]
ggplot(LandsS_O@meta.data) + geom_point(aes(ordi_age_norm,ordi_diff_norm,color=NEUROD6)) + scale_color_gradient2(low="white",mid="grey",high="red",midpoint = 6)

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# GRIDS generation
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
gridlistO <- list(grid_medR=knn_array_medres(LandsS_O$ordi_age_norm,LandsS_O$ordi_diff_norm,k=100L),
                  grid_lowR=knn_array_lowres(LandsS_O$ordi_age_norm,LandsS_O$ordi_diff_norm,k=100L),
                  grid_highR=knn_array_highres(LandsS_O$ordi_age_norm,LandsS_O$ordi_diff_norm,k=100L))
#saveRDS(gridlistO,"out/landscapes_generation/O/gridlistO.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# COMPUTE LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

# MEDIUM RESOLUTION, ALL GENES
L_MR_O <- local({
  m <- LandsS_O@assays$RNA@scale.data
  landscapes_MR_O <- knn_rowMeans(m,gridlistO$grid_medR) ; rownames(landscapes_MR_O) <- rownames(m) ; landscapes_MR_O
})
#saveRDS(L_MR_O,"out/landscapes_generation/O/L_MR_O.rds")

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# STANDARD DEVIATION LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
sd_L_O_df <- data.frame(gene=rownames(L_MR_O),sd_L_O=NA)
sd_L_O <- function(gene){sd(as.vector(L_MR_O[gene,,] ) ) }
a <- mapply(sd_L_O,sd_L_O_df$gene)
sd_L_O_df$sd_L_O <- a[match(sd_L_O_df$gene,names(a))]
#saveRDS(sd_L_O_df,"out/landscapes_generation/O/sd_L_O_df.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

# PNG LANDSCAPES FOR WEBSITE
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

parallel::detectCores()
mclapply(mc.cores=64, rownames(L_MR_O),function(i){
  png(paste0("out/landscapes_generation/pngLs/web_pngLs/O/",i,"_O.png"),width=500,height=500) 
  print(as.array(L_MR_O[i,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% website_ggraster() + ggtitle(paste(i," (Org)")) + 
          xlab("Age") + ylab("Differentiation") + theme(axis.title = element_text(size=40)) )
  dev.off()
})
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# PNG LANDSCAPES FOR DEEPLEARNING
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 

parallel::detectCores()
mclapply(mc.cores=16, rownames(L_MR_O),function(i){
  png(paste0("out/landscapes_generation/pngLs/ml_pngLs/O/",i,"_O.png"),width=224,height=224) 
  print(as.array(L_MR_O[i,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% ML_ggraster())
  dev.off()
})

npngs <- list.files("out/landscapes_generation/pngLs/ml_pngLs/O/")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 


# ENTROPY LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 
pngs_M <- list.files("out/landscapes_generation/pngLs/ml_pngLs/O/",full.names = TRUE,recursive = TRUE,pattern = ".png")
ent_L_M_df <- data.frame(gene=mgsub::mgsub(pattern=c("out/landscapes_generation/pngLs/ml_pngLs/M//",".png"),replacement=c("",""),string=pngs_M),entropy=NA)
for (i in levels(as.factor(pngs_M))){
  ent_L_M_df$entropy[ent_L_M_df$gene==mgsub::mgsub(pattern=c("out/landscapes_generation/pngLs/ml_pngLs/M//",".png"),replacement=c("",""),string=i)] <- entropy::entropy(as.matrix(readPNG(i)),method="CS")
}
#saveRDS(ent_L_M_df,"out/landscapes_generation/O/ent_L_M_df.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- 







