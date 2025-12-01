
# COMPUTER VISION - FEATURE EXTRACTION ON LANDSCAPES USING VGG16 MODEL

# SETUP KERAS
###########################################################################################################################################################################
library(reticulate)
Sys.setenv("RETICULATE_MINICONDA_ENABLED" = TRUE) # make sure that reticulate can use miniconda
reticulate::use_miniconda(condaenv = "minic_py37", required = T) 
reticulate::py_config() # verification
library(keras)
library(tensorflow)
###########################################################################################################################################################################


# LOAD vgg16 MODEL
###########################################################################################################################################################################
model <- keras::application_vgg16(weights="imagenet",include_top = FALSE)
# include_top = FALSE because we don't want to include the last layers that are for classification, we just want the conv layers and the maxpooling for feat extraction
###########################################################################################################################################################################


# load other needed packages
###########################################################################################################################################################################
library(tidyverse) ; library(imager) ; library(magick)  # for preprocessing images
library(raster) ; library(dplyr) ; library(maptools) ; library(png) ; library(rgdal) ; library(microbenchmark)
library(rasterVis) ; library(cowplot) ; library(umap); library(Seurat) ; library(parallel)
###########################################################################################################################################################################


# DEFINE FUNCTION FOR RESHAPING PNGS OF LANDSCAPES TO THE FORMAT NEDDED BY THE vgg16 MODEL 
###########################################################################################################################################################################
image_prep <- function(x) {
  tryCatch(
    {
      img <- readPNG(x)
      x <- image_to_array(img) 
      x <- array_reshape(x, c(1, dim(x)))
    },
    error = function(e){
      return(NA)
    }
  )
} 
###########################################################################################################################################################################


# LAUNCH FEATURE EXTRACTION FOR THE LANDSCAPES USING vgg16 
###########################################################################################################################################################################

# HUMAN
file_list_H <- list.files("out/landscapes_generation/pngLs/ml_pngLs/H/", full.names = TRUE, recursive = FALSE)
vgg16_feature_list_H <- vector(mode="list",length=length(file_list_H))
for (i in 1:length(file_list_H)) {
  print(i)
  vgg16_feature <- predict(model, image_prep(file_list_H[i])) # dim [1,1_7,1:7,1:512]
  flatten <- as.data.frame.table(vgg16_feature, responseName = "value") %>% dplyr::select(value)
  vgg16_feature_list_H[[i]] <-  cbind(i, as.data.frame(t(flatten)))
  if(i %in% c(seq(1000,20000,1000))){gc()}
}
H_vgg16 <- do.call("rbind",vgg16_feature_list_H)
H_vgg16$i <- NULL ; rownames(H_vgg16) <- mgsub::mgsub(file_list_H,pattern=c("out/landscapes_generation/pngLs/ml_pngLs/H//",".png"),replacement=c("",""))
#saveRDS(H_vgg16,"out/computer_vision/H_vgg16.rds")


# MOUSE
file_list_M <- list.files("out/landscapes_generation/pngLs/ml_pngLs/M/", full.names = TRUE, recursive = FALSE)
vgg16_feature_list_M <- vector(mode="list",length=length(file_list_M))
for (i in 1:length(file_list_M)) {
  print(i)
  vgg16_feature <- predict(model, image_prep(file_list_M[i])) # dim [1,1_7,1:7,1:512]
  flatten <- as.data.frame.table(vgg16_feature, responseName = "value") %>% dplyr::select(value)
  vgg16_feature_list_M[[i]] <-  cbind(i, as.data.frame(t(flatten)))
  if(i %in% c(seq(1000,20000,1000))){gc()}
}
M_vgg16 <- do.call("rbind",vgg16_feature_list_M)
M_vgg16$i <- NULL ; rownames(M_vgg16) <- mgsub::mgsub(file_list_M,pattern=c("out/landscapes_generation/pngLs/ml_pngLs/M//",".png"),replacement=c("",""))
#saveRDS(M_vgg16,"out/computer_vision/M_vgg16.rds")


# HORG
file_list_O <- list.files("out/landscapes_generation/pngLs/ml_pngLs/O/", full.names = TRUE, recursive = FALSE)
vgg16_feature_list_O <- vector(mode="list",length=length(file_list_O))
for (i in 1:length(file_list_O)) {
  print(i)
  vgg16_feature <- predict(model, image_prep(file_list_O[i])) # dim [1,1_7,1:7,1:512]
  flatten <- as.data.frame.table(vgg16_feature, responseName = "value") %>% dplyr::select(value)
  vgg16_feature_list_O[[i]] <-  cbind(i, as.data.frame(t(flatten)))
  if(i %in% c(seq(1000,20000,1000))){gc()}
}
O_vgg16 <- do.call("rbind",vgg16_feature_list_O)
O_vgg16$i <- NULL ; rownames(O_vgg16) <- mgsub::mgsub(file_list_O,pattern=c("out/landscapes_generation/pngLs/ml_pngLs/O//",".png"),replacement=c("",""))
#saveRDS(O_vgg16,"out/computer_vision/O_vgg16.rds")

###########################################################################################################################################################################





# LOAD vgg16 EMBEDDINGS AND ASSEMBLE IN SEURAT OBJECT
###########################################################################################################################################################################

Vgg16_S <- local({
  HS <- CreateSeuratObject(t(H_vgg16)) ; HS$dataset <- "human" 
  MS <- CreateSeuratObject(t(M_vgg16)) ; MS$dataset <- "mouse" 
  OS <- CreateSeuratObject(t(O_vgg16)) ; OS$dataset <- "horg" 
  Vgg16_S <- merge(HS,list(MS,OS))
})

Vgg16_S <- Vgg16_S %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA %>% 
  RunUMAP(dims = 1:15,min.dist = 0.01,n.neighbors = 20,n.components=2L,spread=10,local.connectivity=25,verbose=FALSE,return.model=TRUE)
DimPlot(Vgg16_S,group.by = "dataset" )

#saveRDS(Vgg16_S,"out/computer_vision/Vgg16_S.rds")
# NOTE: Vgg16_S.rds contains metadata annotation that can be found at out/computer_vision/Geneinfo.csv

###########################################################################################################################################################################

