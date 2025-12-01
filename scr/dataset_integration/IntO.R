

# INTEGRATION OF ORGHUMAN DATASETS

source("lib/lib_misc.R")
library(Seurat) ; library(dplyr); library(readr) ; library(ggplot2) ; library(plyr) ; library(dplyr) ; library(scales) ; library(purrr)


# load raw org datasets and join them in a list ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
full.list <- list(Av19=readRDS("data/O/Av19.rds"),
                  Kp19=readRDS("data/O/Kp19.rds"),
                  Pt20=readRDS("data/O/Pt20.rds"),
                  Ck19=readRDS("data/O/Ck19.rds"))
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# load HumanInt with selected cells for landscapes ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
IntH <- readRDS("out/ordinals_age_diff/IntH_ordi.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# filter cells on each dataset to keep ages and types of interest, annotate also by dataset and protocol (for integration)
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
DiffTypes_org <- read_csv("data/O//DiffTypes_org.csv")

# Av19
full.list$Av19$age_ek <- DiffTypes_org$age[match(colnames(full.list$Av19),DiffTypes_org$cells)]
full.list$Av19$diff_ek <- DiffTypes_org$DiffTypes[match(colnames(full.list$Av19),DiffTypes_org$cells)]
full.list$Av19$protocol <- DiffTypes_org$proto[match(colnames(full.list$Av19),DiffTypes_org$cells)]
full.list$Av19$dataset <- "Av19"
full.list$Av19$dataset_proto <- paste0("Av19_",full.list$Av19$protocol)
full.list$Av19 <- subset(full.list$Av19,diff_ek %in% c("RG","IPC","N") & age_ek %in% c("1-2 m","2-3 m","3 m","3.5-5 m","6-7 m"))

# Kp19
full.list$Kp19$age_ek <- DiffTypes_org$age[match(colnames(full.list$Kp19),DiffTypes_org$cells)]
full.list$Kp19$diff_ek <- DiffTypes_org$DiffTypes[match(colnames(full.list$Kp19),DiffTypes_org$cells)]
full.list$Kp19$protocol <- DiffTypes_org$proto[match(colnames(full.list$Kp19),DiffTypes_org$cells)]
full.list$Kp19$dataset <- "Kp19"
full.list$Kp19$dataset_proto <- paste0("Kp19_",full.list$Kp19$protocol)
full.list$Kp19 <- subset(full.list$Kp19,diff_ek %in% c("RG","IPC","N") & age_ek %in% c("1-2 m","2-3 m","3 m","3.5-5 m","6-7 m"))

# Ck19
full.list$Ck19$age_ek <- DiffTypes_org$age[match(colnames(full.list$Ck19),DiffTypes_org$cells)]
full.list$Ck19$diff_ek <- DiffTypes_org$DiffTypes[match(colnames(full.list$Ck19),DiffTypes_org$cells)]
full.list$Ck19$protocol <- DiffTypes_org$proto[match(colnames(full.list$Ck19),DiffTypes_org$cells)]
full.list$Ck19$dataset <- "Ck19"
full.list$Ck19$dataset_proto <- paste0("Ck19_",full.list$Ck19$protocol)
full.list$Ck19 <- subset(full.list$Ck19,diff_ek %in% c("RG","IPC","N") & age_ek %in% c("1-2 m","2-3 m","3 m","3.5-5 m","6-7 m"))

# Pt20
full.list$Pt20$age_ek <- DiffTypes_org$age[match(colnames(full.list$Pt20),DiffTypes_org$cells)]
full.list$Pt20$diff_ek <- DiffTypes_org$DiffTypes[match(colnames(full.list$Pt20),DiffTypes_org$cells)]
full.list$Pt20$protocol <- DiffTypes_org$proto[match(colnames(full.list$Pt20),DiffTypes_org$cells)]
full.list$Pt20$dataset <- "Pt20"
full.list$Pt20$dataset_proto <- paste0("Pt20_",full.list$Pt20$protocol)
full.list$Pt20 <- subset(full.list$Pt20,diff_ek %in% c("RG","IPC","N") & age_ek %in% c("1-2 m","2-3 m","3 m","3.5-5 m","6-7 m"))
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #



# CELL SELECTION FOR INTEGRATION ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

mergedO <- merge(full.list$Av19,list(full.list$Kp19,full.list$Pt20,full.list$Ck19))

# SELECTION

mergedO_f1 <- subset(mergedO,  age_ek %in% c("1-2 m","3.5-5 m")  )

mergedO_f2 <- subset(mergedO,  age_ek %in% c("3 m","6-7 m")  )
mergedO_f2$age_dataset <- paste(mergedO_f2$age_ek,mergedO_f2$dataset) ; print(table(mergedO_f2$age_dataset))
mergedO_f2_f1 <- subset(mergedO_f2,age_dataset=="6-7 m Kp19")
mergedO_f2_f2 <- cell.selector(subset(mergedO_f2,age_dataset!="6-7 m Kp19"),
                               colnames(subset(mergedO_f2,age_dataset!="6-7 m Kp19")),
                               subset(mergedO_f2,age_dataset!="6-7 m Kp19")$age_dataset,4322)
mergedO_f2 <- merge(mergedO_f2_f1,mergedO_f2_f2) ; table(mergedO_f2$age_dataset)

mergedO_f3 <- subset(mergedO,  age_ek %in% c("2-3 m"))
mergedO_f3$age_dataset <- paste(mergedO_f3$age_ek,mergedO_f3$dataset) ; print(table(mergedO_f3$age_dataset))
mergedO_f3_f1 <- subset(mergedO_f3,age_dataset %in% c("2-3 m Kp19","2-3 m Ck19"))
mergedO_f3_f2 <- cell.selector(subset(mergedO_f3,age_dataset=="2-3 m Pt20"),
                               colnames(subset(mergedO_f3,age_dataset=="2-3 m Pt20")),
                               subset(mergedO_f3,age_dataset=="2-3 m Pt20")$age_dataset,3639)
mergedO_f3 <- merge(mergedO_f3_f1,mergedO_f3_f2) ; table(mergedO_f3$age_dataset)


# NOW ASSEMBLE ALL SUBSETS
mergedO_f <- merge(mergedO_f1,list(mergedO_f2,mergedO_f3))

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #



# INTEGRATION - HUMAN INT OBJECT AS REFERENCE!!
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# define new seurat object with integrated data assay as counts and data
IntH2 <- CreateSeuratObject(counts=IntH@assays$integrated@data,meta.data = IntH@meta.data)
IntH2@assays$RNA@data <- IntH2@assays$RNA@counts
IntH2$dataset_proto <- "IntH2"
DefaultAssay(IntH2) <- "RNA" ; IntH2[['integrated']] <- NULL

merged <- merge(IntH2,mergedO_f)
merged$dataset_proto2 <- ifelse(grepl("Kp19",merged$dataset_proto),"Kp19",ifelse(grepl("Av19",merged$dataset_proto),"Av19",as.character(merged$dataset_proto)))
table(merged$dataset_proto2,merged$age_ek,useNA = "always")

full.list <- SplitObject(merged,split.by =  "dataset_proto2") 


for (i in 2:length(full.list)) {
  full.list[[i]] <- NormalizeData(full.list[[i]], verbose = FALSE)
}

# find variable features in all (although in human we'll consider only the 2000 used previously - integrated genes)
for (i in 1:length(full.list)) {
  full.list[[i]] <- FindVariableFeatures(full.list[[i]], selection.method = "mvp",nfeatures=10000, verbose = FALSE)
}

# find which genes of human integration are also present in all org datasets
features <- Reduce(intersect,list(rownames(IntH2),rownames(full.list$Av19),rownames(full.list$Ck19),rownames(full.list$Kp19),rownames(full.list$Pt20_Xiang)))

# run integration of integration
IntHO <- double_dataset_integration(list_seuratobjs=full.list,nfeats=10000,FindIntDims=1:15,FindInt_kScore=20,IntKweight=20,reference_dataset="IntH2",features=features)

# scale nFeature_RNA and calculate UMAP 
DefaultAssay(IntHO) <- "integrated" 
IntHO <- IntHO %>% ScaleData(vars.to.regress="nFeature_RNA") %>% RunPCA %>% RunUMAP(dims = 1:25,return.model = TRUE,min.dist = 0.5,n.neighbors = 50,local.connectivity = 6)
DimPlot(IntHO,group.by = "dataset_proto2",pt.size = 0.5,split.by = "dataset" )
DimPlot(IntHO,group.by = "age_ek",pt.size = 0.5,split.by = "dataset_proto2")
DimPlot(IntHO,group.by = "diff_ek",pt.size = 0.5,split.by = "dataset_proto2") 

# save result
#saveRDS(IntHO,"out/dataset_integration/IntHO.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #






