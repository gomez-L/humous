
# INTEGRATION OF MOUSE DATASETS

source("lib/lib_misc.R")
library(Seurat) ; library(dplyr); library(readr) ; library(ggplot2) ; library(plyr) ; library(dplyr) ; library(scales)

      
# load raw mouse datasets and join them in a list ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
full.list <- list(Ad21=readRDS("data/M/Ad21.rds"),
                  My17=readRDS("data/M/My17.rds"),
                  LMl21=readRDS("data/M/LMl21.rds"),
                  Jt19=readRDS("data/M/Jt19.rds"))
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# filter cells on each dataset to keep ages and types of interest, annotate also by dataset (for integration)
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# Ad21
annot_Ad21 <- read_csv("data/M/DiffTypes_Ad21.csv") ; table(annot_Ad21$age,annot_Ad21$DiffTypes,useNA="always") ; dim(annot_Ad21)
full.list$Ad21 <- subset(full.list$Ad21,cells=annot_Ad21$cells[annot_Ad21$DiffTypes != "CR"]) # subset
full.list$Ad21$age_ek <- annot_Ad21$age[match(colnames(full.list$Ad21),annot_Ad21$cells)]
full.list$Ad21$diff_ek <- annot_Ad21$DiffTypes[match(colnames(full.list$Ad21),annot_Ad21$cells)]
table(full.list$Ad21$age_ek,full.list$Ad21$diff_ek,useNA = "always")
full.list$Ad21$dataset <- "Ad21"

# My17
annot_My17 <- read_csv("data//M/DiffTypes_My17.csv") ; table(annot_My17$age,annot_My17$DiffTypes,useNA="always") ; dim(annot_My17)
full.list$My17 <- subset(full.list$My17,cells=annot_My17$cells[ ! annot_My17$DiffTypes %in% c("CR","SP","IPCearly") & annot_My17$age!="11" ] ) # subset
full.list$My17$age_ek <- annot_My17$age[match(colnames(full.list$My17),annot_My17$cells)]
full.list$My17$diff_ek <- annot_My17$DiffTypes[match(colnames(full.list$My17),annot_My17$cells)]
table(full.list$My17$age_ek,full.list$My17$diff_ek,useNA = "always") 
full.list$My17$dataset <- "My17"

# LMl21
annot_LMl21 <- read_csv("data//M/DiffTypes_LMl21.csv")  ; table(annot_LMl21$age,annot_LMl21$DiffTypes,useNA="always") ; dim(annot_LMl21)
full.list$LMl21 <- subset(full.list$LMl21,cells=annot_LMl21$cells[ ! annot_LMl21$DiffTypes %in% c("CR")]) # subset
full.list$LMl21$age_ek <- annot_LMl21$age[match(colnames(full.list$LMl21),annot_LMl21$cells)]
full.list$LMl21$diff_ek <- annot_LMl21$DiffTypes[match(colnames(full.list$LMl21),annot_LMl21$cells)]
table(full.list$LMl21$age_ek,full.list$LMl21$diff_ek,useNA = "always") # all good
full.list$LMl21$dataset <- "LMl21"

# Jt19
annot_Jt19 <- read_csv("data/M/DiffTypes_Jt19.csv")  ; table(annot_Jt19$age,annot_Jt19$DiffTypes,useNA="always") ; dim(annot_Jt19)
# all good, no need to subset
full.list$Jt19$age_ek <- annot_Jt19$age[match(colnames(full.list$Jt19),annot_Jt19$cells)] 
table(full.list$Jt19$age_ek,useNA = "always") # there are some NA cells, remove them
full.list$Jt19$age_ek <- ifelse(is.na(full.list$Jt19$age_ek),"toremove",as.character(full.list$Jt19$age_ek))
full.list$Jt19 <- subset(full.list$Jt19,age_ek!="toremove")
full.list$Jt19$diff_ek <- annot_Jt19$DiffTypes[match(colnames(full.list$Jt19),annot_Jt19$cells)]
table(full.list$Jt19$age_ek,full.list$Jt19$diff_ek,useNA = "always") 
full.list$Jt19$dataset <- "Jt19"

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #



# CELL SELECTION FOR INTEGRATION ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

merged <- merge(full.list$Ad21,list(full.list$My17,full.list$LMl21,full.list$Jt19))
table(merged$age_ek ,useNA = "always") # all good

# SELECT CELLS UP TO 6000 per age group no matter celltype neither dataset

# re-define age groups
merged$age_ek <- plyr::revalue(merged$age_ek,c("12"="12","13"="13","14"="14","15"="15","16"="16-17","17"="16-17"))

# select 6000 per age
merged_f <- cell.selector(merged,colnames(merged),merged$age_ek,6000)
table(merged_f$age_ek ,useNA = "always") # all good

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# Run pairwise Seurat integration ----
# finding anchors between all datasets
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

full.list <- SplitObject(merged_f,split.by =  "dataset")

# run integration

IntM <- dataset_integration(list_seuratobjs=full.list,nfeats=10000,FindIntDims=1:15,FindInt_kScore=20,IntKweight=20)

# run umap
DefaultAssay(IntM) <- "integrated" 
IntM <- IntM %>% ScaleData(vars.to.regress="nFeature_RNA") %>% RunPCA %>% RunUMAP(dims = 1:25,return.model = TRUE,min.dist = 0.8,n.neighbors = 100,local.connectivity = 6)
DimPlot(IntM,group.by = "dataset",pt.size = 0.5 ) 
DimPlot(IntM,group.by = "age_ek",pt.size = 0.5, split.by = "dataset") 
DimPlot(IntM,group.by = "diff_ek",pt.size = 0.5 , split.by = "dataset") 


# save integration
#saveRDS(IntM,"out/dataset_integration/IntM.rds")

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


