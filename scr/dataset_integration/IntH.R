

# INTEGRATION OF HUMAN DATASETS

source("lib/lib_misc.R")
library(Seurat) ; library(dplyr); library(readr) ; library(ggplot2) ; library(plyr) ; library(dplyr) ; library(scales)


# load raw mouse datasets and join them in a list ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
full.list <- list(Kb21=readRDS("data/H/Kb21.rds"),
                  Kn17=readRDS("data/H/Kn17.rds"),
                  Pt21=readRDS("data/H/Pt21.rds"))
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# filter cells on each dataset to keep ages and types of interest, annotate also by dataset (for integration)
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# load annotations prepared by esther, filter and attach age and diff vectors

# Kb21
annot_Kb21 <- read_csv("data/H/DiffTypes_Kb21.csv") ; dim(annot_Kb21) ; table(annot_Kb21$age,annot_Kb21$DiffTypes,annot_Kb21$region_final,useNA="always")
annot_Kb21$DiffTypes <- ifelse(annot_Kb21$DiffTypes=="Neuron","N",as.character(annot_Kb21$DiffTypes))
# no cells to be subseted now, all in age and diff groups of interest
full.list$Kb21$age_ek <- annot_Kb21$age[match(colnames(full.list$Kb21),annot_Kb21$cells)]
full.list$Kb21$diff_ek <- annot_Kb21$DiffTypes[match(colnames(full.list$Kb21),annot_Kb21$cells)]
# add region information 
full.list$Kb21$region_ek <- annot_Kb21$region_final[match(colnames(full.list$Kb21),annot_Kb21$cells)]
table(full.list$Kb21$age_ek,full.list$Kb21$diff_ek,full.list$Kb21$region_ek,useNA = "always") # no nas, so ok
# keep only parietal cells 
full.list$Kb21 <- subset(full.list$Kb21,region_ek=="parietal") # subset
table(full.list$Kb21$age_ek,full.list$Kb21$diff_ek,useNA = "always") # no nas, so ok
full.list$Kb21$dataset <- "Kb21"

# Kn17
annot_Kn17 <- read_csv("data/H/DiffTypes_Kn17.csv") ; dim(annot_Kn17) ; table(annot_Kn17$age,annot_Kn17$DiffTypes,annot_Kn17$region_final,useNA="always")
# some difftypes are not the ones we need, tag them
annot_Kn17$DiffTypes <- ifelse(is.na(annot_Kn17$DiffTypes),"other",as.character(annot_Kn17$DiffTypes))
# N is misspeled as Neuron, correct it
annot_Kn17$DiffTypes <- ifelse(annot_Kn17$DiffTypes=="Neuron","N",as.character(annot_Kn17$DiffTypes))
table(annot_Kn17$age,annot_Kn17$DiffTypes,annot_Kn17$region_final)
annot_Kn17 <- annot_Kn17[annot_Kn17$DiffTypes %in% c("N","IPC","RG") & annot_Kn17$age %in% c("12","16","20") & annot_Kn17$region_final=="parietal",] # subset types, ages and region

full.list$Kn17$age_ek <- annot_Kn17$age[match(colnames(full.list$Kn17),annot_Kn17$cells)]
full.list$Kn17$diff_ek <- annot_Kn17$DiffTypes[match(colnames(full.list$Kn17),annot_Kn17$cells)]
# add region information
full.list$Kn17$region_ek <- annot_Kn17$region_final[match(colnames(full.list$Kn17),annot_Kn17$cells)]
full.list$Kn17 <- subset(full.list$Kn17,cells=annot_Kn17$cells) # subset
table(full.list$Kn17$age_ek,full.list$Kn17$diff_ek,full.list$Kn17$region_ek,useNA = "always")
full.list$Kn17$dataset <- "Kn17"

# Pt21
annot_Pt21 <- read_csv("data/H/DiffTypes_Pt21.csv") ; dim(annot_Pt21) ; table(annot_Pt21$age,annot_Pt21$DiffTypes,useNA="always") 
full.list$Pt21$age_ek <- annot_Pt21$age[match(colnames(full.list$Pt21),annot_Pt21$cells)]
full.list$Pt21$diff_ek <- annot_Pt21$DiffTypes[match(colnames(full.list$Pt21),annot_Pt21$cells)]
table(full.list$Pt21$age_ek,full.list$Pt21$diff_ek,useNA = "always")
# nothing to subset, and area information is not provided anywhere
full.list$Pt21$dataset <- "Pt21"
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #



# CELL SELECTION FOR INTEGRATION ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# as integration will balance by dataset/protocol, and cell types will be redefined based on integration, select based on age

merged <- merge(full.list$Kb21,list(full.list$Kn17,full.list$Pt21))
table(merged$age_ek,merged$region_ek,useNA = "always") # all good


# re-define age groups
merged$age_ek <- plyr::revalue(as.character(merged$age_ek),
                               c("12"="12","15"="15-16","16"="15-16","17"="17-18","18"="17-18","20"="20-21","21"="20-21","23"="23-24","24"="23-24"))
table(merged$age_ek ,useNA = "always") # all good

# select cells
mergedH_f <- local({
  mergedH_f1 <- subset(merged,  age_ek %in% c("15-16","17-18","20-21","23-24")  )
  mergedH_f1 <- cell.selector(mergedH_f1,colnames(mergedH_f1),mergedH_f1$age_ek,6000)
  
  mergedH_f2 <- subset(merged,age_ek %in% c("12"))

  mergedH_f <- merge(mergedH_f1,mergedH_f2)
})

table(mergedH_f$age_ek ,useNA = "always") # all good

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #



# Run pairwise Seurat integration ----
# finding anchors between all datasets
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

full.list <- SplitObject(mergedH_f,split.by =  "dataset")

# run integration

IntH <- pairwise_integration(list_seuratobjs=full.list,nfeats=10000,FindIntDims=1:15,FindInt_kScore=20,IntKweight=20)

DefaultAssay(IntH) <- "integrated" 
IntH <- IntH %>% ScaleData(vars.to.regress="nFeature_RNA") %>% RunPCA %>% RunUMAP(dims = 1:15,return.model = TRUE,min.dist = 0.8,n.neighbors = 100,local.connectivity = 6)

IntH$diff_ek <- ifelse(IntH$diff_ek=="Neuron","N",as.character(IntH$diff_ek))

DimPlot(IntH,group.by = "dataset",pt.size = 0.5 ) 
DimPlot(IntH,group.by = "age_ek",pt.size = 0.5, split.by = "dataset") 
DimPlot(IntH,group.by = "diff_ek",pt.size = 0.5, split.by = "dataset") 


# save integration
#saveRDS(IntH,"out/dataset_integration/IntH.rds")

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #







