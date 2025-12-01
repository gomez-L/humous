

# ORDINAL MODELING OF AGE AND DIFF IN ORGHUMAN

source("lib/lib_ordi.R")
library(Seurat) ; library(dplyr); library(readr) ; library(ggplot2) ; library(plyr) ; library(dplyr) ; library(scales) ; library(purrr) ; library(DescTools)

# LOAD INTEGRATION RESULTS ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
IntHO <- readRDS("out/dataset_integration/IntHO.rds")
# remove human cells and rename object to IntO
IntO <- subset(IntHO,dataset %in% c("Av19","Ck19","Kp19","Pt20"))
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# ANNOTATE CELL TYPES DEFINED UPON INTEGRATION RESULTS 
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
DiffTypes_Org_Int <- read_csv("out/dataset_integration/DiffTypes_IntO.csv")
IntO$diff_ek <- DiffTypes_Org_Int$DiffTypes_final[match(colnames(IntO),DiffTypes_Org_Int$cells)]
IntO$age_ek <- DiffTypes_Org_Int$age_final[match(colnames(IntO),DiffTypes_Org_Int$cells)]
DimPlot(IntO,group.by = "diff_ek")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# CELL SELECTION FOR ORDI TRAINING, BALANCED BY AGE GROUPS AND DIFF GROUPS ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
table(IntO$age_ek,IntO$diff_ek,useNA = "always")

IntO$y_age <- as.integer(as.character(plyr::revalue(as.character(IntO$age_ek),c("1-2 m"="1","2-3 m"="2","3-5 m"="3","6-7 m"="4"))))
IntO$y_diff <- as.integer(as.character(plyr::revalue(as.character(IntO$diff_ek),c("RG"="1","IPC"="2","N"="3")))) ; table(IntO$y_age,IntO$y_diff,useNA = "always")

IntO$grouping_ordi <- paste0(IntO$y_age,"_",IntO$y_diff) ; table(IntO$grouping_ordi)  
IntO_train <- cell.selector(seuratobject=IntO,cellnames=colnames(IntO),grouping=IntO$grouping_ordi,n=50*4) 
table(IntO_train$grouping_ordi,useNA = "always") ; table(IntO_train$y_age,IntO_train$y_diff,useNA = "always") 
#saveRDS(IntO_train,"out/ordinals_age_diff/O/IntO_train.rds")

# annotate in IntO which cells were used for training ordis (those in IntO_train)
IntO$ordi_split <- ifelse(colnames(IntO) %in% colnames(IntO_train),"ordi_train","ordi_test")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# Load cells for landscapes - those for which we want to calculate ordinals
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
Org_4Lands <- read_csv("out/dataset_integration/SelectCells_IntO_forLandscapes.csv")
IntO_L <- subset(IntO,cells=Org_4Lands$x)
table(IntO_L$age_ek,IntO_L$diff_ek)
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# ORDINALS - Reconstruct age and differentiation
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# AGE 
#######################
fullAGE <- training_ordi_full(datamatrix=t(as.matrix(IntO_train@assays$integrated@scale.data)),
                              target=IntO_train$y_age,
                              cost=costM(X=t(as.matrix(IntO_train@assays$integrated@scale.data)),y=IntO_train$y_age),
                              lambda_full=0.1,epsilon_full=1e-7,maxiter=1000) 
#saveRDS(fullAGE,"out/ordinals_age_diff/O/model_fullAGE_O.rds") # save full model

# predict all data and evaluate performance
redAGE <-  custom_red_and_pred(fullmodel=fullAGE,xtrain=t(as.matrix(IntO_train@assays$integrated@scale.data)),target=IntO_train$y_age,
                               xtest=t(as.matrix(IntO_L@assays$integrated@scale.data)),ngenesselect=25,lambda_red=0.05,epsilon_red=1e-7,maxiter=1000,nfolds=20)
#saveRDS(redAGE$redmodel,"out/ordinals_age_diff/O/model_redAGE_O.rds") # save reduced model

# store prediction on df rescaled
predAGE <- data.frame( cellnames=c(names(redAGE$pred)), pred=c(scales::rescale(redAGE$pred,to=c(0,1))))
IntO_L$ordi_age <- predAGE$pred[match(colnames(IntO_L),predAGE$cellnames)]

# visualize prediction
ggplot(IntO_L@meta.data) + geom_density(aes(ordi_age,color=as.character(y_age))) + facet_wrap(~y_diff,ncol=1) + theme_bw() 
ggplot(IntO_L@meta.data) + geom_point(aes(ordi_age,y="y",color=as.character(y_age))) + facet_wrap(~grouping_ordi,ncol=1) + theme_bw() 
ggplot(IntO_L@meta.data) + geom_density(aes(ordi_age,color=as.character(y_age))) + theme_bw() 
ggplot(IntO_L@meta.data,aes(x=ordi_age,y=y_age,color=as.character(y_age)))+ geom_jitter()  + geom_boxplot() + facet_wrap(~y_diff,ncol=1)  + theme_bw() 

DimPlot(IntO_L, group.by =  "age_ek")
FeaturePlot(IntO_L,features = "ordi_age",cols=c("lightblue","blue","orange","red"))
#######################


# DIFF 
#######################
fullDIFF <- training_ordi_full(datamatrix=t(as.matrix(IntO_train@assays$integrated@scale.data)),
                               target=IntO_train$y_diff,
                               cost=costM(X=t(as.matrix(IntO_train@assays$integrated@scale.data)),y=IntO_train$y_diff),
                               lambda_full=0.1,epsilon_full=1e-7,maxiter=1000) 
#saveRDS(fullDIFF,"out/ordinals_age_diff/O/model_fullDIFF_O.rds") # save full model

# predict all data and evaluate performance
redDIFF <-  custom_red_and_pred(fullmodel=fullDIFF,xtrain=t(as.matrix(IntO_train@assays$integrated@scale.data)),target=IntO_train$y_diff,
                                xtest=t(as.matrix(IntO_L@assays$integrated@scale.data)),ngenesselect=25,lambda_red=0.1,epsilon_red=1e-7,maxiter=1000,nfolds=2)
#saveRDS(redDIFF$redmodel,"out/ordinals_age_diff/O/model_redDIFF_O.rds") # save reduced model

# store prediction on df rescaled
predDIFF <- data.frame( cellnames=c(names(redDIFF$pred)), pred=c(scales::rescale(redDIFF$pred,to=c(0,1))))
IntO_L$ordi_diff <- predDIFF$pred[match(colnames(IntO_L),predDIFF$cellnames)]

# visualize prediction
ggplot(IntO_L@meta.data) + geom_density(aes(ordi_diff,color=as.character(y_diff))) + facet_wrap(~y_age,ncol=1) + theme_bw() 
ggplot(IntO_L@meta.data) + geom_point(aes(ordi_diff,y="y",color=as.character(y_diff))) + facet_wrap(~grouping_ordi,ncol=1) + theme_bw() 
ggplot(IntO_L@meta.data) + geom_density(aes(ordi_diff,color=as.character(y_diff))) + theme_bw() 
ggplot(IntO_L@meta.data,aes(x=ordi_diff,y=y_diff,color=as.character(y_diff)))+ geom_jitter()  + geom_boxplot()  + theme_bw() + facet_wrap(~y_age,ncol=1)
#######################
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# NORMALIZE ORDIS
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
ggplot(IntO_L@meta.data) + geom_point(aes(ordi_age,ordi_diff,color=age_ek))
ggplot(IntO_L@meta.data) + geom_point(aes(ordi_age,ordi_diff,color=diff_ek))


IntO_L$ordi_diff_norm <- NA ; IntO_L$ordi_age_norm <- NA 
IntO_L$ordi_diff_norm <- ordi_normalize(ordiscore=IntO_L$ordi_diff,
                                        ordigroups=IntO_L$y_diff,
                                        ordigroups_other=IntO_L$y_age,
                                        limits1=c(0,0.24),limits2=c(0.21,0.5),limits3=c(0.4,1),
                                        applyWinsor=TRUE)

IntO_L$ordi_age_norm <- ordi_normalize(ordiscore=IntO_L$ordi_age,
                                       ordigroups=IntO_L$y_age,
                                       ordigroups_other=IntO_L$y_diff,
                                       applyWinsor=TRUE)


ggplot(IntO_L@meta.data) + geom_point(aes(ordi_age_norm,ordi_diff_norm,color=diff_ek))
ggplot(IntO_L@meta.data) + geom_point(aes(ordi_age_norm,ordi_diff_norm,color=age_ek))
ggplot(IntO_L@meta.data) + geom_point(aes(ordi_age_norm,ordi_diff_norm,color=dataset_proto2))

#saveRDS(IntO_L,"out/ordinals_age_diff/O/IntO_ordi.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #




