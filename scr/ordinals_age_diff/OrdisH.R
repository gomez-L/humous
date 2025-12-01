

# ORDINAL MODELING OF AGE AND DIFF IN HUMAN

source("lib/lib_ordi.R")
library(Seurat) ; library(dplyr); library(readr) ; library(ggplot2) ; library(plyr) ; library(dplyr) ; library(scales) ; library(purrr) ; library(DescTools)

# LOAD INTEGRATION RESULTS ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
IntH <- readRDS("out/dataset_integration/IntH.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# ANNOTATE CELL TYPES DEFINED UPON INTEGRATION RESULTS
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
DiffTypes_Hm_Int <- read_csv("out/IntH/DiffTypes_IntH.csv")
IntH$diff_ek <- DiffTypes_Hm_Int$DiffTypes_final[match(colnames(IntH),DiffTypes_Hm_Int$cells)]
DimPlot(IntH,group.by = "diff_ek")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# CELL SELECTION FOR ORDI TRAINING, BALANCED BY AGE GROUPS AND DIFF GROUPS ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# define integer label on age and diff levels 
IntH$y_age <- as.integer(as.character(plyr::revalue(as.character(IntH$age_ek),c("12"="1","15-16"="2","17-18"="3","20-21"="4","23-24"="5"))))
IntH$y_diff <- as.integer(as.character(plyr::revalue(as.character(IntH$diff_ek),c("RG"="1","IPC"="2","N"="3")))) ; table(IntH$y_age,IntH$y_diff)

# cell selection for training ordinals
IntH$grouping_ordi <- paste0(IntH$y_age,"_",IntH$y_diff) ; table(IntH$grouping_ordi)  
IntH_train <- cell.selector(seuratobject=IntH,cellnames=colnames(IntH),grouping=IntH$grouping_ordi,n=50*5) ; table(IntH_train$grouping_ordi) ; table(IntH_train$y_age,IntH_train$y_diff) 
#saveRDS(IntH_train,"out/ordinals_age_diff/H/IntH_train.rds")

IntH$ordi_split <- ifelse(colnames(IntH) %in% colnames(IntH_train),"ordi_train","ordi_test")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# Load cells selected for landscapes - those for which we want to calculate ordinals
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
Hm_4Lands <- read_csv("out/dataset_integration/SelectCells_IntH_forLandscapes.csv")
IntH_L <- subset(IntH,cells=Hm_4Lands$x)
table(IntH_L$age_ek,IntH_L$diff_ek)
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# ORDINALS - Reconstruct age and differentiation
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# AGE 
#######################
fullAGE <- training_ordi_full(datamatrix=t(as.matrix(IntH_train@assays$integrated@scale.data)),
                              target=IntH_train$y_age,
                              cost=costM(X=t(as.matrix(IntH_train@assays$integrated@scale.data)),y=IntH_train$y_age),
                              lambda_full=0.05,epsilon_full=1e-7,maxiter=1000) 
#saveRDS(fullAGE,"out/ordinals_age_diff/H/model_fullAGE_H.rds") # save full model

# predict all data and evaluate performance
redAGE <-  custom_red_and_pred(fullmodel=fullAGE,xtrain=t(as.matrix(IntH_train@assays$integrated@scale.data)),target=IntH_train$y_age,
                               xtest=t(as.matrix(IntH_L@assays$integrated@scale.data)),ngenesselect=25,lambda_red=0.01,epsilon_red=1e-7,maxiter=1000,nfolds=20)
#saveRDS(redAGE$redmodel,"out/ordinals_age_diff/H/model_redAGE_H.rds") # save reduced model

# store prediction on df rescaled
predAGE <- data.frame( cellnames=c(names(redAGE$pred)), pred=c(scales::rescale(redAGE$pred,to=c(0,1))))
IntH_L$ordi_age <- predAGE$pred[match(colnames(IntH_L),predAGE$cellnames)]


# visualize prediction
ggplot(IntH_L@meta.data) + geom_density(aes(ordi_age,color=as.character(y_age))) + facet_wrap(~y_diff,ncol=1) + theme_bw() 
ggplot(IntH_L@meta.data) + geom_point(aes(ordi_age,y="y",color=as.character(y_age))) + facet_wrap(~grouping_ordi,ncol=1) + theme_bw() 
ggplot(IntH_L@meta.data) + geom_density(aes(ordi_age,color=as.character(y_age))) + theme_bw() 
ggplot(IntH_L@meta.data,aes(x=ordi_age,y=y_age,color=as.character(y_age)))+ geom_jitter()  + geom_boxplot() + facet_wrap(~y_diff,ncol=1)  + theme_bw() 
#######################


# DIFF 
#######################
fullDIFF <- training_ordi_full(datamatrix=t(as.matrix(IntH_train@assays$integrated@scale.data)),
                               target=IntH_train$y_diff,
                               cost=costM(X=t(as.matrix(IntH_train@assays$integrated@scale.data)),y=IntH_train$y_diff),
                               lambda_full=0.1,epsilon_full=1e-7,maxiter=1000) 
#saveRDS(fullDIFF,"out/ordinals_age_diff/H/model_fullDIFF_H.rds") # save full model

# predict all data and evaluate performance
redDIFF <-  custom_red_and_pred(fullmodel=fullDIFF,xtrain=t(as.matrix(IntH_train@assays$integrated@scale.data)),target=IntH_train$y_diff,
                                xtest=t(as.matrix(IntH_L@assays$integrated@scale.data)),ngenesselect=25,lambda_red=0.1,epsilon_red=1e-7,maxiter=1000,nfolds=2)
#saveRDS(redDIFF$redmodel,"out/ordinals_age_diff/H/model_redDIFF_H.rds") # save reduced model

# store prediction on df rescaled
predDIFF <- data.frame( cellnames=c(names(redDIFF$pred)), pred=c(scales::rescale(redDIFF$pred,to=c(0,1))))
IntH_L$ordi_diff <- predDIFF$pred[match(colnames(IntH_L),predDIFF$cellnames)]

# visualize prediction
ggplot(IntH_L@meta.data) + geom_density(aes(ordi_diff,color=as.character(y_diff))) + facet_wrap(~y_age,ncol=1) + theme_bw() 
ggplot(IntH_L@meta.data) + geom_point(aes(ordi_diff,y="y",color=as.character(y_diff))) + facet_wrap(~grouping_ordi,ncol=1) + theme_bw() 
ggplot(IntH_L@meta.data) + geom_density(aes(ordi_diff,color=as.character(y_diff))) + theme_bw() 
ggplot(IntH_L@meta.data,aes(x=ordi_diff,y=y_diff,color=as.character(y_diff)))+ geom_jitter()  + geom_boxplot()  + theme_bw() 
#######################
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# NORMALIZE ORDIS
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
ggplot(IntH_L@meta.data) + geom_point(aes(ordi_age,ordi_diff,color=age_ek))
ggplot(IntH_L@meta.data) + geom_point(aes(ordi_age,ordi_diff,color=diff_ek))

IntH_L$ordi_diff_norm <- NA ; IntH_L$ordi_age_norm <- NA 
IntH_L$ordi_diff_norm <- ordi_normalize(ordiscore=IntH_L$ordi_diff,
                                        ordigroups=IntH_L$y_diff,
                                        ordigroups_other=IntH_L$y_age,
                                        limits1=c(0,0.24),limits2=c(0.21,0.47),limits3=c(0.4,1),
                                        applyWinsor=TRUE)

IntH_L$ordi_age_norm <- ordi_normalize(ordiscore=IntH_L$ordi_age,
                                       ordigroups=IntH_L$y_age,
                                       ordigroups_other=IntH_L$y_diff,
                                       applyWinsor=TRUE)


ggplot(IntH_L@meta.data) + geom_point(aes(ordi_age_norm,ordi_diff_norm,color=diff_ek))
ggplot(IntH_L@meta.data) + geom_point(aes(ordi_age_norm,ordi_diff_norm,color=age_ek))


#saveRDS(IntH_L,"out/ordinals_age_diff/H/IntH_ordi.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #



