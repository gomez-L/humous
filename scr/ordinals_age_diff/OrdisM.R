
# ORDINAL MODELING OF AGE AND DIFF IN MOUSE

source("lib/lib_ordi.R")
library(Seurat) ; library(dplyr); library(readr) ; library(ggplot2) ; library(plyr) ; library(dplyr) ; library(scales) ; library(purrr) ; library(DescTools)

# LOAD INTEGRATION RESULTS ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
IntM <- readRDS("out/dataset_integration/IntM.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# ANNOTATE CELL TYPES DEFINED UPON INTEGRATION RESULTS
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
DiffTypes_Ms_Int <- read_csv("out/dataset_integration/DiffTypes_IntM.csv")
IntM$diff_ek <- DiffTypes_Ms_Int$DiffTypes_final[match(colnames(IntM),DiffTypes_Ms_Int$cells)]
DimPlot(IntM,group.by = "diff_ek")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# CELL SELECTION FOR ORDI TRAINING, BALANCED BY AGE GROUPS AND DIFF GROUPS ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

IntM$y_age <- as.integer(as.character(plyr::revalue(as.character(IntM$age_ek),c("12"="1","13"="2","14"="3","15"="4","16-17"="5"))))
IntM$y_diff <- as.integer(as.character(plyr::revalue(as.character(IntM$diff_ek),c("RG"="1","IPC"="2","N"="3")))) ; table(IntM$y_age,IntM$y_diff)

IntM$grouping_ordi <- paste0(IntM$y_age,"_",IntM$y_diff) ; table(IntM$grouping_ordi)  
IntM_train <- cell.selector(seuratobject=IntM,cellnames=colnames(IntM),grouping=IntM$grouping_ordi,n=50*5) ; table(IntM_train$grouping_ordi) ; table(IntM_train$y_age,IntM_train$y_diff) 
#saveRDS(IntM_train,"out/ordinals_age_diff/M/IntM_train.rds")

# annotate in IntM which cells were used for training ordis (those in IntM_train)
IntM$ordi_split <- ifelse(colnames(IntM) %in% colnames(IntM_train),"ordi_train","ordi_test")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# Load cells for landscapes - those for which we want to calculate ordinals
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
Ms_4Lands <- read_csv("out/dataset_integration/SelectCells_IntM_forLandscapes.csv")
IntM_L <- subset(IntM,cells=Ms_4Lands$x)
table(IntM_L$age_ek,IntM_L$diff_ek)
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# ORDINALS - Reconstruct age and differentiation
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# AGE 
#######################
fullAGE <- training_ordi_full(datamatrix=t(as.matrix(IntM_train@assays$integrated@scale.data)),
                              target=IntM_train$y_age,
                              cost=costM(X=t(as.matrix(IntM_train@assays$integrated@scale.data)),y=IntM_train$y_age),
                              lambda_full=0.05,epsilon_full=1e-7,maxiter=1000) 
#saveRDS(fullAGE,"out/ordinals_age_diff/M/model_fullAGE_M.rds") # save full model

# predict all data and evaluate performance
redAGE <-  custom_red_and_pred(fullmodel=fullAGE,xtrain=t(as.matrix(IntM_train@assays$integrated@scale.data)),target=IntM_train$y_age,
                               xtest=t(as.matrix(IntM_L@assays$integrated@scale.data)),ngenesselect=25,lambda_red=0.05,epsilon_red=1e-7,maxiter=1000,nfolds=20)
#saveRDS(redAGE$redmodel,"out/ordinals_age_diff/M/model_redAGE_M.rds") # save reduced model

# store prediction on df rescaled
predAGE <- data.frame( cellnames=c(names(redAGE$pred)), pred=c(scales::rescale(redAGE$pred,to=c(0,1))))
IntM_L$ordi_age <- predAGE$pred[match(colnames(IntM_L),predAGE$cellnames)]

# visualize prediction
ggplot(IntM_L@meta.data) + geom_density(aes(ordi_age,color=as.character(y_age))) + facet_wrap(~y_diff,ncol=1) + theme_bw() ; #ggsave("humous_v3/pdf/IntM/ordiAGE_density_bydiff.pdf",useDingbats=FALSE)
ggplot(IntM_L@meta.data) + geom_point(aes(ordi_age,y="y",color=as.character(y_age))) + facet_wrap(~grouping_ordi,ncol=1) + theme_bw() ; #ggsave("humous_v3/pdf/IntM/ordiAGE_scatter_byagediff.pdf",useDingbats=FALSE)
ggplot(IntM_L@meta.data) + geom_density(aes(ordi_age,color=as.character(y_age))) + theme_bw() ; #ggsave("humous_v3/pdf/IntM/ordiAGE_density_byage.pdf",useDingbats=FALSE)
ggplot(IntM_L@meta.data,aes(x=ordi_age,y=y_age,color=as.character(y_age)))+ geom_jitter()  + geom_boxplot() + facet_wrap(~y_diff,ncol=1)  + theme_bw() ; #ggsave("humous_v3/pdf/IntM/ordiAGE_jitterboxplot_byage.pdf",useDingbats=FALSE)
#######################


# DIFF 
#######################
fullDIFF <- training_ordi_full(datamatrix=t(as.matrix(IntM_train@assays$integrated@scale.data)),
                               target=IntM_train$y_diff,
                               cost=costM(X=t(as.matrix(IntM_train@assays$integrated@scale.data)),y=IntM_train$y_diff),
                               lambda_full=0.5,epsilon_full=1e-7,maxiter=1000) 
#saveRDS(fullDIFF,"out/ordinals_age_diff/M/model_fullDIFF_M.rds") # save full model

# predict all data and evaluate performance
redDIFF <-  custom_red_and_pred(fullmodel=fullDIFF,xtrain=t(as.matrix(IntM_train@assays$integrated@scale.data)),target=IntM_train$y_diff,
                                xtest=t(as.matrix(IntM_L@assays$integrated@scale.data)),ngenesselect=25,lambda_red=0.5,epsilon_red=1e-7,maxiter=1000,nfolds=2)
#saveRDS(redDIFF$redmodel,"out/ordinals_age_diff/M/model_redDIFF_M.rds") # save reduced model

# store prediction on df rescaled
predDIFF <- data.frame( cellnames=c(names(redDIFF$pred)), pred=c(scales::rescale(redDIFF$pred,to=c(0,1))))
IntM_L$ordi_diff <- predDIFF$pred[match(colnames(IntM_L),predDIFF$cellnames)]

# visualize prediction
ggplot(IntM_L@meta.data) + geom_density(aes(ordi_diff,color=as.character(y_diff))) + facet_wrap(~y_age,ncol=1) + theme_bw() ; #ggsave("humous_v3/pdf/IntM/ordiDIFF_density_byage.pdf",useDingbats=FALSE)
ggplot(IntM_L@meta.data) + geom_point(aes(ordi_diff,y="y",color=as.character(y_diff))) + facet_wrap(~grouping_ordi,ncol=1) + theme_bw() ; #ggsave("humous_v3/pdf/IntM/ordiDIFF_scatter_byagediff.pdf",useDingbats=FALSE)
ggplot(IntM_L@meta.data) + geom_density(aes(ordi_diff,color=as.character(y_diff))) + theme_bw() ; #ggsave("humous_v3/pdf/IntM/ordiDIFF_density_bydiff.pdf",useDingbats=FALSE)
ggplot(IntM_L@meta.data,aes(x=ordi_diff,y=y_diff,color=as.character(y_diff)))+ geom_jitter()  + geom_boxplot()  + theme_bw() ; #ggsave("humous_v3/pdf/IntM/ordiDIFF_jitterboxplot_bydiff.pdf",useDingbats=FALSE)
#######################
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# NORMALIZE ORDIS
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

ggplot(IntM_L@meta.data) + geom_point(aes(ordi_age,ordi_diff,color=diff_ek))

IntM_L$ordi_diff_norm <- NA ; IntM_L$ordi_age_norm <- NA 
IntM_L$ordi_diff_norm <- ordi_normalize(ordiscore=IntM_L$ordi_diff,
                                        ordigroups=IntM_L$y_diff,
                                        ordigroups_other=IntM_L$y_age,
                                        limits1=c(0,0.26),limits2=c(0.21,0.5),limits3=c(0.4,1),
                                        applyWinsor=TRUE)

IntM_L$ordi_age_norm <- ordi_normalize(ordiscore=IntM_L$ordi_age,
                                       ordigroups=IntM_L$y_age,
                                       ordigroups_other=IntM_L$y_diff,
                                       applyWinsor=TRUE)


ggplot(IntM_L@meta.data) + geom_point(aes(ordi_age_norm,ordi_diff_norm,color=diff_ek))

#saveRDS(IntM_L,"out/ordinals_age_diff/M/IntM_ordi.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #




