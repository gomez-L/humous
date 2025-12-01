

# RAW PERMUTED HUMOUS - subsample 1 - ORDINAL MODELING OF AGE AND DIFF IN HUMAN

source("humous_v3/lib/lib_misc.R")
library(Seurat) ; library(dplyr); library(readr) ; library(ggplot2) ; library(plyr) ; library(dplyr) ; library(scales) ; library(purrr) ; library(DescTools)

# load integration results, annotate cells, and select subset for permutation 1
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# load integrated object - with cell selection only per age
IntH <- readRDS("~/humous/humous_v3/out/IntH/IntH.rds")

# select cells permutation 1 and re-annotate age and diff groups
DiffTypes_Hm_Int_annotated <- read_csv("humous_v4/out/raw_permuted_humous/cellselection_RPH/DiffTypes_Hm_Int_annotated.csv")
IntH$age_ek <- DiffTypes_Hm_Int_annotated$age[match(colnames(IntH),DiffTypes_Hm_Int_annotated$cells)]
IntH$diff_ek <- DiffTypes_Hm_Int_annotated$DiffTypes_final[match(colnames(IntH),DiffTypes_Hm_Int_annotated$cells)]

IntH_pe1 <- subset(IntH,cells=DiffTypes_Hm_Int_annotated$cells[DiffTypes_Hm_Int_annotated$sel1==TRUE])
table(IntH_pe1$age_ek,IntH_pe1$diff_ek)
table(IntH_pe1$diff_ek,IntH_pe1$age_ek,IntH_pe1$dataset)
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# CELL SELECTION FOR ORDI TRAINING, BALANCED BY AGE GROUPS (X5) AND DIFF GROUPS (X3) ----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# define integer label on age and diff levels 
IntH_pe1$y_age <- as.integer(as.character(plyr::revalue(as.character(IntH_pe1$age_ek),c("12"="1","15-16"="2","17-18"="3","20-21"="4","23-24"="5"))))
IntH_pe1$y_diff <- as.integer(as.character(plyr::revalue(as.character(IntH_pe1$diff_ek),c("RG"="1","IPC"="2","N"="3")))) 
table(IntH_pe1$y_age,IntH_pe1$y_diff)

# cell selection for training ordinals (balanced by diff and age, not by dataset because this is regressed out in the integration)
# Paste y_age and y_diff for selecting a balanced n of cells per group
IntH_pe1$grouping_ordi <- paste0(IntH_pe1$y_age,"_",IntH_pe1$y_diff) ; table(IntH_pe1$grouping_ordi)  
IntH_pe1_train <- cell.selector(seuratobject=IntH_pe1,cellnames=colnames(IntH_pe1),grouping=IntH_pe1$grouping_ordi,n=50*5)
table(IntH_pe1_train$grouping_ordi) ; table(IntH_pe1_train$y_age,IntH_pe1_train$y_diff) 
#saveRDS(IntH_pe1_train,"humous_v4/out/raw_permuted_humous/RPH_pe1/ordiH_pe1/IntH_pe1_train.rds")

# annotate in IntH_pe1 which cells were used for training ordis (those in IntH_pe1_train)
IntH_pe1$ordi_split <- ifelse(colnames(IntH_pe1) %in% colnames(IntH_pe1_train),"ordi_train","ordi_test")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #



# ORDINALS - Reconstruct age and differentiation
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# AGE 
#######################
fullAGE_pe1 <- training_ordi_full(datamatrix=t(as.matrix(IntH_pe1_train@assays$integrated@scale.data)),
                              target=IntH_pe1_train$y_age,
                              cost=costM(X=t(as.matrix(IntH_pe1_train@assays$integrated@scale.data)),y=IntH_pe1_train$y_age),
                              lambda_full=0.05,epsilon_full=1e-7,maxiter=1000) 
#saveRDS(fullAGE_pe1,"humous_v4/out/raw_permuted_humous/RPH_pe1/ordiH_pe1/model_fullAGE_pe1_H.rds") # save full model

# predict all data and evaluate performance
redAGE_pe1 <-  custom_red_and_pred(fullmodel=fullAGE_pe1,xtrain=t(as.matrix(IntH_pe1_train@assays$integrated@scale.data)),target=IntH_pe1_train$y_age,
                               xtest=t(as.matrix(IntH_pe1@assays$integrated@scale.data)),ngenesselect=25,lambda_red=0.01,epsilon_red=1e-7,maxiter=1000,nfolds=20)
#saveRDS(redAGE_pe1$redmodel,"humous_v4/out/raw_permuted_humous/RPH_pe1/ordiH_pe1/model_redAGE_H.rds") # save reduced model

# store prediction on df rescaled
predAGE_pe1 <- data.frame( cellnames=c(names(redAGE_pe1$pred)), pred=c(scales::rescale(redAGE_pe1$pred,to=c(0,1))))
IntH_pe1$ordi_age <- predAGE_pe1$pred[match(colnames(IntH_pe1),predAGE_pe1$cellnames)]

# visualize prediction
ggplot(IntH_pe1@meta.data) + geom_density(aes(ordi_age,color=as.character(y_age))) + facet_wrap(~y_diff,ncol=1) + theme_bw() 
ggplot(IntH_pe1@meta.data) + geom_point(aes(ordi_age,y="y",color=as.character(y_age))) + facet_wrap(~grouping_ordi,ncol=1) + theme_bw() 
ggplot(IntH_pe1@meta.data) + geom_density(aes(ordi_age,color=as.character(y_age))) + theme_bw() 
ggplot(IntH_pe1@meta.data,aes(x=ordi_age,y=y_age,color=as.character(y_age)))+ geom_jitter()  + geom_boxplot() + facet_wrap(~y_diff,ncol=1)  + theme_bw() 
#######################


# DIFF 
#######################
fullDIFF_pe1 <- training_ordi_full(datamatrix=t(as.matrix(IntH_pe1_train@assays$integrated@scale.data)),
                               target=IntH_pe1_train$y_diff,
                               cost=costM(X=t(as.matrix(IntH_pe1_train@assays$integrated@scale.data)),y=IntH_pe1_train$y_diff),
                               lambda_full=0.1,epsilon_full=1e-7,maxiter=1000) 
#saveRDS(fullDIFF_pe1,"humous_v4/out/raw_permuted_humous/RPH_pe1/ordiH_pe1/model_fullDIFF_pe1_H.rds") # save full model

# predict all data and evaluate performance
redDIFF_pe1 <-  custom_red_and_pred(fullmodel=fullDIFF_pe1,xtrain=t(as.matrix(IntH_pe1_train@assays$integrated@scale.data)),target=IntH_pe1_train$y_diff,
                                xtest=t(as.matrix(IntH_pe1@assays$integrated@scale.data)),ngenesselect=25,lambda_red=0.1,epsilon_red=1e-7,maxiter=1000,nfolds=2)
#saveRDS(redDIFF_pe1$redmodel,"humous_v4/out/raw_permuted_humous/RPH_pe1/ordiH_pe1/model_redDIFF_pe1_H.rds") # save reduced model

# store prediction on df rescaled
predDIFF_pe1 <- data.frame( cellnames=c(names(redDIFF_pe1$pred)), pred=c(scales::rescale(redDIFF_pe1$pred,to=c(0,1))))
IntH_pe1$ordi_diff <- predDIFF_pe1$pred[match(colnames(IntH_pe1),predDIFF_pe1$cellnames)]

# visualize prediction
ggplot(IntH_pe1@meta.data) + geom_density(aes(ordi_diff,color=as.character(y_diff))) + facet_wrap(~y_age,ncol=1) + theme_bw() 
ggplot(IntH_pe1@meta.data) + geom_point(aes(ordi_diff,y="y",color=as.character(y_diff))) + facet_wrap(~grouping_ordi,ncol=1) + theme_bw() 
ggplot(IntH_pe1@meta.data) + geom_density(aes(ordi_diff,color=as.character(y_diff))) + theme_bw() 
ggplot(IntH_pe1@meta.data,aes(x=ordi_diff,y=y_diff,color=as.character(y_diff)))+ geom_jitter()  + geom_boxplot()  + theme_bw() 
#######################
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# NORMALIZE ORDIS
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
ggplot(IntH_pe1@meta.data) + geom_point(aes(ordi_age,ordi_diff,color=age_ek))
ggplot(IntH_pe1@meta.data) + geom_point(aes(ordi_age,ordi_diff,color=diff_ek))

IntH_pe1$ordi_diff_norm <- NA ; IntH_pe1$ordi_age_norm <- NA 
IntH_pe1$ordi_diff_norm <- ordi_normalize(ordiscore=IntH_pe1$ordi_diff,
                                        ordigroups=IntH_pe1$y_diff,
                                        ordigroups_other=IntH_pe1$y_age,
                                        limits1=c(0,0.24),limits2=c(0.21,0.47),limits3=c(0.4,1),
                                        applyWinsor=TRUE)

IntH_pe1$ordi_age_norm <- ordi_normalize(ordiscore=IntH_pe1$ordi_age,
                                       ordigroups=IntH_pe1$y_age,
                                       ordigroups_other=IntH_pe1$y_diff,
                                       applyWinsor=TRUE)


ggplot(IntH_pe1@meta.data) + geom_point(aes(ordi_age_norm,ordi_diff_norm,color=diff_ek))
ggplot(IntH_pe1@meta.data) + geom_point(aes(ordi_age_norm,ordi_diff_norm,color=age_ek))


#saveRDS(IntH_pe1,"humous_v4/out/raw_permuted_humous/RPH_pe1/ordiH_pe1/IntH_pe1_ordi.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #





