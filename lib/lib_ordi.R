

# LIBRARY FOR ORDINALS
################################################################################################################################################################

library(bmrm) ; library(Seurat) ; library(scales) ; library(parallel)

# CELL SELECTION FUNCTION
#####################################
cell.selector <- function(seuratobject,cellnames,grouping,n){set.seed(1234);
  df <- data.frame(cellnames=cellnames,celltypes=grouping) ; library(dplyr) ; df_f <- df %>% group_by(celltypes) %>% sample_n(size = n)
  S_f <- subset(seuratobject,cells=df_f$cellnames)}
#####################################


# COST MATRIX FUNCTION
#####################################
costM <- function(X,y,...) {
  C <- costMatrix(y)
  C <- C / tabulate(y)[col(C)] / tabulate(y)
  C <- C/sum(C)
}
#####################################


# TRAINING FUNCTION - FULL MODEL
#####################################
training_ordi_full <- function(datamatrix,target,cost,lambda_full,epsilon_full,maxiter){
  set.seed(1234)
  # train full model
  full <- nrbm(ordinalRegressionLoss(datamatrix,target,C=cost),LAMBDA=lambda_full,EPSILON_TOL=epsilon_full,MAX_ITER=maxiter)
}
#####################################



# CUSTOM RED-TRAINING AND PREDICTION
#####################################
# ordinal SVM reduced model training-testing CV function ----
custom_red_and_pred <- function(fullmodel,xtrain,target,xtest,ngenesselect,
                                lambda_red,epsilon_red,maxiter,nfolds){
  set.seed(1234)
  print("extract feat weights from full model and select top X depending on the dataset for testing")
  w <- t(attr(fullmodel,"gradient")) 
  # check which genes are common between w and xtest and, based on this, select top X model genes (most negative or positive weights)
  print(paste0("select top ",ngenesselect," model genes"))
  common <- intersect(colnames(xtest),rownames(w))
  w_f <- as.matrix(w[common,])
  
  # feat selection on weights
  w_f <- w_f[order(w_f),] ; w_f <- w_f[c(1:ngenesselect,(length(w_f)-(ngenesselect-1)):length(w_f))] ; print(length(w_f))
  # subset xtrain and xtest by feats selected
  xtrain_red <- xtrain[,names(w_f)] ; print(dim(xtrain_red))
  xtest_red <- xtest[,names(w_f)] ; print(dim(xtest_red))
  
  costM <- function(X,y,...) {
    C <- costMatrix(y)
    C <- C / tabulate(y)[col(C)] / tabulate(y)
    C <- C/sum(C)
  }
  cost=costM(X=xtrain_red,y=target) ; print(length(target)) ; print(dim(cost))
  
  # train reduced model without crossval (for future predictions)
  redmodel <- nrbm(ordinalRegressionLoss(xtrain_red,target,C=cost),LAMBDA = lambda_red,EPSILON_TOL =epsilon_red ,MAX_ITER = maxiter)
  
  # run reduced model with crossval and predict
  print("train reduced model with CV and predict")
  folds <- balanced.cv.fold(target,nfolds)
  pred <- simplify2array(mclapply(levels(folds),mc.cores=16,function(f) {
    w <- nrbm(ordinalRegressionLoss(xtrain_red[folds!=f,],target[folds!=f],C=cost),LAMBDA = lambda_red,EPSILON_TOL =epsilon_red ,MAX_ITER = maxiter)
    Y <- predict(w,xtest_red) ; Y[folds!=f] <- NA ; Y
  }))
  # prediction result aggregated across CVs
  pred <- rowSums(pred,na.rm=TRUE)
  # scale prediction from 0 to 1
  #pred <- rescale(pred,to=c(0,1))
  
  # return reduced model (without fold) and prediction vector
  return(list(redmodel=redmodel,pred=pred))
  
}


# ACTIVATION OF ORDINALS RESULTS
#####################################
ordi_activate <- function(ordiscore){
  
  # detect if not unimodal, if so, enter the loop
  if(diptest::dip.test(ordiscore,B=100)$p.value < 0.05){
    
      print("is not unimodal, check if is uniform")
      
      if( ks.test( ordiscore, "punif")$p.value > 0.05 ) {
        print("is uniform, nothing to do")
        return(ordiscore)
        
      }else{
        print("not uniform, smooth tanh*1")
        ordiscore <- tanh((ordiscore - median(ordiscore))*1)
      }
      
  }else{
    print("is unimodal, strong tahn")
    sd <- sd(ordiscore) ; print(sd)
    if( sd<0.15 ){
      print("sd<0.15, strong tanh*5")
      ordiscore <- tanh((ordiscore - median(ordiscore))*5)  
    }else if(sd<0.175){
      print("sd<0.175, strong tanh*2.5")
      ordiscore <- tanh((ordiscore - median(ordiscore))*2.5)
    }else{
      print("sd>0.175, strong tanh*2")
      ordiscore <- tanh((ordiscore - median(ordiscore))*2)
    }
  }
  
}
    

# NORMALIZATION OF ORDINALS RESULTS
#####################################
ordi_normalize <- function(ordiscore,ordigroups=NULL,ordigroups_other=NULL,limits1=NULL,limits2=NULL,limits3=NULL,applyWinsor=TRUE){
  
  ordiscore <- rescale(ordi_activate(ordiscore),to=c(0,1))
  
  if(is.null(limits1)){
    print("limits not set, letting default and applying tanh activation")
  }else{
    print("applying limits")
    for (i in levels(as.factor(ordigroups))){
      print(i)
      if(i==1){
        ordiscore[ordigroups==i] <- rescale(ordi_activate(ordiscore[ordigroups==i]),to=limits1)
      }else if(i==2){
        ordiscore[ordigroups==i] <- rescale(ordi_activate(ordiscore[ordigroups==i]),to=limits2)
      }else{
        ordiscore[ordigroups==i] <- rescale(ordi_activate(ordiscore[ordigroups==i]),to=limits3)
      }
    }
  }
  
  if(is.null(ordigroups_other)){
    print("2D correction not needed")
  }else{
    print("applying 2D correction")
    for (i in levels(as.factor(ordigroups_other))){ordiscore[ordigroups_other==i] <- rescale(ordi_activate(ordiscore[ordigroups_other==i]),to=c(0,1))}
  }

  if(applyWinsor==TRUE){
    print("applying winsor transformation")
    ordiscore <- Winsorize(ordiscore,probs=c(0.01,0.99),type=1)
  }else{
    print("not applying winsor transformation, letting default")
    return(ordiscore)
  }
  
  return(ordiscore)
}
    


