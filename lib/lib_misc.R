

# LIBRARY FOR MISCELANEOUS HUMOUS FUNCTIONS

library(Seurat) ; library(SeuratObject) ; library(bmrm) ; library(parallel) ; library(scales) ; library(diptest) ; library(fitdistrplus)


# gene selection for landscapes, expressed in all groups (datasets) an more than 10 cells overall
gene_selector <- function(seuratobj,grouping,threshold_grouping,threshold_ncells){
  
  df <- data.frame(cells=colnames(seuratobj),grouping=grouping)
  
  crossgroup_genes <- local({
    crossgroup_genes <- list()
    for(i in levels(as.factor(df$grouping))){
      
      sumsvect <- rowSums(seuratobj@assays$RNA@counts[,df$cells[df$grouping==i]])
      crossgroup_genes[[i]] <- names(sumsvect)[which(rowSums(seuratobj@assays$RNA@counts[,df$cells[df$grouping==i]]) > threshold_grouping)] 
    }
    crossgroup_genes <- Reduce(intersect,crossgroup_genes) 
  })
  
  sumsvect <- rowSums(seuratobj@assays$RNA@counts)
  expr_genes10 <- names(sumsvect)[which(rowSums(seuratobj@assays$RNA@counts) > threshold_ncells)]
  
  selected_genes <- intersect(crossgroup_genes,expr_genes10)
  
  return(selected_genes)
  
}


# dataset integration parwise (if reference_dataset is not provided) or with reference ----
dataset_integration <- function(list_seuratobjs,nfeats,FindIntDims,FindInt_kScore,IntKweight,reference_dataset=NULL){
  
  if(is.null(reference_dataset)){
    
    for (i in 1:length(list_seuratobjs)) {
      list_seuratobjs[[i]] <- NormalizeData(list_seuratobjs[[i]], verbose = FALSE)
      list_seuratobjs[[i]] <- FindVariableFeatures(list_seuratobjs[[i]], selection.method = "vst",nfeatures=nfeats, verbose = FALSE)
    }
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = list_seuratobjs)
    mer.anchors <- FindIntegrationAnchors(object.list = list_seuratobjs, anchor.features=features,dims=FindIntDims,k.score = FindInt_kScore)
    mer.integrated <- IntegrateData(anchorset = mer.anchors, features = features,k.weight = IntKweight)
    return(mer.integrated)
    
  }else{
    
  for (i in 1:length(list_seuratobjs)) {
    list_seuratobjs[[i]] <- NormalizeData(list_seuratobjs[[i]], verbose = FALSE)
    list_seuratobjs[[i]] <- FindVariableFeatures(list_seuratobjs[[i]], selection.method = "vst",nfeatures=nfeats, verbose = FALSE)
  }
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = list_seuratobjs)
  reference_data <- which(names(full.list) == reference_dataset)
  mer.anchors <- FindIntegrationAnchors(object.list = list_seuratobjs, reference=reference_data,anchor.features=features,dims=FindIntDims,k.score = FindInt_kScore)
  mer.integrated.ref <- IntegrateData(anchorset = mer.anchors, features = features,k.weight = IntKweight)
  return(mer.integrated.ref)
  }
}


# integrate datasets into an already integrated object ----
double_dataset_integration <- function(list_seuratobjs,nfeats,FindIntDims,FindInt_kScore,IntKweight,reference_dataset,features){
  reference_data <- which(names(full.list) == reference_dataset)
  mer.anchors <- FindIntegrationAnchors(object.list = list_seuratobjs, reference=reference_data,anchor.features=features,dims=FindIntDims,k.score = FindInt_kScore)
  mer.integrated.ref <- IntegrateData(anchorset = mer.anchors, features = features,k.weight = IntKweight)
  return(mer.integrated.ref)
}



# cell selection depending on grouping ----
cell.selector <- function(seuratobject,cellnames,grouping,n){set.seed(1234);
  df <- data.frame(cellnames=cellnames,celltypes=grouping) ; library(dplyr) ; df_f <- df %>% group_by(celltypes) %>% sample_n(size = n)
  S_f <- subset(seuratobject,cells=df_f$cellnames)}



# cost matrix function ----
costM <- function(X,y,...) {
  C <- costMatrix(y)
  C <- C / tabulate(y)[col(C)] / tabulate(y)
  C <- C/sum(C)
}














