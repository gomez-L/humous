

# PERMUTED LANSCAPES HUMAN, COMPARISON WITH NONPERMUTED BY LANDSCAPE AVERAGE VS SINGLE LANDSCAPE THROUGH SUBSTRATION MAPS

source("humous_v3/lib/lib_misc.R")
source("humous_v3/lib/lib_lands.R")
library(Seurat) ; library(SeuratObject) ; library(dplyr) ; library(ggplot2) ; library(png) ; library(dplyr); library(purrr) ; library(parallel) ; library(pbmcapply) ; library(matrixStats) ; library(magrittr) ; library(entropy)

# load seurat object used for calculating humous_v4 landscapes: included gene selection, scaling, and normalization
LandsS_H <- readRDS("humous_v4/out/landsH/LandsS_H.rds")
LandsS_H$cellnames <- colnames(LandsS_H)

# define and apply permutation strategy
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
table(LandsS_H$y_diff,LandsS_H$y_age) # 490 per age diff group except for one, 320

###############

# SAMPLING STRATEGY
  # 25% of cells on each subsample, 10 subsamples with replacement, likelyhood of 27 cells being repeated in any of the samples 

# generate subsamples - retrieve metadata with cellnames
permutedLands_H_dfs <- mclapply(1:10, function(i) {
  set.seed(1234 + i)
  # Step 1: Sample metadata (25% of full set)
  sampled_df <- LandsS_H@meta.data %>%
    group_by(y_age, y_diff) %>%
    slice_sample(n = 122, replace = TRUE) %>%
    ungroup() %>%
    mutate(sample_id = i)
}, mc.cores = 100)
#saveRDS(permutedLands_H_dfs,"humous_v4/out/Lands_permuted_H/permutedLands_H_dfs.rds")


# visuzalize the sets
cell_sets <- lapply(permutedLands_H_dfs, function(df) df$cellnames)
names(cell_sets) <- paste0("sample_", seq_along(cell_sets))
all_cells <- unique(unlist(cell_sets))
binary_matrix <- sapply(cell_sets, function(set) all_cells %in% set) ; rownames(binary_matrix) <- all_cells
upset_df <- as.data.frame(binary_matrix)


# generate grid with indexes for the cells in each subsample
permutedLands_H_grids <- mclapply(permutedLands_H_dfs, function(i) {
  gridlist <- knn_array_medres(i$ordi_age_norm, i$ordi_diff_norm, k = 100)
}, mc.cores = 100)
#saveRDS(permutedLands_H_grids,"humous_v4/out/Lands_permuted_H/permutedLands_H_grids.rds")


# generate expression matrices for each subsample
permutedLands_H_exprmats <- mclapply(permutedLands_H_dfs, function(i) {
  expr_mat <- LandsS_H@assays$RNA@scale.data[, i$cellnames]
}, mc.cores = 100)
#saveRDS(permutedLands_H_exprmats,"humous_v4/out/Lands_permuted_H/permutedLands_H_exprmats.rds")


# calculate knn_maps for each subsample, executing serially cause paralelization fails
permutedLands_H_rowmeans <- list()
for (i in seq_along(permutedLands_H_exprmats)){
  row_means <- knn_rowMeans(permutedLands_H_exprmats[[i]], permutedLands_H_grids[[i]])
  rownames(row_means) <- rownames(permutedLands_H_exprmats[[i]])
  permutedLands_H_rowmeans[[i]] <- row_means
}
#saveRDS(permutedLands_H_rowmeans,"humous_v4/out/Lands_permuted_H/permutedLands_H_rowmeans.rds")

# visualize one permuted landscape (1 out of 10 permutations of the same landscape)
permutedLands_H_rowmeans[[1]]["SOX2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2")

# average knn_maps across subsamples (element-wise)
mean_array_LandsH <- DelayedArray(0 * permutedLands_H_rowmeans[[1]])
for (arr in permutedLands_H_rowmeans) { mean_array_LandsH <- mean_array_LandsH + arr} # Sum all arrays
mean_array_LandsH <- mean_array_LandsH / length(permutedLands_H_rowmeans) # Divide by number of arrays

# visualize one averaged landscape
mean_array_LandsH["SOX2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# COMPARE AVERAGED LANDSCAPES WITH ORIGINAL LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# load original landscapes knn object
L_MR_H <- readRDS("~/humous/humous_v4/out/landsH/L_MR_H.rds")

# plot one gene for the averaged vs the original version of the landscapes
gridExtra::grid.arrange(mean_array_LandsH["SOX2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2"),
                        L_MR_H["SOX2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2"),ncol=2)

# APPROACH 1 (fast result) - matrix differences
#################
# flatten the arrays (faster computation) to matrices with 62500 columns instead of 3D array of 13564 250 250
A1 <- matrix(aperm(L_MR_H, c(1,2,3)), nrow = dim(L_MR_H)[1]) 
A2 <- matrix(aperm(mean_array_LandsH, c(1,2,3)), nrow = dim(L_MR_H)[1])

n_chunks <- as.integer(ceiling(nrow(A1)/500))
diff_list_H <- pbmclapply(seq_len(n_chunks), function(i) { # number of chunks is ceiling(nrow(A1) / 500) = 26 - 26 chunks of 500 rows each
  message("Processing chunk ", i) # to print to console
  idx <- ((i - 1) * 500 + 1):min(i * 500, nrow(A1))
  list(A1[idx, , drop = FALSE], A2[idx, , drop = FALSE]) %>%
    lapply(function(mat) {
      mat %>%
        { (.-rowMins(.)) / (rowMaxs(.) - rowMins(.)) } %>%
        { replace(., is.nan(.), 0) } %>%
        { (.-rowMeans(.)) / rowSds(.) }
    }) %>%
    { .[[1]] - .[[2]] }
}, mc.cores = 10)
# diff_list_H is a list of length 28, where each element is a matrix of 500 × 62500
# reshape it
diff_matrix_H <- do.call(rbind, diff_list_H) ; rownames(diff_matrix_H) <- rownames(L_MR_H)   # 13400 × 62500
diff_array_H <- array(diff_matrix_H, dim = dim(L_MR_H)) ; dimnames(diff_array_H) <- list(rownames(L_MR_H), NULL, NULL) # 13564   250   250, like the shape of L_MR_H
saveRDS(diff_array_H,"humous_v4/out/Lands_permuted_H/diff_array_H.rds")

# diff_array_H_numeric <- diff_array_H
# dimnames(diff_array_H_numeric) <- list(NULL, NULL, NULL) ; rownames(diff_array_H_numeric) <- NULL
# diff_array_H_numeric <- array(as.numeric(diff_array_H_numeric), dim = dim(diff_array_H_numeric))
# dimnames(diff_array_H_numeric) <- list(rownames(diff_array_H), NULL, NULL)

# df with stats to see if the diff arrays are mostly made of zeros or there is any permuted landscape that changed significantly respect to the original one
df_stats_H <- data.frame(
  slice = 1:dim(diff_array_H)[1],
  mean_diff = apply(diff_array_H, 1, mean, na.rm = TRUE), # shows if the difference is centered around zero
  sd_diff = apply(diff_array_H, 1, sd, na.rm = TRUE), # captures the spread of the differences
  mad_diff = apply(diff_array_H, 1, function(x) mean(abs(x - mean(x, na.rm = TRUE)), na.rm = TRUE)), # robust alternative to SD, not sensitive to outliers
  max_abs_diff = apply(diff_array_H, 1, function(x) max(abs(x), na.rm = TRUE)), # tells if the landscape pair has large spikes in differences
  prop_near_zero = apply(diff_array_H, 1, function(x) mean(abs(x) < 0.05, na.rm = TRUE)) # how much the landscape pair is essentially unchanged, if near zero, mostly noise
)
saveRDS(df_stats_H,"humous_v4/out/Lands_permuted_H/df_stats_H.rds")

# visualize the results
# most important metrics
ggplot(df_stats_H) + geom_point(aes(mean_diff,mad_diff)) + theme_classic() +  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(0, 1))
# see the landscapes and sustration map for one high MAD gene ZNF324B
gridExtra::grid.arrange(mean_array_LandsH["GRWD1", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("GRWD1"),
                        L_MR_H["GRWD1", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("GRWD1"),ncol=2)
plot(raster::raster(diff_array_H_numeric["GRWD1",,]),breaks=seq(-10, 10, length.out=10),col=rev(colorRampPalette(c("green", "white", "red"))(length(seq(-10, 10, length.out=10))-1)))
# see them for one with low MAD
gridExtra::grid.arrange(mean_array_LandsO["MEF2C", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("MEF2C"),
                        L_MR_H["MEF2C", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("MEF2C"),ncol=2)
plot(raster::raster(diff_array_H["MEF2C",,]),breaks=seq(-10, 10, length.out=10),col=rev(colorRampPalette(c("green", "white", "red"))(length(seq(-10, 10, length.out=10))-1)))

#################

































# OTHER STUFF


# generate substraction maps for each landscape on each subsample and the
SMplotter <- function(gene,Lspecies1,Lspecies2,palette,colors,breaks_colors){
  substracted_raster <- raster::overlay( (as.array(Lspecies1[gene,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() ) ,
                                         ( as.array(Lspecies2[gene,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() ) ,fun=function(x,y){ (scale(x)-scale(y)) })
  # the following line is not needed! we don't want to rescale the substractions cause else all the genes look super different between species!
  #scaled_substracted_raster <- t(matrix(as.array(substracted_raster@data@values) %>% scales::rescale(to=c(0,1)), nrow = 250, ncol = 250))  %>% raster::raster()
  par(mar = rep(0, 4))
  plot(flip(substracted_raster),breaks=breaks_colors,col=colors,legend = FALSE, axes = FALSE, box=FALSE)
}

# try the function SMplotter on one gene
SMplotter(gene="SOX2",Lspecies1=L_MR_H,Lspecies2=L_MR_M,breaks=seq(-10, 10, length.out=10),col=rev(colorRampPalette(c("green", "white", "red"))(length(seq(-10, 10, length.out=10))-1)))








# generate landscapes for each subsample (deeplearning version, no axis)
mclapply(mc.cores=60, rownames(permutedLands_H_rowmeans[[i]]),function(i){
  png(paste0("humous_v4/out/pngLs/ml_pngLs/H/",i,"_H.png"),width=224,height=224) 
  print(as.array(permutedLands_H_rowmeans[[i]][i,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% ML_ggraster())
  dev.off()
})



row_means_1 <- knn_rowMeans(permutedLands_H_exprmats[[1]], permutedLands_H_grids[[1]])
rownames(row_means_1) <- rownames(permutedLands_H_exprmats[[1]])



permutedLands_H_rowmeans <- mclapply(seq_along(permutedLands_H_exprmats), function(i) {
  row_means <- knn_rowMeans(permutedLands_H_exprmats[[i]], permutedLands_H_grids[[i]])
  rownames(row_means) <- rownames(expr_mat)
  return(row_means)
}, mc.cores = 80)


permutedLands_H_rowmeans <- mclapply(seq_along(permutedLands_H_exprmats), function(i) {
  tryCatch({
    row_means <- knn_rowMeans(permutedLands_H_exprmats[[i]], permutedLands_H_grids[[i]])
    rownames(row_means) <- rownames(permutedLands_H_exprmats[[i]])
    return(row_means)
  }, error = function(e) {
    message(sprintf("Error on element %d: %s", i, e$message))
    return(NULL)
  })
}, mc.cores = 64)


row_means_1["EOMES", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2")







sampled_knn_maps_H <- mclapply(1:10, function(i) {
  
  set.seed(1234)
  # Step 1: Sample metadata (25% of full set)
  sampled_df <- LandsS_H@meta.data %>%
    group_by(y_age, y_diff) %>%
    slice_sample(n = 122, replace = TRUE) %>%
    ungroup() %>%
    mutate(sample_id = i)
  
  # Step 2: Compute knn grid (adjusted to 25% of full set)
  gridlist <- knn_array_medres(sampled_df$ordi_age_norm, sampled_df$ordi_diff_norm, k = 100)
  
  # Step 3: Subset expression matrix for selected genes and sampled cells
  expr_mat <- LandsS_H@assays$RNA@scale.data[, sampled_df$cellnames]
  
  # Step 4: Compute rowMeans over the gridlist
  row_means <- knn_rowMeans(expr_mat, gridlist)
  rownames(row_means) <- rownames(expr_mat)
  
  # Output: list with the rowMeans and maybe metadata
  list(
    sample_id = i,
    knn_means = row_means,
    gridlist=gridlist,
    sampled_df = sampled_df
  )
  
}, mc.cores = 100)

saveRDS(sampled_knn_maps_H,"humous_v4/out/Lands_permuted_H/sampled_knn_maps_grids_H.rds")

# check some genes in one subsample
plots <- map(c("SOX2", "EOMES", "NEUROD6", "JUNB", "APBB2", "MEF2C"), ~
               sampled_knn_maps[[1]]$row_means[.x, , ] %>%
               as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle(.x))

# check just one gene in one subsample
sampled_knn_maps_grids_H[[6]]$knn_means["SOX2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2")















# calculate subsamples
set.seed(1234)

samplingH <- map(1:10,function(i){
  cell.selector(seuratobject=LandsS_H,cellnames=colnames(LandsS_H),grouping=LandsS_H$y_age_y_diff,n=122)
})
scaled_matrices <- lapply(samplingH, function(seurat_obj) {seurat_obj@assays$RNA@scale.data})
avg_matrix <- Reduce("+", scaled_matrices) / length(scaled_matrices)

LandsS_H_f <- subset(LandsS_H,cells=colnames(avg_matrix))
gridlistH_sampling <- knn_array_medres(LandsS_H_f$ordi_age_norm,LandsS_H_f$ordi_diff_norm,k=100L)

m <- avg_matrix[c("SOX2","EOMES","NEUROD6","JUNB","APBB2","MEF2C"),]
knn_rowMeans <- function(t, i) {
  M <- m[, as.vector(i)]
  M <- array(M, c(nrow(m), dim(i)))
  DelayedArray::DelayedArray(rowMeans(M, dims = length(dim(i))))
}
L_MR_H <- knn_rowMeans(m,gridlistH_sampling)  ; rownames(L_MR_H) <- rownames(m) 



sampling_risky <- map(1:10,function(i){
  LandsS_H@meta.data %>%
    group_by(y_age,y_diff) %>%
    slice_sample(n=122,replace=FALSE) %>%
    ungroup() %>%
    mutate(sample_id=i)
})

gridlistH_samplingrisky <- mclapply(sampling_risky,
                                    function(df){
                                      knn_array_medres(df$ordi_age_norm,df$ordi_diff_norm,k=100) # adapted, before was 100, by rule of three, keeping 25%
                                    },mc.cores=80)

range(as.vector(gridlistH_samplingrisky[[1]]))
m <- LandsS_H@assays$RNA@scale.data[c("SOX2","EOMES","NEUROD6","JUNB","APBB2","MEF2C"),sampling_risky[[1]]$cellnames] # first for a subsample of genes
L_MR_H <- knn_rowMeans(m,gridlistH_samplingrisky[[1]])  ; rownames(L_MR_H) <- rownames(m) 











L_MR_H <- local({
  m <- avg_matrix[c("SOX2","EOMES","NEUROD6","JUNB","APBB2","MEF2C"),]
  landscapes_MR_H <- knn_rowMeans(m,gridlistH_sampling) ; rownames(landscapes_MR_H) <- rownames(m) ; landscapes_MR_H
})

knn_rowMeans <- function(m,i) {
  B <- local({
    batches <- split(seq(nrow(m)),seq(nrow(m)) %/% 100L)
    B <- lapply(batches,function(b){
      print(paste("batch",b,"out of", (length(batches)*100) )) # progress report
      M <- m[b,as.vector(i)]
      M <- array(M,c(length(b),dim(i)))
      DelayedArray(rowMeans(M,dims=length(dim(i))))
    })
    do.call(rbind,B)
  })       
}





sampling_risky <- map(1:10,function(i){
  LandsS_H@meta.data %>%
    group_by(y_age,y_diff) %>%
    slice_sample(n=320,replace=FALSE) %>%
    ungroup() %>%
    mutate(sample_id=i)
})
# check it went as expected
table(sampling_risky[[1]]$y_age,sampling_risky[[1]]$y_diff)

# calculate the gridlist for each subsample, at medium resolution
detectCores()
gridlistH_samplingrisky <- mclapply(sampling_risky,
                                    function(df){
                                      knn_array_medres(df$ordi_age_norm,df$ordi_diff_norm,k=100) # adapted, before was 100, by rule of three, keeping 25%
                                    },mc.cores=80)

# Compute the element-wise average
m <- LandsS_H@assays$RNA@scale.data[c("SOX2","EOMES","NEUROD6","JUNB","APBB2","MEF2C"),] # first for a subsample of genes

avg_matrix <- Reduce("+", mat_list) / length(mat_list)


L_MR_H <- local({
  m <- LandsS_H@assays$RNA@scale.data[c("SOX2","EOMES","NEUROD6","JUNB","APBB2","MEF2C"),]
  landscapes_MR_H <- knn_rowMeans(m,gridlistH_samplingrisky[[6]]) ; rownames(landscapes_MR_H) <- rownames(m) ; landscapes_MR_H
})
gridExtra::grid.arrange(as.array(L_MR_H["SOX2",,]) %>% #scales::rescale(to=c(0,1)) %>% 
                          raster::raster() %>% base_ggraster() + ggtitle("SOX2"),
                        as.array(L_MR_H["EOMES",,]) %>% #scales::rescale(to=c(0,1)) %>% 
                          raster::raster() %>% base_ggraster()+ ggtitle("EOMES") ,
                        as.array(L_MR_H["NEUROD6",,]) %>% #scales::rescale(to=c(0,1)) %>% 
                          raster::raster() %>% base_ggraster()+ ggtitle("NEUROD6"), 
                        as.array(L_MR_H["JUNB",,]) %>% #scales::rescale(to=c(0,1)) %>% 
                          raster::raster() %>% base_ggraster() + ggtitle("JUNB"),
                        as.array(L_MR_H["APBB2",,]) %>% #scales::rescale(to=c(0,1)) %>% 
                          raster::raster() %>% base_ggraster()+ ggtitle("APBB2") ,
                        as.array(L_MR_H["MEF2C",,]) %>% #scales::rescale(to=c(0,1)) %>% 
                          raster::raster() %>% base_ggraster()+ ggtitle("MEF2C"),ncol=3 )



##########


# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
