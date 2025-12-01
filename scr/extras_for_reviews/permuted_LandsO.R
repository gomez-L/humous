


# PERMUTED LANSCAPES ORG, COMPARISON WITH NONPERMUTED BY LANDSCAPE AVERAGE VS SINGLE LANDSCAPE THROUGH SUBSTRATION MAPS

source("humous_v3/lib/lib_misc.R")
source("humous_v3/lib/lib_lands.R")
library(Seurat) ; library(SeuratObject) ; library(dplyr) ; library(ggplot2) ; library(png) ; library(dplyr); library(purrr) ; library(parallel) ; library(pbmcapply) ; library(matrixStats) ; library(magrittr) ; library(entropy)

# load seurat object used for calculating humous_v4 landscapes: included gene selection, scaling, and normalization
LandsS_O <- readRDS("humous_v4/out/landsO/LandsS_O.rds")
LandsS_O$cellnames <- colnames(LandsS_O)

# define and apply permutation strategy
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
table(LandsS_O$y_diff,LandsS_O$y_age) # 450 per age diff group (113 is the 25%) 

# SAMPLING STRATEGY: 25% of cells on each subsample, 10 subsamples with replacement, likelyhood of 27 cells being repeated in any of the samples 

# generate subsamples - retrieve metadata with cellnames
permutedLands_O_dfs <- mclapply(1:10, function(i) {
  set.seed(1234 + i)
  # Step 1: Sample metadata (25% of full set)
  sampled_df <- LandsS_O@meta.data %>%
    group_by(y_age, y_diff) %>%
    slice_sample(n = 113, replace = TRUE) %>%
    ungroup() %>%
    mutate(sample_id = i)
}, mc.cores = 100)
#saveRDS(permutedLands_O_dfs,"humous_v4/out/Lands_permuted_O/permutedLands_O_dfs.rds")


# visuzalize the sets
cell_sets <- lapply(permutedLands_H_dfs, function(df) df$cellnames)
names(cell_sets) <- paste0("sample_", seq_along(cell_sets))
all_cells <- unique(unlist(cell_sets))
binary_matrix <- sapply(cell_sets, function(set) all_cells %in% set) ; rownames(binary_matrix) <- all_cells
upset_df <- as.data.frame(binary_matrix)


# generate grid with indexes for the cells in each subsample
permutedLands_O_grids <- mclapply(permutedLands_O_dfs, function(i) {
  gridlist <- knn_array_medres(i$ordi_age_norm, i$ordi_diff_norm, k = 100)
}, mc.cores = 100)
#saveRDS(permutedLands_O_grids,"humous_v4/out/Lands_permuted_O/permutedLands_O_grids.rds")


# generate expression matrices for each subsample
permutedLands_O_exprmats <- mclapply(permutedLands_O_dfs, function(i) {
  expr_mat <- LandsS_O@assays$RNA@scale.data[, i$cellnames]
}, mc.cores = 100)
#saveRDS(permutedLands_O_exprmats,"humous_v4/out/Lands_permuted_O/permutedLands_O_exprmats.rds")


# calculate knn_maps for each subsample, executing serially cause paralelization fails
permutedLands_O_rowmeans <- list()
for (i in seq_along(permutedLands_O_exprmats)){
  row_means <- knn_rowMeans(permutedLands_O_exprmats[[i]], permutedLands_O_grids[[i]])
  rownames(row_means) <- rownames(permutedLands_O_exprmats[[i]])
  permutedLands_O_rowmeans[[i]] <- row_means
}
#saveRDS(permutedLands_O_rowmeans,"humous_v4/out/Lands_permuted_O/permutedLands_O_rowmeans.rds")

# visualize one permuted landscape (1 out of 10 permutations of the same landscape)
permutedLands_O_rowmeans[[1]]["SOX2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2")

# average knn_maps across subsamples (element-wise)
mean_array_LandsO <- DelayedArray(0 * permutedLands_O_rowmeans[[1]])
for (arr in permutedLands_O_rowmeans) { mean_array_LandsO <- mean_array_LandsO + arr} # Sum all arrays
mean_array_LandsO <- mean_array_LandsO / length(permutedLands_O_rowmeans) # Divide by number of arrays
#saveRDS(mean_array_LandsO,"humous_v4/out/Lands_permuted_O/mean_array_LandsO.rds")

# visualize one averaged landscape
mean_array_LandsO["SOX2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# COMPARE AVERAGED LANDSCAPES WITH ORIGINAL LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# load original landscapes knn object
L_MR_O <- readRDS("~/humous/humous_v4/out/landsO/L_MR_O.rds")

# plot one gene for the averaged vs the original version of the landscapes
gridExtra::grid.arrange(mean_array_LandsO["SOX2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2"),
                        L_MR_O["SOX2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2"),ncol=2)

# APPROACH 1 (fast result) - matrix differences
#################
  # flatten the arrays (faster computation) to matrices with 62500 columns instead of 3D array of 13564 250 250
A1 <- matrix(aperm(L_MR_O, c(1,2,3)), nrow = dim(L_MR_O)[1]) 
A2 <- matrix(aperm(mean_array_LandsO, c(1,2,3)), nrow = dim(L_MR_O)[1])

diff_list_O <- pbmclapply(seq_len(28), function(i) { # number of chunks is ceiling(nrow(A1) / 500) = 28 - 28 chunks of 500 rows each
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
}, mc.cores = 20)
# diff_list_O is a list of length 28, where each element is a matrix of 500 × 62500
  # reshape it
diff_matrix_O <- do.call(rbind, diff_list_O) ; rownames(diff_matrix_O) <- rownames(L_MR_O)   # 13400 × 62500
diff_array_O <- array(diff_matrix_O, dim = dim(L_MR_O)) ; dimnames(diff_array_O) <- list(rownames(L_MR_O), NULL, NULL) # 13564   250   250, like the shape of L_MR_O
#saveRDS(diff_array_O,"humous_v4/out/Lands_permuted_O/diff_array_O.rds")


# df with stats to see if the diff arrays are mostly made of zeros or there is any permuted landscape that changed significantly respect to the original one
df_stats_O <- data.frame(
  slice = 1:dim(diff_array_O)[1],
  mean_diff = apply(diff_array_O, 1, mean, na.rm = TRUE), # shows if the difference is centered around zero
  sd_diff = apply(diff_array_O, 1, sd, na.rm = TRUE), # captures the spread of the differences
  mad_diff = apply(diff_array_O, 1, function(x) mean(abs(x - mean(x, na.rm = TRUE)), na.rm = TRUE)), # robust alternative to SD, not sensitive to outliers
  max_abs_diff = apply(diff_array_O, 1, function(x) max(abs(x), na.rm = TRUE)), # tells if the landscape pair has large spikes in differences
  prop_near_zero = apply(diff_array_O, 1, function(x) mean(abs(x) < 0.05, na.rm = TRUE)) # how much the landscape pair is essentially unchanged, if near zero, mostly noise
)
#saveRDS(df_stats_O,"humous_v4/out/Lands_permuted_O/df_stats_O.rds")

# visualize the results
# most important metrics
ggplot(df_stats_O) + geom_point(aes(mean_diff,mad_diff)) + theme_classic() +  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(0, 1))
# see the landscapes and sustration map for one high MAD gene ZNF324B
gridExtra::grid.arrange(mean_array_LandsO["ZNF324B", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("ZNF324B"),
                        L_MR_O["ZNF324B", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("ZNF324B"),ncol=2)
plot(raster::raster(diff_array_O["ZNF324B",,]),breaks=seq(-10, 10, length.out=10),col=rev(colorRampPalette(c("green", "white", "red"))(length(seq(-10, 10, length.out=10))-1)))
# see them for one with low MAD
gridExtra::grid.arrange(mean_array_LandsO["MEF2C", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("MEF2C"),
                        L_MR_O["MEF2C", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("MEF2C"),ncol=2)
plot(flip(raster::raster(diff_array_O["MEF2C",,]),direction="y"),breaks=seq(-10, 10, length.out=10),col=rev(colorRampPalette(c("green", "white", "red"))(length(seq(-10, 10, length.out=10))-1)))

#################














# APPROACH 2 - SUBSTRACTION MAPS

SMplotter <- function(gene,Lspecies1,Lspecies2,palette,colors,breaks_colors){
  substracted_raster <- raster::overlay( (as.array(Lspecies1[gene,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() ) ,
                                         ( as.array(Lspecies2[gene,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() ) ,fun=function(x,y){ (scale(x)-scale(y)) })
  par(mar = rep(0, 4))
  plot(raster::flip(substracted_raster),breaks=breaks_colors,col=colors,legend = FALSE, axes = FALSE, box=FALSE)
}

SMplotter(gene="SOX2",Lspecies1=L_MR_O,Lspecies2=mean_array_LandsO,breaks=seq(-10, 10, length.out=10),col=rev(colorRampPalette(c("green", "white", "red"))(length(seq(-10, 10, length.out=10))-1)))


trial_diff <- raster::overlay( (as.array(L_MR_O["SOX2",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() ) ,
                 ( as.array(mean_array_LandsO["NEUROD2",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() ) ,fun=function(x,y){ (scale(x)-scale(y)) })
mean(trial_diff@data@values)
trial_sim <- raster::overlay( (as.array(L_MR_O["SOX2",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() ) ,
                              ( as.array(mean_array_LandsO["SOX2",,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() ) ,fun=function(x,y){ (scale(x)-scale(y)) })
mean(trial_sim@data@values)

plot(raster::flip(trial),breaks=seq(-10, 10, length.out=10),col=rev(colorRampPalette(c("green", "white", "red"))(length(seq(-10, 10, length.out=10))-1)))

# APPROACH 2 -  Structural similarity (SSIM) — for image-like comparison

library(SSIMmap) # https://cran.r-project.org/web/packages/SSIMmap/vignettes/Introduction_to_SSIMmap.html

image1 <- L_MR_O["SOX2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster()


image1 <- raster::raster(as.array(L_MR_O["SOX2", , ]) %>% rescale(to = c(0, 1)))
image2 <- terra::rast(raster::raster(mean_array_LandsO["SOX2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster())) #%>%  base_ggraster() +  ggtitle("SOX2")
ssim_raster(image1,image2) 

trial <- OpenImageR::image_similarity(L_MR_O["SOX2", , ], mean_array_LandsO["SOX2", , ],method="SSIM")
  
compute_ssim <- function(mat1, mat2) { OpenImageR::SSIM(mat1, mat2)}

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #








