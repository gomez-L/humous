

# PERMUTED LANSCAPES MOUSE, COMPARISON WITH NONPERMUTED BY LANDSCAPE AVERAGE VS SINGLE LANDSCAPE THROUGH SUBSTRATION MAPS

source("humous_v3/lib/lib_misc.R")
source("humous_v3/lib/lib_lands.R")
library(Seurat) ; library(SeuratObject) ; library(dplyr) ; library(ggplot2) ; library(png) ; library(dplyr); library(purrr) ; library(parallel) ; library(pbmcapply) ; library(matrixStats) ; library(magrittr) ; library(entropy)

# load seurat object used for calculating humous_v4 landscapes: included gene selection, scaling, and normalization
LandsS_M <- readRDS("humous_v4/out/landsM/LandsS_M.rds")
LandsS_M$cellnames <- colnames(LandsS_M)

# define and apply permutation strategy
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
table(LandsS_M$y_diff,LandsS_M$y_age) # 485 per age diff group (122 is the 25%) 

###############

# SAMPLING STRATEGY
# 25% of cells on each subsample, 10 subsamples with replacement, likelyhood of 27 cells being repeated in any of the samples 

# generate subsamples - retrieve metadata with cellnames
permutedLands_M_dfs <- mclapply(1:10, function(i) {
  set.seed(1234 + i)
  # Step 1: Sample metadata (25% of full set)
  sampled_df <- LandsS_M@meta.data %>%
    group_by(y_age, y_diff) %>%
    slice_sample(n = 122, replace = TRUE) %>%
    ungroup() %>%
    mutate(sample_id = i)
}, mc.cores = 100)
#saveRDS(permutedLands_M_dfs,"humous_v4/out/Lands_permuted_M/permutedLands_M_dfs.rds")


# visuzalize the sets
cell_sets <- lapply(permutedLands_H_dfs, function(df) df$cellnames)
names(cell_sets) <- paste0("sample_", seq_along(cell_sets))
all_cells <- unique(unlist(cell_sets))
binary_matrix <- sapply(cell_sets, function(set) all_cells %in% set) ; rownames(binary_matrix) <- all_cells
upset_df <- as.data.frame(binary_matrix)


# generate grid with indexes for the cells in each subsample
permutedLands_M_grids <- mclapply(permutedLands_M_dfs, function(i) {
  gridlist <- knn_array_medres(i$ordi_age_norm, i$ordi_diff_norm, k = 100)
}, mc.cores = 100)
#saveRDS(permutedLands_M_grids,"humous_v4/out/Lands_permuted_M/permutedLands_M_grids.rds")


# generate expression matrices for each subsample
permutedLands_M_exprmats <- mclapply(permutedLands_M_dfs, function(i) {
  expr_mat <- LandsS_M@assays$RNA@scale.data[, i$cellnames]
}, mc.cores = 100)
#saveRDS(permutedLands_M_exprmats,"humous_v4/out/Lands_permuted_M/permutedLands_M_exprmats.rds")


# calculate knn_maps for each subsample, executing serially cause paralelization fails
permutedLands_M_rowmeans <- list()
for (i in seq_along(permutedLands_M_exprmats)){
  row_means <- knn_rowMeans(permutedLands_M_exprmats[[i]], permutedLands_M_grids[[i]])
  rownames(row_means) <- rownames(permutedLands_M_exprmats[[i]])
  permutedLands_M_rowmeans[[i]] <- row_means
}
#saveRDS(permutedLands_M_rowmeans,"humous_v4/out/Lands_permuted_M/permutedLands_M_rowmeans.rds")
permutedLands_M_rowmeans <- readRDS("~/humous/humous_v4/out/Lands_permuted_M/permutedLands_M_rowmeans.rds")

# visualize one permuted landscape (1 out of 10 permutations of the same landscape)
permutedLands_M_rowmeans[[1]]["SOX2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2")

# average knn_maps across subsamples (element-wise)
mean_array_LandsM <- DelayedArray(0 * permutedLands_M_rowmeans[[1]])
for (arr in permutedLands_M_rowmeans) { mean_array_LandsM <- mean_array_LandsM + arr} # Sum all arrays
mean_array_LandsM <- mean_array_LandsM / length(permutedLands_M_rowmeans) # Divide by number of arrays

# visualize one averaged landscape
mean_array_LandsM["Sox2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("Sox2")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# COMPARE AVERAGED LANDSCAPES WITH ORIGINAL LANDSCAPES
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# load original landscapes knn object
L_MR_M <- readRDS("~/humous/humous_v4/out/landsM/L_MR_M.rds")

# plot one gene for the averaged vs the original version of the landscapes
gridExtra::grid.arrange(mean_array_LandsM["Sox2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("Sox2"),
                        L_MR_M["Sox2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("Sox2"),ncol=2)

# APPROACH 1 (fast result) - matrix differences
#################
# flatten the arrays (faster computation) to matrices with 62500 columns instead of 3D array of 13564 250 250
A1 <- matrix(aperm(L_MR_M, c(1,2,3)), nrow = dim(L_MR_M)[1]) 
A2 <- matrix(aperm(mean_array_LandsM, c(1,2,3)), nrow = dim(L_MR_M)[1])

diff_list_M <- pbmclapply(seq_len(26), function(i) { # number of chunks is ceiling(nrow(A1) / 500) = 26 - 26 chunks of 500 rows each
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
# diff_list_M is a list of length 26, where each element is a matrix of 500 × 62500
# reshape it
diff_matrix_M <- do.call(rbind, diff_list_M) ; rownames(diff_matrix_M) <- rownames(L_MR_M)   # 13400 × 62500
diff_array_M <- array(diff_matrix_M, dim = dim(L_MR_M)) ; dimnames(diff_array_M) <- list(rownames(L_MR_M), NULL, NULL) # 13564   250   250, like the shape of L_MR_M
#saveRDS(diff_array_M,"humous_v4/out/Lands_permuted_M/diff_array_M.rds")


# df with stats to see if the diff arrays are mostly made of zeros or there is any permuted landscape that changed significantly respect to the original one
df_stats_M <- data.frame(
  slice = 1:dim(diff_array_M)[1],
  mean_diff = apply(diff_array_M, 1, mean, na.rm = TRUE), # shows if the difference is centered around zero
  sd_diff = apply(diff_array_M, 1, sd, na.rm = TRUE), # captures the spread of the differences
  mad_diff = apply(diff_array_M, 1, function(x) mean(abs(x - mean(x, na.rm = TRUE)), na.rm = TRUE)), # robust alternative to SD, not sensitive to outliers
  max_abs_diff = aply(diff_array_M, 1, function(x) max(abs(x), na.rm = TRUE)), # tells if the landscape pair has large spikes in differences
  prop_near_zero = apply(diff_array_M, 1, function(x) mean(abs(x) < 0.05, na.rm = TRUE)) # how much the landscape pair is essentially unchanged, if near zero, mostly noise
)
#saveRDS(df_stats_M,"humous_v4/out/Lands_permuted_M/df_stats_M.rds")

# visualize the results
# most important metrics
ggplot(df_stats_M) + geom_point(aes(mean_diff,mad_diff)) + theme_classic() +  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(0, 1))
# see the landscapes and sustration map for one high MAD gene Txlng
gridExtra::grid.arrange(mean_array_LandsM["Txlng", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("Txlng"),
                        L_MR_M["Txlng", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("Txlng"),ncol=2)
plot(raster::raster(diff_array_M["Txlng",,]),breaks=seq(-10, 10, length.out=10),col=rev(colorRampPalette(c("green", "white", "red"))(length(seq(-10, 10, length.out=10))-1)))
# see them for one with low MAD
gridExtra::grid.arrange(mean_array_LandsM["Neurod2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("Neurod2"),
                        L_MR_M["Neurod2", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("Neurod2"),ncol=2)
plot(raster::raster(diff_array_M["Neurod2",,]),breaks=seq(-10, 10, length.out=10),col=rev(colorRampPalette(c("green", "white", "red"))(length(seq(-10, 10, length.out=10))-1)))

#################



