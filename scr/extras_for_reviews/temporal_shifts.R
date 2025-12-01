

# ASSESSING TEMPORAL SHIFTS IN TIME-SPECIFIC GENES BETWEEN SPECIES USING SUSTRACTION MAPS
  # SMs -> VGG16 embs -> clustering
    # if random noise -> biology // if consistent bias --> missalignment

source("humous_v3/lib/lib_misc.R")
source("humous_v3/lib/lib_lands.R")
library(Seurat) ; library(SeuratObject) ; library(dplyr) ; library(ggplot2) ; library(png) ; library(dplyr); library(purrr) ; library(parallel) ;
library(pbmcapply) ; library(matrixStats) ; library(magrittr) ; library(entropy) ; library(readr) ; library(raster)


# LOAD DATA
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# landscapes arrays
L_MR_H <- readRDS("~/humous/humous_v4/out/landsH/L_MR_H.rds")
L_MR_M <- readRDS("~/humous/humous_v4/out/landsM/L_MR_M.rds")
L_MR_O <- readRDS("~/humous/humous_v4/out/landsO/L_MR_O.rds")
# df with annotations on age and diff specificity (vgg16 clustering)
lands_clust <- read_csv("humous_v4/out/vgg16/GeneInfo.csv")
# indexing table
indexing_table_all <- read_csv("humous_v4/out/indexing_table_all.csv")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# to start - find which genes have age specificity in human and see if they have consistent temporal shifts in the other species
  # if consistent temporal shifts - there is a bias in data (for example, all human early are org late)
  # if inconsistent temporal shifts - there is no bias, just random noise due to biological variations (some shifts in one direction, others in another direction of age)
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# find all human genes that have specificity in age
table(lands_clust$cond,lands_clust$AnnotAge)
H_agespec <- lands_clust[lands_clust$cond=="H" & lands_clust$AnnotAge %in% c("e","l"),] # 2136 genes

# find those 2136 genes in M and O (not all of them have the same name)
#O_H_agespec <- lands_clust[lands_clust$gene %in% H_agespec$gene & lands_clust$cond=="O",] ; dim(O_H_agespec) # 2013
common_Hagespec_O <- indexing_table_all[indexing_table_all$hsapiens_H %in% H_agespec$gene & indexing_table_all$O==TRUE , ] # 2018
#M_H_agespec <- lands_clust[lands_clust$gene %in% tools::toTitleCase(tolower(H_agespec$gene)) & lands_clust$cond=="M",] ; dim(M_H_agespec) # 1793
common_Hagespec_M <- indexing_table_all[indexing_table_all$hsapiens_H %in% H_agespec$gene & indexing_table_all$M==TRUE , ] # 1818

# which ones are not matched exactly by name?
#diffs_H_O <- setdiff(H_agespec$gene,O_H_agespec$gene) # for example ATP6, we just have a landscape for human
# see how one human age specific looks like in other species
#L_MR_H["AASDHPPT", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2") # (e)
#L_MR_O["AASDHPPT", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2") # (l)
#L_MR_M["Aasdppt", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("SOX2") # (l)

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #



# HUMAN VS ORGANOID - do substraction maps for the common genes that are age specific in human
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# the list of common gene names will be used to reorder arrays so we ensure that rows are cross species correspondances
  
# filter and flatten the arrays ((faster computation) to matrices with 62500 columns instead of 3D arrays)
L_MR_H_f <- L_MR_H[rownames(L_MR_H) %in% common_Hagespec_O$hsapiens_H , , ] ; dim(L_MR_H_f)
L_MR_H_filt_flat <- matrix(aperm(L_MR_H_f, c(1,2,3)), nrow = dim(L_MR_H_f)[1]) ; dim(L_MR_H_filt_flat)
rownames(L_MR_H_filt_flat) <- rownames(L_MR_H_f)
# reorder to ensure that Human and Org matrices have the same order
L_MR_H_filt_flat <- L_MR_H_filt_flat[common_Hagespec_O$hsapiens_H, ]
head(rownames(L_MR_H_filt_flat))

L_MR_O_f <- L_MR_O[rownames(L_MR_O) %in% common_Hagespec_O$hsapiens_O, , ] ; dim(L_MR_O_f)
L_MR_O_filt_flat <- matrix(aperm(L_MR_O_f, c(1,2,3)), nrow = dim(L_MR_O_f)[1]) ; dim(L_MR_O_filt_flat)
rownames(L_MR_O_filt_flat) <- rownames(L_MR_O_f)
# reorder to ensure that Human and Org matrices have the same order
L_MR_O_filt_flat <- L_MR_O_filt_flat[common_Hagespec_O$hsapiens_O, ]
head(rownames(L_MR_O_filt_flat))

# calculate matrix differences, gene to gene (thats why the order was important)
n_chunks <- ceiling(nrow(L_MR_O_filt_flat) / 500)
diff_list_H_O <- pbmclapply(seq_len(n_chunks), function(i) { 
  message("Processing chunk ", i) # to print to console
  idx <- ((i - 1) * 500 + 1):min(i * 500, nrow(L_MR_O_filt_flat))
  list(L_MR_H_filt_flat[idx, , drop = FALSE], L_MR_O_filt_flat[idx, , drop = FALSE]) %>%
    lapply(function(mat) {
      mat %>%
        { (.-rowMins(.)) / (rowMaxs(.) - rowMins(.)) } %>%
        { replace(., is.nan(.), 0) } %>%
        { (.-rowMeans(.)) / rowSds(.) }
    }) %>%
    { .[[1]] - .[[2]] }
}, mc.cores = 20)
# diff_list_H_O is a list of length 5, where each element is a matrix of 500 × 62500
# reshape it
diff_matrix_H_O <- do.call(rbind, diff_list_H_O) ; rownames(diff_matrix_H_O) <- rownames(L_MR_H_filt_flat)   # 2013 × 62500
diff_array_H_O <- array(diff_matrix_H_O, dim = dim(L_MR_H_f)) ; dimnames(diff_array_H_O) <- list(rownames(L_MR_H_filt_flat), NULL, NULL) # 13564   250   250, like the shape of L_MR_O
#saveRDS(diff_array_H_O,"humous_v4/out/temporal_shifts/diff_array_H_O.rds")

# do and save substraction map plots - for vgg16, so without axes

# trial for one gene
  # first see the landscapes for each species (to make sure the SM makes sense)
gridExtra::grid.arrange(L_MR_H_f["AASDHPPT", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("AASDHPPT_H"),
                        L_MR_O_f["AASDHPPT", , ] %>% as.array() %>% rescale(to = c(0, 1)) %>%  raster::raster() %>%  base_ggraster() +  ggtitle("AASDHPPT_O"),ncol=2)
  # with legend and axes
plot(flip(raster::raster(diff_array_H_O["AASDHPPT",,]),direction="y"),breaks=seq(-10, 10, length.out=10),col=rev(colorRampPalette(c("green", "white", "red"))(length(seq(-10, 10, length.out=10))-1)))
  # with no legend nor axes and vgg16 dimensionsr <- raster(diff_array_H_O["AASDHPPT", , ])

r <- raster(diff_array_H_O["AASDHPPT", , ])
# Set breaks and colors
breaks <- seq(-10, 10, length.out = 10)
colors <- rev(colorRampPalette(c("green", "white", "red"))(length(breaks) - 1))
col_index <- cut(values(r), breaks = breaks, include.lowest = TRUE, labels = FALSE)
color_matrix <- matrix(colors[col_index], nrow = nrow(r), byrow = TRUE)
# Convert hex colors to RGB array
rgb <- col2rgb(color_matrix) / 255
dim(rgb) <- c(3, nrow(r), ncol(r))
rgb <- aperm(rgb, c(2, 3, 1))  # reshape to [height, width, channels]
rgb <- rgb[nrow(rgb):1, , ]    # flip vertically
#writePNG(rgb, target = "humous_v4/out/temporal_shifts/SMs_vgg16/SMs_HvsO/AASDHPPT_HvsO.png") # Save directly as PNG without any plotting


# now loop to plot and save all SMs HvsO (2013)
  # I am reusing breaks and colors cause they are always the same scale
for (i in rownames(diff_array_H_O)) {
  print(i)
  r <- raster(diff_array_H_O[i, , ])
  breaks <- seq(-10, 10, length.out = 10)
  colors <- rev(colorRampPalette(c("green", "white", "red"))(length(breaks) - 1))
  col_index <- cut(values(r), breaks = breaks, include.lowest = TRUE, labels = FALSE)
  color_matrix <- matrix(colors[col_index], nrow = nrow(r), byrow = TRUE)
  rgb <- col2rgb(color_matrix) / 255
  dim(rgb) <- c(3, nrow(r), ncol(r))
  rgb <- aperm(rgb, c(2, 3, 1))  # reshape to [height, width, channels]
  rgb <- rgb[nrow(rgb):1, , ]    # flip vertically
  writePNG(rgb, target = paste0("humous_v4/out/temporal_shifts/SMs_vgg16/SMs_HvsO/", i, "_HvsO.png"))
}


# besides SMs, see if there is a numerical indication of bias vs noise
# mean
mean_map <- apply(diff_array_H_O, c(2, 3), mean) # mean pixel value across all genes - see if some areas are consistently negative or positive in average
mean_map <- reshape2::melt(mean_map)
colnames(mean_map) <- c("x", "y", "mean_value")
ggplot(mean_map, aes(x = x, y = y, fill = mean_value)) +
  geom_raster() +
  scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0,limits = c(-10, 10)) + # SMs and the diff matrix is in range -10 to +10
  theme_minimal() +
  ggtitle("Mean value across all genes")
#saveRDS(mean_map,"humous_v4/out/temporal_shifts/SMs_vgg16/mean_map_HvsO.rds")

# sd
sd_map <- apply(diff_array_H_O, c(2, 3), sd)
sd_map <- reshape2::melt(sd_map)
colnames(sd_map) <- c("x", "y", "sd_value")
ggplot(sd_map, aes(x = x, y = y, fill = sd_value)) +
  geom_raster() +
  scale_fill_gradient(low = "white", high = "red") + 
  theme_minimal() +
  ggtitle("Mean value across all genes")
#saveRDS(sd_map,"humous_v4/out/temporal_shifts/SMs_vgg16/sd_map_HvsO.rds")


# signal to noise ratio (combine mean map and sd map) - signed SNR = mean / sd
snr_signed_map_HvsO <- apply(diff_array_H_O, c(2, 3), mean) / (apply(diff_array_H_O, c(2, 3), sd) + 1e-8) # add a small epsilon to avoid division by zero
snr_signed_map_HvsO <- reshape2::melt(snr_signed_map_HvsO)
colnames(snr_signed_map_HvsO) <- c("x", "y", "snr_signed_value")
#saveRDS(snr_signed_map_HvsO,"humous_v4/out/temporal_shifts/SMs_vgg16/snr_signed_map_HvsO.rds")

ggplot(snr_signed_map_HvsO, aes(x = x, y = y, fill = snr_signed_value)) +
  geom_raster() +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red", midpoint = 0) +
  theme_minimal() + ggtitle("Signed Signal-to-Noise Ratio (mean / SD) map")
#ggsave("humous_v4/out/temporal_shifts/SMs_vgg16/SignalToNoise_map_SMs_HvsO.pdf",useDingbats=FALSE)

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #



# HUMAN VS MOUSE - do substraction maps for the 1818 common genes that are age specific in human
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# find the common genes in H and O and create vector - use it to reorder so the order of rows is preserved
common_Hagespec_M <- indexing_table_all[indexing_table_all$hsapiens_H %in% H_agespec$gene & indexing_table_all$M==TRUE , ]

# filter and flatten the arrays ((faster computation) to matrices with 62500 columns instead of 3D arrays)
L_MR_H_f <- L_MR_H[rownames(L_MR_H) %in% common_Hagespec_M$hsapiens_H, , ] ; dim(L_MR_H_f)
L_MR_H_filt_flat <- matrix(aperm(L_MR_H_f, c(1,2,3)), nrow = dim(L_MR_H_f)[1]) ; dim(L_MR_H_filt_flat)
rownames(L_MR_H_filt_flat) <- rownames(L_MR_H_f)
# reorder to ensure that Human and Org matrices have the same order
L_MR_H_filt_flat <- L_MR_H_filt_flat[common_Hagespec_M$hsapiens_H, ]
head(rownames(L_MR_H_filt_flat))

L_MR_M_f <- L_MR_M[rownames(L_MR_M) %in% common_Hagespec_M$mmusculus, , ] ; dim(L_MR_M_f)
L_MR_M_filt_flat <- matrix(aperm(L_MR_M_f, c(1,2,3)), nrow = dim(L_MR_M_f)[1]) ; dim(L_MR_M_filt_flat)
rownames(L_MR_M_filt_flat) <- rownames(L_MR_M_f)
# reorder to ensure that Human and Org matrices have the same order
L_MR_M_filt_flat <- L_MR_M_filt_flat[common_Hagespec_M$mmusculus, ]
head(rownames(L_MR_M_filt_flat))

# calculate matrix differences, gene to gene (thats why the order was important)
n_chunks <- ceiling(nrow(L_MR_M_filt_flat) / 500)
diff_list_H_M <- pbmclapply(seq_len(n_chunks), function(i) { 
  message("Processing chunk ", i) # to print to console
  idx <- ((i - 1) * 500 + 1):min(i * 500, nrow(L_MR_M_filt_flat))
  list(L_MR_H_filt_flat[idx, , drop = FALSE], L_MR_M_filt_flat[idx, , drop = FALSE]) %>%
    lapply(function(mat) {
      mat %>%
        { (.-rowMins(.)) / (rowMaxs(.) - rowMins(.)) } %>%
        { replace(., is.nan(.), 0) } %>%
        { (.-rowMeans(.)) / rowSds(.) }
    }) %>%
    { .[[1]] - .[[2]] }
}, mc.cores = 20)
# diff_list_H_O is a list of length 5, where each element is a matrix of 500 × 62500
# reshape it
diff_matrix_H_M <- do.call(rbind, diff_list_H_M) ; rownames(diff_matrix_H_M) <- rownames(L_MR_H_filt_flat)   # 2013 × 62500
diff_array_H_M <- array(diff_matrix_H_M, dim = dim(L_MR_H_f)) ; dimnames(diff_array_H_M) <- list(rownames(L_MR_H_filt_flat), NULL, NULL) # 13564   250   250, like the shape of L_MR_O
#saveRDS(diff_array_H_M,"humous_v4/out/temporal_shifts/diff_array_H_M.rds")
diff_array_H_M <- readRDS("~/humous/humous_v4/out/temporal_shifts/diff_array_H_M.rds")

# do and save substraction map plots - for vgg16, so without axes

# now loop to plot and save all SMs HvsM 
# I am reusing breaks and colors cause they are always the same scale
for (i in rownames(diff_array_H_M)) {
  print(i)
  r <- raster(diff_array_H_M[i, , ])
  breaks <- seq(-10, 10, length.out = 10)
  colors <- rev(colorRampPalette(c("green", "white", "red"))(length(breaks) - 1))
  col_index <- cut(values(r), breaks = breaks, include.lowest = TRUE, labels = FALSE)
  color_matrix <- matrix(colors[col_index], nrow = nrow(r), byrow = TRUE)
  rgb <- col2rgb(color_matrix) / 255
  dim(rgb) <- c(3, nrow(r), ncol(r))
  rgb <- aperm(rgb, c(2, 3, 1))  # reshape to [height, width, channels]
  rgb <- rgb[nrow(rgb):1, , ]    # flip vertically
  writePNG(rgb, target = paste0("humous_v4/out/temporal_shifts/SMs_vgg16/SMs_HvsM/", i, "_HvsM.png"))
}


# besides SMs, see if there is a numerical indication of bias vs noise
# mean
mean_map_HvsM <- apply(diff_array_H_M, c(2, 3), mean) # mean pixel value across all genes - see if some areas are consistently negative or positive in average
mean_map_HvsM <- reshape2::melt(mean_map_HvsM)
colnames(mean_map_HvsM) <- c("x", "y", "mean_value")
#saveRDS(mean_map_HvsM,"humous_v4/out/temporal_shifts/SMs_vgg16/mean_map_HvsM.rds")
mean_map_HvsM <- readRDS("~/humous/humous_v4/out/temporal_shifts/SMs_vgg16/mean_map_HvsM.rds")

ggplot(mean_map_HvsM, aes(x = x, y = y, fill = mean_value)) +
  geom_raster() + scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0) + 
  theme_minimal() + ggtitle("Mean value HvsM substraction across all landscapes that are age specific in H")
#ggsave("humous_v4/out/temporal_shifts/SMs_vgg16/mean_map_SMs_HvsM.pdf",useDingbats=FALSE)

# sd
sd_map_HvsM <- apply(diff_array_H_M, c(2, 3), sd)
sd_map_HvsM <- reshape2::melt(sd_map_HvsM)
colnames(sd_map_HvsM) <- c("x", "y", "sd_value")
#saveRDS(sd_map,"humous_v4/out/temporal_shifts/SMs_vgg16/sd_map_HvsM.rds")
sd_map_HvsM <- readRDS("~/humous/humous_v4/out/temporal_shifts/SMs_vgg16/sd_map_HvsM.rds")

ggplot(sd_map_HvsM, aes(x = x, y = y, fill = sd_value)) +
  geom_raster() +  scale_fill_gradient(low = "white", high = "red") +  theme_minimal() + 
  ggtitle("SD value HvsM substraction across all landscapes that are age specific in H")
#ggsave("humous_v4/out/temporal_shifts/SMs_vgg16/sd_map_SMs_HvsM.pdf",useDingbats=FALSE)


# signal to noise ratio (combine mean map and sd map) - signed SNR = mean / sd
snr_signed_map_HvsM <- apply(diff_array_H_M, c(2, 3), mean) / (apply(diff_array_H_M, c(2, 3), sd) + 1e-8) # add a small epsilon to avoid division by zero
snr_signed_map_HvsM <- reshape2::melt(snr_signed_map_HvsM)
colnames(snr_signed_map_HvsM) <- c("x", "y", "snr_signed_value")
#saveRDS(snr_signed_map_HvsM,"humous_v4/out/temporal_shifts/SMs_vgg16/snr_signed_map_HvsM.rds")
snr_signed_map_HvsM <- readRDS("~/humous/humous_v4/out/temporal_shifts/SMs_vgg16/snr_signed_map_HvsM.rds")

ggplot(snr_signed_map_HvsM, aes(x = x, y = y, fill = snr_signed_value)) +
  geom_raster() +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red", midpoint = 0) +
  theme_minimal() + ggtitle("Signed Signal-to-Noise Ratio (mean / SD) value HvsM substraction across all landscapes that are age specific in H")
#ggsave("humous_v4/out/temporal_shifts/SMs_vgg16/SignalToNoise_map_SMs_HvsM.pdf",useDingbats=FALSE)

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #







