##3 Jan 2024


read_seurat <- function(file_path) {
  seurat <- readRDS(file_path)
  return(seurat)
}

# Function to process Change-O database files for each lane
process_lane_data <- function(dir, lane, metadata) {
  h <- readChangeoDb(file.path(dir, "filtered_contig_heavy_productive-T.tsv"))
  l <- readChangeoDb(file.path(dir, "filtered_contig_light_productive-T.tsv"))
  
  comb <- bind_rows(h, l)
  
  meta <- filter(metadata, Lane == lane)
  
  comb$orig.barcode <- comb$cell_id
  comb$cell_id <- paste0(unlist(lapply(strsplit(comb$cell_id, split = "-"), function(x) x[1])), "-", lane)
  
  comb <- filter(comb, cell_id %in% rownames(meta))
  meta$mcell_id <- rownames(meta)
  
  m <- match(comb$cell_id, meta$mcell_id)
  
  comb <- bind_cols(comb, meta[m,])
  
  return(comb)
}

# Function to filter productive data
filter_productive_data <- function(data) {
  data <- filter(data, productive == TRUE)
  return(data)
}

# Function to write processed data to Change-O database file
write_changeo_db <- function(data, file_path) {
  writeChangeoDb(data, file = file_path)
}

# Function to calculate SHM frequency and write germline data to file
calculate_shm_frequency <- function(comb_germline, file_path) {
  comb_germline <- observedMutations(comb_germline, sequenceColumn = "sequence_alignment",
                                     germlineColumn = "germline_alignment_d_mask",
                                     regionDefinition = IMGT_V,
                                     frequency = TRUE,
                                     combine = TRUE, 
                                     nproc = 3)
  writeChangeoDb(comb_germline, file_path)
}

# Function to generate plots for distance analysis
generate_distance_plots <- function(dist_cross) {
  pdf("intermediates/crossDistance.pdf", height = 20, width = 8)
  ggplot(subset(dist_cross, !is.na(cross_dist_nearest)), 
         aes(x = cross_dist_nearest)) + 
    theme_bw() + 
    xlab("Cross-sample_id Hamming distance") + 
    ylab("Count") +
    geom_histogram(color = "white", binwidth = 0.02) +
    geom_vline(xintercept = 0.12, color = "firebrick", linetype = 2) +
    facet_grid(Lane ~ ., scales = "free_y")
  dev.off()
}

# Function to perform clonal analysis
perform_clonal_analysis <- function(data) {
  clones <- data %>%
    group_by(Lane) %>%
    do(as.data.frame(
      hierarchicalClones(., 0.1, cell_id = "cell_id", locus = "locus",
                         only_heavy = FALSE, split_light = TRUE, cdr3 = TRUE, nproc = 2,
                         verbose = FALSE, log = NULL, summarize_clones = FALSE)))
  
  clones <- ungroup(clones)
  clones$clone_id <- paste0(clones$Donor, "-", clones$clone_id)
  
  return(clones)
}

# Function to create germlines
create_germlines <- function(comb, references) {
  h <- createGermlines(filter(comb, locus == "IGH"), references)
  k <- createGermlines(filter(comb, locus == "IGK"), references)
  l <- createGermlines(filter(comb, locus == "IGL"), references)
  
  comb_germline <- bind_rows(h, l, k)
  
  return(comb_germline)
}



