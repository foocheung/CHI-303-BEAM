library(alakazam)
library(shazam)
library(dplyr)
library(tidyr)
library(scoper)
library(dowser)
library(Biostrings)
source("./BCR_functions.R")

setwd("/Users/cheungf/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/R/BEAM/IMMCANTATION")

seurat = read_seurat("seurat.RDS")
metadata = seurat@meta.data

lanes = paste0(c(1,2))

data = tibble()
for(lane in lanes){
  dir = paste0("./sample_", lane, sep = "")
  print(dir)
  lane_data = process_lane_data(dir, lane, metadata)
  data = bind_rows(data, lane_data)
}
data = filter_productive_data(data)
write_changeo_db(data, file = "processed/alldata.tsv")

alldata_path <- "./processed/alldata.tsv"
data = readChangeoDb(file = alldata_path)
data = filter_productive_data(data)

dist_cross = distToNearest(filter(data, locus == "IGH"),
                           sequenceColumn = "junction", 
                           vCallColumn = "v_call", jCallColumn = "j_call",
                           model = "ham", normalize = "len", nproc = 1,
                           cross = "Lane")

generate_distance_plots(dist_cross)

multi_heavy = table(filter(data, locus == "IGH")$cell_id)
multi_heavy_cells = names(multi_heavy)[multi_heavy > 1]
data = filter(data, !cell_id %in% multi_heavy_cells)

clones = perform_clonal_analysis(data)

references = readIMGT(dir = "vdj")
comb_germline = create_germlines(clones, references)

calculate_shm_frequency(comb_germline, "processed/all_cloned_data.tsv")
