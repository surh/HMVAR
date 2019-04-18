###########################
# read_midas_data
library(HMVAR)
library(tidyverse)

# Get file paths
midas_dir <- system.file("toy_example/merged.snps/", package = "HMVAR")
map <- read_tsv(system.file("toy_example/map.txt", package = "HMVAR"),
                col_types = cols(.default = col_character())) %>%
  select(sample = ID, Group)

# Read data
midas_data <- read_midas_data(midas_dir = midas_dir, map = map, cds_only = TRUE)
midas_data

########################
# midas_mktest
library(HMVAR)

# Get paths
midas_dir <- system.file("toy_example/merged.snps/", package = "HMVAR")
map_file <- system.file("toy_example/map.txt", package = "HMVAR")

# Process map yourself
map <- readr::read_tsv(map_file,
                       col_types = readr::cols(.default = readr::col_character())) %>%
  dplyr::select(sample = ID, Group)
midas_mktest(midas_dir = midas_dir,
             map = map)

# Give a map file path to midas_mktest
midas_mktest(midas_dir = midas_dir,
             map = map_file)

############################
# determine_snp_effect, determine_snp_dist, mkvalues
library(HMVAR)

# Get file paths
midas_dir <- system.file("toy_example/merged.snps/", package = "HMVAR")
map <- readr::read_tsv(system.file("toy_example/map.txt", package = "HMVAR"),
                       col_types = readr::cols(.default = readr::col_character())) %>%
  dplyr::select(sample = ID, Group)

# Read data
midas_data <- read_midas_data(midas_dir = midas_dir, map = map, cds_only = TRUE)

info <- determine_snp_effect(midas_data$info) %>%
  determine_snp_dist(freq = midas_data$freq,
                     depth = midas_data$depth, map = map,
                     depth_thres = 1, freq_thres = 0.5)
info

mktable <- info %>%
  split(.$gene_id) %>%
  purrr::map_dfr(mkvalues,
                 .id = "gene_id")
mktable
                 