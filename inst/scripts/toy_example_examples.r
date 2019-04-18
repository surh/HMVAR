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
library(magrittr)

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



############

m <- dplyr::tibble(sample = letters[1:4], Group = rep(LETTERS[1:2], each = 2))
m
info
i <- dplyr::tibble(site_id = 's1')
i

d <- dplyr::tibble(site_id = i$site_id, a = 1, b = 1, c = 1 , d = 1)
d

f <- dplyr::tibble(site_id = i$site_id, a = 0, b = 0, c = 1, d = 1)
f


determine_snp_dist(info = i, freq = f, depth = d, map = m, depth_thres = 1, freq_thres = 0.5)

e <- tibble(site_id = i$site_id, distribution = "Fixed")
e
