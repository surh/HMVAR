#!/usr/bin/env Rscript

# (C) Copyright 2018 Sur Herrera Paredes
# 
# This file is part of HMVAR.
# 
# HMVAR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# HMVAR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with HMVAR.  If not, see <http://www.gnu.org/licenses/>.

library(HMVAR)
library(tidyverse)
library(argparser)


# 
# mktest_file <- "../2018-12-14.hmp_mktest/results/Granulicatella_adiacens_61980_mktest.txt"
# Genes <- read_tsv(mktest_file)
# gene <- Genes$gene[2]
genes <- c("638301.3.peg.1",
           "638301.3.peg.444",
           "638301.3.peg.884",
           "638301.3.peg.944",
           "638301.3.peg.955",
           "638301.3.peg.1091",
           "638301.3.peg.1273",
           "638301.3.peg.1346",
           "638301.3.peg.1361")

# Parameters
midas_dir <- "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/Granulicatella_adiacens_61980/"
depth_thres <- 1
map_file <- "/godot/users/sur/exp/fraserv/2018/2018-12-14.hmp_mktest/hmp_SPvsTD_map.txt"
freq_thres <- 0.5


mktest <- midas_mktest(midas_dir = midas_dir,
                       map_file = map_file,
                       genes = genes,
                       depth_thres = depth_thres,
                       freq_thres = 0.5)
mktest

######### Placing function code ############
map <- readr::read_tsv(map_file)
map <- map %>% dplyr::select(sample = ID, Group)
Dat <- read_midas_data(midas_dir = midas_dir, map = map, 
                       genes = genes)
Dat$info <- determine_snp_effect(Dat$info)
# Works up to here

# Dat$info <- determine_snp_dist(info = Dat$info, freq = Dat$freq, 
#                                depth = Dat$depth, map = map,
#                                depth_thres = depth_thres, 
#                                freq_thres = freq_thres)


info <- Dat$info
freq <- Dat$freq
depth <- Dat$depth


freq_thres <- 0.2

if (freq_thres < 0 || freq_thres > 1) 
  stop("ERROR: freq_thres must have values in [0, 1]", 
       call. = TRUE)
freq_thres <- min(freq_thres, 1 - freq_thres)
depth <- depth %>% tidyr::gather(key = "sample", value = "depth", 
                                 -site_id)
freq <- freq %>% tidyr::gather(key = "sample", value = "freq", 
                               -site_id)

dat <- depth %>%
  dplyr::inner_join(freq, by = c("site_id", "sample")) %>% 
  dplyr::left_join(map, by = "sample") %>% 
  dplyr::filter(depth >= depth_thres)
# dat

dat <- dat %>%
  dplyr::mutate(allele = replace(freq, freq < freq_thres, "major")) %>%
  dplyr::mutate(allele = replace(allele, freq >= (1 - freq_thres), "minor"))
# dat

dat <- dat %>%
  dplyr::mutate(allele = replace(allele, (freq >= freq_thres) & (freq < (1 - freq_thres)), NA))
dat

  dplyr::filter(!is.na(allele))
dat

site_dist <- dat %>% split(.$site_id) %>% purrr::map_chr(get_site_dist)
site_dist <- tibble::tibble(site_id = names(site_dist), distribution = factor(site_dist, 
                                                                              levels = c("Fixed", "Invariant", "Polymorphic")))
info <- info %>% dplyr::inner_join(site_dist, by = "site_id")
return(info)








Res <- Dat$info %>% split(.$gene_id) %>% purrr::map_dfr(mkvalues, 
                                                        depth_thres = depth_thres, .id = "gene_id")
return(Res)






vmwa <- read_tsv("hmp.vmwa.pvals.genes/Granulicatella_adiacens_61980_vmwa.genes.txt")
vmwa <- vmwa %>% filter(gene_id %in% genes)
vmwa %>% arrange(q.value)
