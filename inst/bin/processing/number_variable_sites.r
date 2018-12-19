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

# setwd("/godot/users/sur/exp/fraserv/2018/today3/")

# Parameters
args <- list(midas_dir = "/godot/users/sur/exp/fraserv/2018/2018-12-14.hmp_mktest/genomes/Granulicatella_adiacens_61980/",
             depth_thres = 1,
             map_file = "/godot/users/sur/exp/fraserv/2018/2018-12-14.hmp_mktest/hmp_SPvsTD_map.txt",
             genes = NULL,
             depth_thres = 1,
             freq_thres = 0.5)

# midas_dir <- "../2018-12-14.hmp_mktest/genomes/Granulicatella_adiacens_61980/"
# depth_thres <- 1
# map_file <- "../2018-12-14.hmp_mktest/hmp_SPvsTD_map.txt"
# genes <- "638301.3.peg.283"

map <- read_tsv(args$map_file)
map <- map %>%
  select(sample = ID,
         everything())

Dat <- read_midas_data(midas_dir = args$midas_dir,
                       map = map,
                       genes = args$genes)

# Calcualate snp effect
Dat$info <- determine_snp_effect(Dat$info)
# Calculate snp dist
Dat$info <- determine_snp_dist(info = Dat$info,
                               freq = Dat$freq,
                               depth = Dat$depth,
                               map = map,
                               depth_thres = args$depth_thres,
                               freq_thres = args$freq_thres)


# Match freqs and depth
depth <- Dat$depth %>% gather(key = "sample", value = 'depth', -site_id)
freq <- Dat$freq %>% gather(key = "sample", value = 'freq', -site_id)
meta <- Dat$info %>% select(site_id, ref_id, ref_pos, snp_effect, distribution)

dat <- depth %>%
  inner_join(freq, by = c("site_id", "sample")) %>%
  left_join(map, by = "sample") %>%
  filter(depth >= depth_thres) %>%
  left_join(meta, by = "site_id") %>%
  mutate(allele = replace(freq, freq < 0.5, 'major')) %>%
  mutate(allele = replace(allele, freq >= 0.5, 'minor')) %>%
  filter(distribution != "Invariant")
dat


# Count number of sites and number of variable per sample
varsites <- dat %>%
  filter(depth >= 1) %>%
  split(.$sample) %>%
  map_dfr(function(d){
    nsites <- nrow(d)
    sites <- d$freq
    nfixed <- sum(sites >= 1 | sites <= 0)
    return(tibble(nsites = nsites,
                  nfixed = nfixed,
                  nvariable = nsites - nfixed))
  }, .id = "sample") %>%
  inner_join(map, by = "sample")
varsites

p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*nvariable/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
ggsave("Granulicatella_adiacens_61980_varsites_depth1.png", p1, width = 5, height = 4, dpi = 150)

p1 <- varsites %>%
  ggplot(aes(x = Group, y = nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
ggsave("Granulicatella_adiacens_61980_totsites.png", p1, width = 5, height = 4, dpi = 150)

# varsites %>%
#   ggplot(aes(x = Group, y = nvariable, col = Group)) +
#   geom_boxplot() +
#   geom_point(position = position_jitter(width = 0.2)) +
#   theme(panel.background = element_blank(),
#         panel.grid = element_blank())


varsites <- dat %>%
  filter(depth >= 15) %>%
  split(.$sample) %>%
  map_dfr(function(d){
    nsites <- nrow(d)
    sites <- d$freq
    nfixed <- sum(sites >= 1 | sites <= 0)
    return(tibble(nsites = nsites,
                  nvariable = nsites - nfixed))
  }, .id = "sample") %>%
  inner_join(map, by = "sample")
varsites

p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*nvariable/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
p1
ggsave("Granulicatella_adiacens_61980_varsites_depth15.png",
       p1, width = 5, height = 4, dpi = 150)

p1 <- varsites %>%
  ggplot(aes(x = Group, y = nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
p1
ggsave("Granulicatella_adiacens_61980_totsites_depth15.png",
       p1, width = 5, height = 4, dpi = 150)
