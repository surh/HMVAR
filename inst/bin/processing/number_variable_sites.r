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
             map_file = "/godot/users/sur/exp/fraserv/2018/2018-12-14.hmp_mktest/hmp_SPvsTD_map.txt",
             genes = NA,
             depth_thres = 1,
             freq_thres = 0.5,
             cds_only = FALSE,
             plot_positin = TRUE,
             prefix = "Granulicatella_adiacens_61980",
             outdir = "results/")

# Read genes
if(is.na(args$genes)){
  genes <- NULL
}else{
  genes <- read_tsv(args$genes, col_names = FALSE, col_types = 'c')
  genes <- genes$X1
}

# Read map
map <- read_tsv(args$map_file)
map <- map %>%
  select(sample = ID,
         everything())

# Read MIDAS data
info <- read_tsv(paste0(args$midas_dir, "/snps_info.txt"),
                 col_types = 'ccncccnnnnnccccc',
                 na = 'NA')
depth <- read_midas_abun(paste0(args$midas_dir, "/snps_depth.txt"))
freq <- read_midas_abun(paste0(args$midas_dir, "/snps_freq.txt"))

# Process data
# Clean info
info <- info %>% 
  select(-locus_type, -starts_with("count_"), -ref_allele)
# Clean depth and freq
depth <- depth %>%
  select(site_id, intersect(map$sample, colnames(depth)) )
freq <- freq %>%
  select(site_id, intersect(map$sample, colnames(freq)) )
# Clean map
map <- map %>% 
  filter(sample %in% colnames(depth))

# Select gene/cds data
if(args$cds_only){
  info <- info %>% 
    filter(!is.na(gene_id))
}
if(!is.null(genes)){
  info <- info %>%
    filter(gene_id %in% genes)
}
freq <- freq %>% 
  filter(site_id %in% info$site_id)
depth <- depth %>% 
  filter(site_id %in% info$site_id)


# Calcualate snp effect
info <- determine_snp_effect(info)
# # Calculate snp dist
# info <- determine_snp_dist(info = info,
#                            freq = freq,
#                            depth = depth,
#                            map = map,
#                            depth_thres = args$depth_thres,
#                            freq_thres = args$freq_thres)

# Subsitution type
info <- info %>%
  add_column(substitution = info %>%
               pmap_chr(function(major_allele, minor_allele, ...){
                 base_type <- c(A = "purine", C = "pyrimidine", G = "purine", T = "pyrimidine")
                 if(base_type[major_allele] == base_type[minor_allele]){
                   substitution <- "transition" 
                 }else{
                   substitution <- "transversion"
                 }
                 return(substitution)
               }))



# Match freqs and depth
depth <- depth %>% gather(key = "sample", value = 'depth', -site_id)
freq <- freq %>% gather(key = "sample", value = 'freq', -site_id)
# meta <- Dat$info %>% select(site_id, ref_id, ref_pos, snp_effect, distribution)

dat <- depth %>%
  inner_join(freq, by = c("site_id", "sample")) %>%
  left_join(map, by = "sample") %>%
  filter(depth >= args$depth_thres) %>%
  left_join(info, by = "site_id") %>%
  filter(distribution != "Invariant")
# dat


# Count number of sites and number of variable per sample
varsites <- dat %>%
  filter(depth >= 1) %>%
  split(.$sample) %>%
  map_dfr(function(d){
    nsites <- nrow(d)
    sites <- d$freq
    nfixed <- sum(sites >= 1 | sites <= 0)
    ntransitions <- sum(d$substitution == "transition")
    nsynonymous <- sum(d$snp_effect == "synonymous")
    nnonsynonymous <- sum(d$snp_effect == "non-synonymous")
    return(tibble(nsites = nsites,
                  fixed = nfixed,
                  variable = nsites - nfixed,
                  transitions = ntransitions,
                  transversions = nsites - ntransitions,
                  synonymous = nsynonymous,
                  non.synonymous = nnonsynonymous,
                  IGR = nsites - nsynonymous - nnonsynonymous))
  }, .id = "sample") %>%
  inner_join(map, by = "sample")
varsites

# Plot total. variable and fixed sites
p1 <- varsites %>%
  ggplot(aes(x = Group, y = nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
p1
# ggsave("Granulicatella_adiacens_61980_totsites.png", p1, width = 5, height = 4, dpi = 150)

p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*variable/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1
# ggsave("Granulicatella_adiacens_61980_varsites_depth1.png", p1, width = 5, height = 4, dpi = 150)

p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*fixed/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1

# Plot transitions and transversions
p1 <- varsites %>%
  ggplot(aes(x = Group, y = transitions, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1

p1 <- varsites %>%
  ggplot(aes(x = Group, y = transversions, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1

p1 <- varsites %>%
  ggplot(aes(x = Group, y = transitions/transversions, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_log10() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1


# Plot synonymous and non-synonymous
p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*synonymous/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_continuous(limits = c(0,100)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1

p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*non.synonymous/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_continuous(limits = c(0,100)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1

p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*IGR/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_continuous(limits = c(0,100)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1

p1 <- varsites %>%
  ggplot(aes(x = Group, y = non.synonymous/synonymous, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_log10() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1

#### Plot position






ggplot(dat, aes(x = ref_pos ))



pos.varsites <- dat %>%
  split(.$site_id) %>%
  map_dfr(function(d){
    res <- d %>%
      split(.$Group) %>% map_dfr(function(d){
      nsamples <- nrow(d)
      nfixed <- sum(d$freq == 0 | d$freq == 1)
      return(tibble(site_id = unique(d$site_id),
                    nsamples = nsamples,
                    fixed = nfixed))
    }, .id = "Group")
  }) %>% left_join(info, by = "site_id")

pos.varsites


p1 <- ggplot(pos.varsites, aes(x = ref_pos, y = fixed / nsamples)) +
  facet_grid(~ ref_id, space = "free_x") +
  # geom_line(aes(color = Group)) +
  geom_point(aes(color = Group), size = 0.1) +
  geom_smooth(aes(color = Group), se = FALSE, method = "gam", formula = y ~ s(x, bs = "cs")) +
  scale_y_continuous(limits = c(0,1)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line())
p1





