#!/usr/bin/env Rscript

# (C) Copyright 2019 Sur Herrera Paredes
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

# setwd("/cashew/users/sur/exp/fraserv/2019/today2")
library(HMVAR)
library(tidyverse)

metawas_file <- commandArgs(trailingOnly = TRUE)[1]
info_file <- commandArgs(trailingOnly = TRUE)[2]

# metawas_file <- "/cashew/users/sur/data/gathered_results/2019a.hmp.subsite/metawas/Supragingival.plaque/metawas/lmm/Actinomyces_odontolyticus_57475_lmm.assoc.txt"
# info_file <- "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/Actinomyces_odontolyticus_57475/snps_info.txt"

info <- readr::read_tsv(paste0(info_file), 
                        col_types = readr::cols(.default = readr::col_character(), 
                                                ref_pos = readr::col_number(), count_samples = readr::col_number(), 
                                                count_a = readr::col_number(), count_c = readr::col_number(), 
                                                count_g = readr::col_number(), count_t = readr::col_number()), 
                        na = "NA") %>%
  select(-starts_with("count_"), -locus_type)

metawas <- read_tsv(metawas_file,
                    col_types = cols(.default = col_character(),
                                     ps = col_number(),
                                     n_miss = col_number(),
                                     af = col_number(),
                                     logl_H1 = col_double(),
                                     l_mle = col_number(),
                                     p_lrt = col_number()),
                    na = c("", "na", "nan", "inf", "NA", "-nan"))
metawas <- metawas %>%
  select(site_id = rs, ref_id = chr, ref_pos = ps, everything()) %>%
  left_join(info, by = c("site_id", "ref_id", "ref_pos")) %>%
  select(-allele1, -allele0, -logl_H1, -l_mle) %>%
  determine_snp_effect %>%
  mutate(snp_effect = replace(as.character(snp_effect), is.na(snp_effect), 'non-coding')) %>%
  mutate(snp_effect = factor(snp_effect, levels = c("non-coding", 'synonymous', 'non-synonymous')))
# metawas

# p1 <- ggplot(metawas, aes(x = af, y = p_lrt)) +
#   geom_point() +
#   scale_y_continuous(trans = scales::trans_new(name = 'minus_log',
#                                                transform = function(x) -log10(x),
#                                                inverse = function(x) 10^(-x),
#                                                breaks = function(x){
#                                                  rng <- range(x)
#                                                  breaks <- seq(rng[1], rng[2], length.out = 5)
#                                                  10^(-breaks)
#                                                })) +
#   AMOR::theme_blackbox()

p1 <- ggplot(metawas, aes(x = af, y = p_lrt, col = snp_effect)) +
  geom_point() +
  scale_y_continuous(trans = scales::trans_new(name = 'minus_log',
                                               transform = function(x) -log10(x),
                                               inverse = function(x) 10^(-x),
                                               breaks = function(x){
                                                 rng <- range(x)
                                                 rng[is.infinite(rng)] <- 0
                                                 breaks <- seq(rng[1], rng[2], length.out = 5)
                                                 10^(-breaks)
                                               })) +
  AMOR::theme_blackbox()
ggsave("freq_vs_pval.png", p1, width = 5, height = 4, dpi = 150)

p1 <- ggplot(metawas, aes(x = af, y = p_lrt, col = snp_effect)) +
  facet_wrap(~ snp_effect) +
  geom_point() +
  scale_y_continuous(trans = scales::trans_new(name = 'minus_log',
                                               transform = function(x) -log10(x),
                                               inverse = function(x) 10^(-x),
                                               breaks = function(x){
                                                 rng <- range(x)
                                                 rng[is.infinite(rng)] <- 0
                                                 breaks <- seq(rng[1], rng[2], length.out = 5)
                                                 10^(-breaks)
                                               })) +
  AMOR::theme_blackbox()
ggsave("freq_vs_pval.facets.png", p1, width = 6, height = 3, dpi = 150)

p1 <- metawas %>%
  split(.$snp_effect) %>%
  map_dfr(function(d){
    tibble(n_total = nrow(d),
           n_sig = sum(d$p_lrt < 1e-5))
  }, .id = "snp_effect") %>%
  mutate(snp_effect = factor(snp_effect, levels = c("non-coding", "synonymous", "non-synonymous"))) %>%
  ggplot(aes(x = snp_effect, y = 100*n_sig/n_total, fill = snp_effect)) +
  geom_bar(stat = 'identity') +
  AMOR::theme_blackbox()
ggsave("sig_freqs.png", p1, width = 5, height = 4, dpi = 150)
write_tsv(p1$data, "sig_nums.txt")


