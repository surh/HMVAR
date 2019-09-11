setwd("/cashew/users/sur/exp/fraserv/2019/today2")
library(HMVAR)
library(tidyverse)

metawas_file <- "/cashew/users/sur/data/gathered_results/2019a.hmp.subsite/metawas/Supragingival.plaque/metawas/lmm/Actinomyces_odontolyticus_57475_lmm.assoc.txt"
midas_dir <- "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/Actinomyces_odontolyticus_57475/"
map_file <- "/cashew/users/sur/data/gathered_results/2019a.hmp.subsite/hmp.subsite_map.txt"

map <- read_tsv(map_file) %>%
  select(sample = ID, Group = Group)

Dat <- read_midas_data(midas_dir = midas_dir, map = map, cds_only = FALSE)


metawas <- read_tsv(metawas_file,
                    col_types = cols(.default = col_character(),
                                     ps = col_number(),
                                     n_miss = col_number(),
                                     af = col_number(),
                                     logl_H1 = col_number(),
                                     l_mle = col_number(),
                                     p_lrt = col_number()))
metawas <- metawas %>%
  select(site_id = rs, ref_id = chr, ref_pos = ps, everything()) %>%
  left_join(Dat$info, by = c("site_id", "ref_id", "ref_pos")) %>%
  select(-allele1, -allele0, -logl_H1, -l_mle) %>%
  determine_snp_effect %>%
  mutate(snp_effect = replace(as.character(snp_effect), is.na(snp_effect), 'non-coding')) %>%
  mutate(snp_effect = factor(snp_effect, levels = c("non-coding", 'synonymous', 'non-synonymous')))
metawas


p1 <- ggplot(metawas, aes(x = af, y = p_lrt)) +
  geom_point() +
  scale_y_continuous(trans = scales::trans_new(name = 'minus_log',
                                               transform = function(x) -log10(x),
                                               inverse = function(x) 10^(-x),
                                               breaks = function(x){
                                                 rng <- range(x)
                                                 breaks <- seq(rng[1], rng[2], length.out = 5)
                                                 10^(-breaks)
                                               })) +
  AMOR::theme_blackbox()
p1

p1 <- ggplot(metawas, aes(x = af, y = p_lrt, col = snp_effect)) +
  geom_point() +
  scale_y_continuous(trans = scales::trans_new(name = 'minus_log',
                                               transform = function(x) -log10(x),
                                               inverse = function(x) 10^(-x),
                                               breaks = function(x){
                                                 rng <- range(x)
                                                 breaks <- seq(rng[1], rng[2], length.out = 5)
                                                 10^(-breaks)
                                               })) +
  AMOR::theme_blackbox()
p1

p1 <- ggplot(metawas, aes(x = af, y = p_lrt, col = snp_effect)) +
  facet_wrap(~ snp_effect) +
  geom_point() +
  scale_y_continuous(trans = scales::trans_new(name = 'minus_log',
                                               transform = function(x) -log10(x),
                                               inverse = function(x) 10^(-x),
                                               breaks = function(x){
                                                 rng <- range(x)
                                                 breaks <- seq(rng[1], rng[2], length.out = 5)
                                                 10^(-breaks)
                                               })) +
  AMOR::theme_blackbox()
p1


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
p1
p1$data


