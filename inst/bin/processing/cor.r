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




args <- list(vmwa = "/godot/users/sur/exp/fraserv/2018/2018-12-17.plot_genes/hmp.vmwa.pvals.genes/Granulicatella_adiacens_61980_vmwa.genes.txt",
             mktest = "/godot/users/sur/exp/fraserv/2018/2018-12-14.hmp_mktest/results/Granulicatella_adiacens_61980_mktest.txt",
             prefix = "Granulicatella_adiacens_61980",
             stat1 = 'OR',
             stat2 = 'ratio',
             outdir = "out/",
             log10 = TRUE,
             stat_thres = 2)


# Read data
dat1 <- read_tsv(args$vmwa)
dat1
dat2 <- read_tsv(args$mktest)
dat2    

# Match
dat2 <- dat2 %>%
  select(gene_id = gene, everything())
dat <- dat1 %>% inner_join(dat2, by = "gene_id")
dat

# Filter
ii <- is.na(dat[,args$stat1,drop = TRUE]) | is.na(dat[,args$stat2,drop = TRUE])
dat <- dat %>%
  filter(!ii)

# Check "significance"
res <- 1*(dat[,args$stat1] > args$stat_thres) + 2*(dat[,args$stat2] > args$stat_thres)
dat <- dat %>%
  add_column(significant = as.vector(res)) %>%
  mutate(significant = factor(significant)) %>%
  mutate(significant = recode_factor(significant,
                                     `0` = "none",
                                     `1` = "vmwa",
                                     `2` = "mktest",
                                     `3` = "both"))

p1 <- ggplot(dat, aes_string(x = args$stat1, y = args$stat2)) +
  geom_point(aes(col = significant)) +
  # scale_y_log10() +
  # scale_x_log10() +
  AMOR::theme_blackbox()
p1

p1 <- ggplot(Full,
             aes(x = -log10(vmwa.p.value),
                 y = -log10(mk.p.value))) +
  geom_point(aes(col = significant)) +
  # scale_y_log10() +
  # scale_x_log10() +
  AMOR::theme_blackbox()
p1
ggsave("wvma_mk_cor.png", p1, width = 5, height = 5, dpi = 150)
