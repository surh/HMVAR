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


# args <- list(dir = "qin2012.mktest.DvsND//",
#              threshold = 0.01,
#              ni_file = "qin2012.DvsND.NI.txt",
#              sig_genes = "qin2012.mktest.DvsND_significant_genes.txt",
#              outdir = "qin2012.all")


args <- list(vmwa.genes.dir = "hmp.vmwa.pvals.genes/",
             mkres.dir = "../2018-12-14.hmp_mktest/results/",
             q.value.thres = 0.01)

# List of files from mktest
mkres_files <- list.files(args$mkres.dir)

Res <- NULL
for(f in mkres_files){
  # f <- mkres_files[1]
  species <- str_replace(string = f, pattern = "_mktest.txt$", replacement = "")
  vmwa_file <- paste0(species, "_vmwa.genes.txt")
  cat(species, "\n", "\t", f, "\n\t", vmwa_file, "\n")
  
  # Read data
  vmwa <- read_tsv(paste0(args$vmwa.genes.dir, "/", vmwa_file))
  vmwa <- vmwa %>% select(everything(),vmwa.OR = OR, vmwa.p.value = p.value, vmwa.q.value = q.value)
  mkres <- read_tsv(paste0(args$mkres.dir, "/", f))
  mkres <- mkres %>% select(gene_id = gene, everything()) %>% select(-ratio.pval)
  
  # MK test
  mkpval <- mkres %>%
    pmap_dbl(.f = function(Dn, Ds, Pn, Ps, ...){
      mat <- matrix(c(Dn, Ds, Pn, Ps), ncol = 2)
      res <- fisher.test(mat)
      return(res$p.value)})
  mkres <- mkres %>% add_column(mk.p.value = mkpval)
  
  # Join
  Full <- vmwa %>% full_join(mkres, by = "gene_id") %>% select(-ref_id, -start, -end)
  Full <- Full %>%
    add_column(significant = 1*(Full$vmwa.p.value < 0.01) + 2*(Full$mk.p.value < 0.01)) %>%
    mutate(significant = factor(significant)) %>%
    mutate(significant = recode_factor(significant,
                                       `0` = "none",
                                       `1` = "vmwa",
                                       `2` = "mktest",
                                       `3` = "both"))
  
  # p1 <- ggplot(Full, aes(x = vmwa.OR, y = ratio)) +
  #   geom_point(aes(col = significant)) + 
  #   scale_y_log10() + 
  #   scale_x_log10() +
  #   AMOR::theme_blackbox()
  # p1
  
  # Select
  Full <- Full %>%
    filter(vmwa.OR > 1 & vmwa.p.value < 0.01 & !is.na(ratio)) %>%
    add_column(species = species)
  
  Res <- Res %>% bind_rows(Full)
}

write_tsv(Res, "hmp.vmwaselected.mktest.txt")

p1 <- ggplot(Res, aes(x = vmwa.OR, y = ratio)) +
  geom_point(aes(col = significant)) +
  # scale_y_log10() +
  # scale_x_log10() +
  AMOR::theme_blackbox()
p1


#####################


args <- list(vmwa.genes.dir = "qin2012.vmwa.pvals.genes/",
             mkres.dir = "../2018-12-14.qin2012_mktest/results/",
             q.value.thres = 0.01)

# List of files from mktest
mkres_files <- list.files(args$mkres.dir)

Res <- NULL
for(f in mkres_files){
  # f <- mkres_files[67]
  species <- str_replace(string = f, pattern = "_mktest.txt$", replacement = "")
  vmwa_file <- paste0(species, "_vmwa.genes.txt")
  cat(species, "\n", "\t", f, "\n\t", vmwa_file, "\n")
  
  # Read data
  vmwa <- read_tsv(paste0(args$vmwa.genes.dir, "/", vmwa_file))
  vmwa <- vmwa %>% select(everything(),vmwa.OR = OR, vmwa.p.value = p.value, vmwa.q.value = q.value)
  mkres <- read_tsv(paste0(args$mkres.dir, "/", f))
  mkres <- mkres %>% select(gene_id = gene, everything()) %>% select(-ratio.pval)
  
  # MK test
  mkpval <- mkres %>%
    pmap_dbl(.f = function(Dn, Ds, Pn, Ps, ...){
      mat <- matrix(c(Dn, Ds, Pn, Ps), ncol = 2)
      res <- fisher.test(mat)
      return(res$p.value)})
  mkres <- mkres %>% add_column(mk.p.value = mkpval)
  
  # Join
  Full <- vmwa %>% full_join(mkres, by = "gene_id") %>% select(-ref_id, -start, -end)
  Full <- Full %>%
    add_column(significant = 1*(Full$vmwa.p.value < 0.01) + 2*(Full$mk.p.value < 0.01)) %>%
    mutate(significant = factor(significant)) %>%
    mutate(significant = recode_factor(significant,
                                       `0` = "none",
                                       `1` = "vmwa",
                                       `2` = "mktest",
                                       `3` = "both"))
  
  # p1 <- ggplot(Full, aes(x = vmwa.OR, y = ratio)) +
  #   geom_point(aes(col = significant)) + 
  #   scale_y_log10() + 
  #   scale_x_log10() +
  #   AMOR::theme_blackbox()
  # p1
  
  # Select
  Full <- Full %>%
    filter(vmwa.OR > 1 & vmwa.p.value < 0.01 & !is.na(ratio)) %>%
    add_column(species = species)
  
  Res <- Res %>% bind_rows(Full)
}

write_tsv(Res, "hmp.vmwaselected.mktest.txt")
write_tsv(Res, "qin2012.vmwaselected.mktest.txt")

p1 <- ggplot(Res, aes(x = vmwa.OR, y = ratio)) +
  geom_point(aes(col = significant)) +
  # scale_y_log10() +
  # scale_x_log10() +
  AMOR::theme_blackbox()
p1


