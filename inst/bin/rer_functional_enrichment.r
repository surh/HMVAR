#!/usr/bin/env Rscript

# (C) Copyright 2020 Sur Herrera Paredes
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

library(tidyverse)
library(HMVAR)

# setwd("/cashew/users/sur/exp/fraserv/2020/today/gut")
args <- list(input = "rerperms/Bacteroides_ovatus_58035.rer.fdr.txt",
             annots = "annots/Bacteroides_ovatus_58035.emapper.annotations",
             outdir = "output/")


# Read data
RER <- read_tsv(args$input,
                col_types = cols(gene_id = col_character(),
                                 Rho = col_number(),
                                 P = col_number(),
                                 p.adj = col_number(),
                                 FDR = col_number()))
RER

annots <- read_eggnog(args$annots) %>%
  mutate(query_name = str_replace(string = query_name,
                                  pattern = "\\([+-]\\)_[0-9]+$",
                                  replacement = "")) %>%
  select(gene_id = query_name, predicted_gene_name,
         GO_terms, KEGG_KOs, BiGG_reactions,
         Annotation_tax_scope, OGs,
         'bestOG|evalue|score', COG_cat, eggNOG_annot)
annots

p1 <- ggplot(RER, aes(x = -log10(P), y = -log10(FDR))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  AMOR::theme_blackbox()
p1

pval_qqplot(RER$P)
pval_qqplot(RER$FDR)

p1 <- ggplot(RER, aes(x = Rho, y = -log10(FDR))) +
  geom_point(aes(size = N), shape = 19, alpha = 0.1) +
  scale_radius(trans="sqrt") +
  AMOR::theme_blackbox()
p1

# Match genes and annots
RER <- RER %>%
  filter(!is.na(FDR)) %>%
  left_join(annots, by = "gene_id")

RER %>%
  select(-p.adj,
         -BiGG_reactions,
         -Annotation_tax_scope,
         -"bestOG|evalue|score") %>%
  print(n = 30)

RER 
  


terms_enrichment(dat = RE)

