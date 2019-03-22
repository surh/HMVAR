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

library(HMVAR)
library(tidyverse)
source("functions.r")

# args <- list(dir = "hmp.mktest.SPvsTD/",
#              threshold = 0.01,
#              ni_file = "hmp.SPvsTD.NI.txt",
#              sig_genes = "hmp.mktest.SPvsTD_significant_genes.txt",
#              outdir = "hmp.all")

args <- list(dir = "qin2012.mktest.DvsND//",
             threshold = 0.01,
             ni_file = "qin2012.DvsND.NI.txt",
             sig_genes = "qin2012.mktest.DvsND_significant_genes.txt",
             outdir = "qin2012.all")

if(!dir.exists(args$outdir))
  dir.create(args$outdir)


files <- list.files(args$dir)

Res <- NULL
for(f in files){
  name <- str_replace(string = f, pattern = "_mktest.txt$", replacement = "")
  cat("\t", name, "\n")
  mkres <- process_mkres(f, args$dir)
  Res <- mkres %>% mutate(name = name) %>% bind_rows(Res)
}

# Get NI and alpha
NI <- read_tsv(args$ni_file)
NI <- NI %>% mutate(name = str_replace(string = File, pattern = "_mktest.txt$", replacement = "")) %>%
  arrange(desc(alpha))
NI %>% head(20)

# Match NI with number of genes
gene_numbers <- Res %>% split(.$name) %>% map_int(nrow)
gene_numbers <- tibble(name = names(gene_numbers), ngenes = gene_numbers)
NI <- NI %>% inner_join(gene_numbers, by = "name")
NI %>% head(20)

# Read sig genes
Sig <- read_tsv(args$sig_genes)
Sig

# Match NI with signigicant genes
gene_numbers <- Sig %>% split(.$name) %>% map_int(nrow)
gene_numbers <- tibble(name = names(gene_numbers), sig.genes = gene_numbers)
NI <- NI %>% left_join(gene_numbers, by = "name")
NI %>% head(20)

# Correct for multiple testing by genome
Res <- Res %>% split(.$name) %>%
  map_dfr(~ mutate(., q.value = p.adjust(p.value, 'fdr')))

##### Start plotting

### Number of tested genes vs alpha
p1 <- ggplot(NI, aes(x = alpha, y = ngenes)) +
  geom_point(aes(color = alpha > 0)) +
  geom_smooth()
p1
filename <- paste0(args$outdir, "/", basename(args$outdir), "_testedVSalpha.points.svg")
ggsave(filename, p1, width = 6, height = 4)


p1 <- ggplot(NI, aes(x = alpha > 0, y = ngenes, color = alpha > 0)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme_blackbox()
p1
filename <- paste0(args$outdir, "/", basename(args$outdir), "_testedVSalpha.boxplot.svg")
ggsave(filename, p1, width = 6, height = 4)


### Number of sig.genes vs alpha
p1 <- ggplot(NI, aes(x = alpha, y = sig.genes)) +
  geom_point(aes(color = alpha > 0)) +
  geom_smooth()
p1
filename <- paste0(args$outdir, "/", basename(args$outdir), "_sigVSalpha.points.svg")
ggsave(filename, p1, width = 6, height = 4)

p1 <- ggplot(NI, aes(x = alpha > 0, y = sig.genes, color = alpha > 0)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme_blackbox()
p1
filename <- paste0(args$outdir, "/", basename(args$outdir), "_sigVSalpha.boxplot.svg")
ggsave(filename, p1, width = 6, height = 4)

### significant genes vs ngenes
p1 <- ggplot(NI, aes(x = ngenes, y = sig.genes)) +
  geom_point(aes(color = alpha > 0)) +
  geom_smooth()
p1
filename <- paste0(args$outdir, "/", basename(args$outdir), "_genesVSsig.svg")
ggsave(filename, p1, width = 6, height = 4)

# ### pvalue qqplot
# pval_qqplot(Res$p.value)
plot_mkres(mkres = Res, name = basename(args$outdir), dir = args$outdir, threshold = args$threshold)

# 
# pos_genomes <- filter(NI, alpha > 0)$name
# 
# ### Only genomes with positive alpha
# Res.pos <- Res %>% filter(name %in% pos_genomes)
# 
# pval_qqplot(Res.pos$p.value)
# 
# summary(qvalue::qvalue(Res$p.value))

# pval_qqplot(Res$p.value[Res$p.value < 1])
# pval_qqplot(Res.pos$p.value[Res.pos$p.value < 1])
# summary(qvalue::qvalue(Res.pos$p.value[Res.pos$p.value < 1]))
