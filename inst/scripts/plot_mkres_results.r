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

args <- list(dir = "hmp.mktest.SPvsTD/",
             outdir = "hmp.mktest.SPvsTD.plots",
             threshold = 0.01,
             outfile = "hmp.mktest.SPvsTD_significant_genes.txt")

if(!dir.exists(args$outdir))
  dir.create(args$outdir)

files <- list.files(args$dir)

Res <- NULL
for(f in files){
  name <- str_replace(string = f, pattern = "_mktest.txt$", replacement = "")
  cat("\t", name, "\n")
  mkres <- process_mkres(f, args$dir)
  plot_mkres(mkres = mkres, name = name, dir = args$outdir, threshold = args$threshold)
  Res <- mkres %>% filter(p.value < args$threshold) %>% mutate(name = name) %>% bind_rows(Res)
}

write_tsv(Res, args$outfile)

# Repeat for the other comparison
args <- list(dir = "qin2012.mktest.DvsND/",
             outdir = "qin2012.mktest.DvsND.plots",
             threshold = 0.01,
             outfile = "qin2012.mktest.DvsND_significant_genes.txt")

if(!dir.exists(args$outdir))
  dir.create(args$outdir)

files <- list.files(args$dir)

Res <- NULL
for(f in files){
  name <- str_replace(string = f, pattern = "_mktest.txt$", replacement = "")
  cat("\t", name, "\n")
  mkres <- process_mkres(f, args$dir)
  plot_mkres(mkres = mkres, name = name, dir = args$outdir, threshold = args$threshold)
  Res <- mkres %>% filter(p.value < args$threshold) %>% mutate(name = name) %>% bind_rows(Res)
}

write_tsv(Res, args$outfile)
