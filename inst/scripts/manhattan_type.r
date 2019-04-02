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


library(tidyverse)

indir <- "lmmres/"
outdir <- "manhattans/"
dir.create(outdir)

files <- list.files(indir)
for(lmm_file in files){
  # lmm_file <- "Actinomyces_odontolyticus_57475_lmm.results.txt"
  genome <- basename(lmm_file) %>% str_replace("_lmm.results.txt$", "")
  cat(genome, "\n")
  # lmm <- read_tsv(file.path(indir,lmm_file),
  #                 col_types = 'ccnnnnnnnc')

  lmm <- read_tsv(file.path(indir,lmm_file),
		  col_types = cols(chr = col_character(),
				   rs = col_character(),
				   allele1 = col_character(),
				   allele0 = col_character(),
				   type = col_character(),
				   .default = col_double()))

  # lmm
  
  p1 <- ggplot(lmm,
               aes(x = ps, y = -log10(p_lrt.lmm))) +
    facet_grid(~chr, scales = "free_x", space = "free_x") +
    geom_point(aes(color = type)) +
    geom_hline(yintercept = -log10(1e-6), color = "red", size = 2) +
    theme_classic()
  # print(p1)
  filename <- paste0(outdir, "/", genome, ".manhattan.lmm.png")
  ggsave(filename, p1, width = 16, height = 4)
  
  
  p1 <- ggplot(lmm,
               aes(x = ps, y = -log10(p_lrt.lmmpcs))) +
    facet_grid(~chr, scales = "free_x", space = "free_x") +
    geom_point(aes(color = type)) +
    geom_hline(yintercept = -log10(1e-6), color = "red", size = 2) +
    theme_classic()
  # print(p1)
  filename <- paste0(outdir, "/", genome, ".manhattan.lmmpcs.png")
  ggsave(filename, p1, width = 16, height = 4)
}
