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
library(mice)

map <- read_tsv("map.txt", col_types = cols(.default = col_character())) %>%
  select(sample = ID, Group = Group)
map
Dat <- midas_to_bimbam(midas_dir = "midas_output_small/", outdir = "bimbam", map = map, focal_group = "Supragingival.plaque", prefix = "test")


imputed <- mice_impute(geno = Dat$Dat$geno,
                       snp = Dat$Dat$snp,
                       outdir = "imputed/",
                       m = 1,
                       verbose = FALSE,
                       prefix = "imputed",
                       return_table = TRUE,
                       seed = 234567)



imp <- Dat$Dat$geno %>% split(Dat$Dat$snp$chr) %>%
  map_dfr(~tidy_mice(.), m1 = 1, seed = 234567)

all.equal(imp, imputed$imp)




Dat$Dat$geno 

colnames(imp)[-1] == colnames(Dat$Dat$geno)[-(1:3)]
imp
Dat$Dat$geno %>% select(site_id, minor_allele)


table(is.na(Dat$Dat$geno[,-(1:3)]))
table(is.na(imp[,-1]))

summary(as.matrix(imp[,-1])[ is.na(Dat$Dat$geno[,-(1:3)]) & !is.na(imp[,-1]) ])
as.matrix(imp[,-1])[ is.na(Dat$Dat$geno[,-(1:3)]) & !is.na(imp[,-1]) ] %>% tibble(x = .) %>%
  ggplot(aes(x = x)) +
  geom_histogram(bins = 20) +
  AMOR::theme_blackbox()
