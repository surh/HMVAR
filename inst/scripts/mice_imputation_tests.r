library(HMVAR)
library(tidyverse)
library(mice)


# data("airquality")
# airquality
# imp <- mice(airquality, m = 1)
# airquality_impute <- complete(imp)
# imp

map <- read_tsv("map.txt", col_types = cols(.default = col_character())) %>%
  select(sample = ID, Group = Group)
map
Dat <- midas_to_bimbam(midas_dir = "midas_output_small/", outdir = "bimbam", map = map, focal_group = "Supragingival.plaque", prefix = "test")


myfun <- function(d, m = 5){
  d %>%
    select(-site_id, -minor_allele, -major_allele) %>%
    filter_all(any_vars(!is.na(.))) %>%
    mice(m = 5) %>% 
    mice::complete() %>%
    as_tibble
}

imp <- Dat$Dat$geno %>% split(Dat$Dat$snp$chr) %>%
  map_dfr(~myfun(.), .id = "site_id", m1 = 1)
imp


i <- d[[5]] %>%
  select(-site_id, -minor_allele, -major_allele) %>%
  filter_all(any_vars(!is.na(.))) %>%
  mice(m = 5) %>% 
  mice::complete()

table(is.na(d[[5]][,-(1:3)]))
table(is.na(i))




i[is.na(d[[5]][,-(1:3)]) & !is.na(i)]
as.matrix(d[[5]][,-(1:3)])[is.na(d[[5]][,-(1:3)]) & !is.na(i)]



rowSums(!is.na(d[[5]]))
colSums(!is.na(d[[5]]))

Dat$Dat$geno %>%
  select(-site_id, -minor_allele, -major_allele) %>%
  filter_all(any_vars(!is.na(.)))

tibble(x = c(NA,1,NA), y = c(1,1,NA)) %>% filter_all(any_vars(!is.na(.)))
tibble(x = c(0,1,0), y = c(1,1,0)) %>% filter_all(all_vars(!is.na(.)))
