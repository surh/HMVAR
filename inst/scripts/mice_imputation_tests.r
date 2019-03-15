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


myfun <- function(d, m = 5, verbose = FALSE){
  d %>%
    select(-site_id, -minor_allele, -major_allele) %>%
    filter_all(any_vars(!is.na(.))) %>%
    mice(m = 5, printFlag = verbose) %>% 
    mice::complete() %>%
    as_tibble
}

imp <- Dat$Dat$geno %>% split(Dat$Dat$snp$chr) %>%
  map_dfr(~myfun(.), .id = "site_id", m1 = 1)
imp

table(is.na(Dat$Dat$geno[,-(1:3)]))
table(is.na(imp[,-1]))

summary(as.matrix(imp[,-1])[ is.na(Dat$Dat$geno[,-(1:3)]) & !is.na(imp[,-1]) ])
as.matrix(imp[,-1])[ is.na(Dat$Dat$geno[,-(1:3)]) & !is.na(imp[,-1]) ] %>% tibble(x = .) %>%
  ggplot(aes(x = x)) +
  geom_histogram(bins = 20) +
  AMOR::theme_blackbox()
