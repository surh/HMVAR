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
