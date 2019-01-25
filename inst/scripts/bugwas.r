library(HMVAR)
library(tidyverse)
library(bugwas)

Sys.setenv(LD_LIBRARY_PATH="/opt/modules/pkgs/eqtlbma/git/lib/")


args <- list(midas_dir = "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/Actinomyces_odontolyticus_57475/",
             map_file = "hmp_SPvsTD_map.txt",
             outdir = "testout/")


map <- read_tsv(args$map_file, col_types = 'cc')
map <- map %>% select(sample = ID, Group = Group)
Dat <- read_midas_data(midas_dir = args$midas_dir,
                       map = map,
                       genes = NULL)



# Match freqs and depth
Dat$depth <- Dat$depth %>% gather(key = "sample", value = 'depth', -site_id)
Dat$freq <- Dat$freq %>% gather(key = "sample", value = 'freq', -site_id)
Dat$info <- Dat$info %>% select(site_id, ref_id, major_allele, minor_allele)

# Find missing data
dat <- Dat$depth %>%
  inner_join(Dat$freq, by = c("site_id", "sample"))
dat$freq[ dat$depth < 1 ] <- NA
Dat$freq <- dat %>% select(-depth) %>% spread(sample, freq)

geno <- Dat$info %>%
  select(site_id, minor_allele, major_allele) %>%
  left_join(Dat$freq, by = "site_id")  

pheno <- map %>%
  filter(sample %in% colnames(geno)) %>%
  arrange(factor(sample, levels = colnames(geno)[-(1:3)])) %>%
  mutate(phenotype = 1*(Group == "Supragingival.plaque")) %>%
  select(id = sample, phenotype)




%>%
  left_join(map, by = "sample") %>%
  filter(depth >= depth_thres) %>%
  left_join(meta, by = "site_id") %>%
  mutate(allele = replace(freq, freq < 0.5, 'major')) %>%
  mutate(allele = replace(allele, freq >= 0.5, 'minor')) %>%
  filter(distribution != "Invariant")
dat

dat <- dat %>%
  mutate(site_id = factor(site_id,
                          levels = as.character(unique(sort(as.numeric(dat$site_id))))))
