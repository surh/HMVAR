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
                       genes = NULL,
                       cds_only = FALSE)



# Match freqs and depth
Dat$depth <- Dat$depth %>% gather(key = "sample", value = 'depth', -site_id)
Dat$freq <- Dat$freq %>% gather(key = "sample", value = 'freq', -site_id)
Dat$info <- Dat$info %>% select(site_id, ref_id, ref_pos, major_allele, minor_allele)

# Find missing data
dat <- Dat$depth %>%
  inner_join(Dat$freq, by = c("site_id", "sample"))
dat$freq[ dat$depth < 1 ] <- NA
Dat$freq <- dat %>% select(-depth) %>% spread(sample, freq)

# Create BIMBAM tables
geno <- Dat$info %>%
  select(site_id, minor_allele, major_allele) %>%
  left_join(Dat$freq, by = "site_id")  

pheno <- map %>%
  filter(sample %in% colnames(geno)) %>%
  arrange(factor(sample, levels = colnames(geno)[-(1:3)])) %>%
  mutate(phenotype = 1*(Group == "Supragingival.plaque")) %>%
  select(id = sample, phenotype)

snp <- Dat$info %>% select(ID = site_id, pos = ref_pos, chr = ref_id)

# Write bimbam tables
dir.create(args$outdir)
dir.create(paste0(args$outdir, "/bimbam/"))
gen_file <- paste0(args$outdir, "/bimbam/geno.bimbam")
write_tsv(geno, path = gen_file, col_names = FALSE)

phen_file <- paste0(args$outdir, "/bimbam/pheno.bimbam")
write_tsv(pheno, path = phen_file)

snp_file <- paste0(args$outdir, "/bimbam/snp.bimbam")
write_tsv(snp, path = snp_file, col_names = FALSE)


