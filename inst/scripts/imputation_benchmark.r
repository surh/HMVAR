library(HMVAR)
library(tidyverse)

indir <- commandArgs(trailingOnly = TRUE)[1]
spec <- commandArgs(trailingOnly = TRUE)[2]
indir <- "./"
spec <- "midas_output_small/"

# Eventually replace this with argparse
args <- list(midas_dir = file.path(indir, spec),
             map_file = "map.txt",
             outdir = "benchmark_imputation",
             prefix = "test",
             focal_group = "Supragingival.plaque",
             hidden_proportion = 0.1,
             seed = 76543,
             m = 1)

# Main output directory
dir.create(args$outdir)

# Create list for filenames
Files <- list(Dirs = list(),
              Files = list())

### Read and process original data
# Read map
map <- read_tsv(args$map_file, col_types = 'cc')
map <- map %>% select(sample = ID, Group = Group)

# Convert to bimbam
Files$Dirs$bimbam_dir <- file.path(args$outdir, "original")
midas_bimbam <- midas_to_bimbam(midas_dir = args$midas_dir,
                                map = map,
                                outdir = Files$Dirs$bimbam_dir,
                                focal_group = args$focal_group,
                                prefix = NULL)
Files$Files$midas_geno_file <- midas_bimbam$filenames$geno_file
Files$Files$pheno_file <- midas_bimbam$filenames$pheno_file
Files$Files$snp_file <- midas_bimbam$filenames$snp_file
rm(map)

### Hide 10% of data
Files$Dirs$data_hidden_dir <- file.path(args$outdir, "data_hidden_geno_files/")
dir.create(Files$Dirs$data_hidden_dir)

# Select positions to impute
gen_only <- midas_bimbam$Dat$geno %>% select(-site_id, -minor_allele, -major_allele)
set.seed(args$seed)
set.seed(12345)
res <- which(!is.na(gen_only), arr.ind = TRUE) %>%
  as_tibble() %>%
  bind_cols(hide = sample(c(0,1),
                      prob = c(1 - args$hidden_proportion, args$hidden_proportion),
                      size = nrow(.), replace = TRUE)) %>%
  filter(hide == 1)

# Collect observations
gen_only <- gen_only %>% as.matrix
ii <- res %>% select(row, col) %>% as.matrix
res$observed <- gen_only[ii]

# Hide data
gen_only[ii] <- NA
hidden_bimbam <- midas_bimbam$Dat$geno %>% select(site_id, minor_allele, major_allele) %>% bind_cols(gen_only %>% as_tibble)

# Impute
date()
imp <- mice_impute(geno = hidden_bimbam,
                   snp = midas_bimbam$Dat$snp,
                   outdir = file.path(args$outdir, "imputed_hidden/"),
                   m = args$m,
                   verbose = FALSE,
                   prefix = "imputed",
                   return_table = TRUE,
                   seed = args$seed)
date()

# Collect imputed values
if( any(imp$imp$site_id != hidden_bimbam$site_id)){
  stop("ERROR")
}
res$imputed <- (imp$imp %>% select(-site_id, -minor_allele, -major_allele) %>% as.matrix)[ii]
res$path <- args$outdir

# Compare
cor(res$observed, res$imputed, use = "complete.obs")
hist(res$imputed)
hist(res$observed)
p1 <- ggplot(res, aes(x = observed, y = imputed)) +
  geom_point() +
  geom_smooth(method = "lm")
p1
summary(res$imputed)

