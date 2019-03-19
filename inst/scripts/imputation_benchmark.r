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
             m = 5)

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

### Benchmark imputation
Files$Dirs$data_hidden_dir <- file.path(args$outdir, "data_hidden_geno_files/")
Res <- benchmark_imputation(geno = midas_bimbam$Dat$geno,
                            snp = midas_bimbam$Dat$snp,
                            outdir = Files$Dirs$data_hidden_dir,
                            p = args$hidden_proportion,
                            m = args$m,
                            verbose = FALSE,
                            seed = args$seed)
Files$Files$imputed_geno_file <- Res$imputed_geno_file
