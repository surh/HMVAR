library(HMVAR)
library(tidyverse)

benchmark_imputation <- function(geno, snp, outdir, p = 0.1 ,m = 5, verbose = FALSE, seed = NA){
  dir.create(outdir)
  
  # Select positions to impute
  # geno <- midas_bimbam$Dat$geno
  gen_only <- geno %>% dplyr::select(-site_id, -minor_allele, -major_allele)
  
  if (!is.na(seed)){
    set.seed(seed)
  }
  res <- dplyr::as_tibble(which(!is.na(gen_only), arr.ind = TRUE))
  res <- res %>%
    dplyr::bind_cols(hide = sample(c(0,1),
                                   prob = c(1 - p, p),
                                   size = nrow(.), replace = TRUE)) %>%
    filter(hide == 1)
  
  # Collect observations
  gen_only <- gen_only %>% as.matrix
  ii <- res %>%
    dplyr::select(row, col) %>%
    as.matrix
  res$observed <- gen_only[ii]
  
  # Hide data
  gen_only[ii] <- NA
  geno_hidden <- geno %>%
    dplyr::select(site_id, minor_allele, major_allele) %>%
    dplyr::bind_cols(dplyr::as_tibble(gen_only))
  
  # Impute
  t <- system.time(imp <- mice_impute(geno = geno_hidden,
                                      snp = snp,
                                      outdir = outdir,
                                      m = m,
                                      verbose = verbose,
                                      prefix = "imputed",
                                      return_table = TRUE,
                                      seed = seed))
  
  # Check imputation output
  if( any(imp$imp$site_id != geno_hidden$site_id)){
    stop("ERROR")
  }
  if( any(colnames(imp$imp) != colnames(geno_hidden))){
    stop("ERROR")
  }
  
  # Collect imputed values
  res$imputed <- (imp$imp %>%
                    select(-site_id, -minor_allele, -major_allele) %>%
                    as.matrix)[ii]
  res$path <- outdir
  
  # Calculate correlation
  r <- cor(res$observed, res$imputed, use = "complete.obs")
  p.imputed <- 1 - (is.na(res$imputed) / nrow(res))
  
  # Plot
  p1 <- ggplot(res, aes(x = observed, y = imputed)) +
    geom_point() +
    geom_smooth(method = "lm") +
    AMOR::theme_blackbox()
  filename <- filepath(outdir, "observed_vs_imputed.svg")
  ggsave(filename, p1, width = 5, height = 5)
  
  p1 <- res %>%
    gather(key = "Type", value = "allele_frequency", observed, imputed) %>%
    ggplot(aes(x=allele_frequency)) +
    facet_grid(Type ~ .) +
    geom_histogram(bins = 20) +
    AMOR::theme_blackbox()
  filename <- filepath(outdir, "alllele_freq_histograms.svg")
  ggsave(filename, p1, width = 12, height = 5)
  
  return(list(r = r, p.imputed = p.imputed, res = res, imputed_geno_file = imp$imputed_file))
}

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









# 
# # Compare
# 
# hist(res$imputed)
# hist(res$observed)
# p1 <- ggplot(res, aes(x = observed, y = imputed)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# p1

