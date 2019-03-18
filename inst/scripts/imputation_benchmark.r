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
             hidden_proportion = 0.1)

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
set.seed(12345)
positions <- is.na(midas_bimbam$Dat$geno) %>% which(arr.ind = TRUE) %>%
  as_tibble() %>%
  filter(col > 3) %>%
  bind_cols(hide = sample(c(0,1),
                          prob = c(1 - args$hidden_proportion, args$hidden_proportion),
                          size = length(.$col), replace = TRUE))
positions %>% filter(hide == 1)

geno_partial <- midas_bimbam$Dat$geno

geno_partial

gen_only <- midas_bimbam$Dat$geno %>% select(-site_id, -minor_allele, -major_allele)

res <- which(!is.na(gen_only), arr.ind = TRUE) %>%
  as_tibble() %>%
  bind_cols(hide = sample(c(0,1),
                      prob = c(1 - args$hidden_proportion, args$hidden_proportion),
                      size = nrow(.), replace = TRUE)) %>%
  filter(hide == 1)
res
gen_only
head(res)
dim(res)
sum(res$hide)

gen_only <- gen_only %>% as.matrix
ii <- res %>% select(row, col) %>% as.matrix
res$observed <- gen_only[ii]
# Hide data
gen_only[ii] <- NA
hidden_bimbam <- midas_bimbam$Dat$geno %>% select(site_id, minor_allele, major_allele) %>% bind_cols(gen_only %>% as_tibble)


imp <- mice_impute(geno = hidden_bimbam,
                   snp = midas_bimbam$Dat$snp,
                   outdir = file.path(args$outdir, "imputed_hidden/"),
                   m = 1,
                   verbose = FALSE,
                   prefix = "imputed",
                   return_table = TRUE,
                   seed = 76543)

# Collect imputed values
if( any(imp$imp$site_id != hidden_bimbam$site_id)){
  stop("ERROR")
}
res$imputed <- (imp$imp %>% select(-site_id, -minor_allele, -major_allele) %>% as.matrix)[ii]
res
# 
# geno_partial[c(1,3,4,5),c(1,3,4,5)]
# as.matrix(geno_partial)[c(1,2,3,4,5), c(2,3,4,5)]
# 
# matrix(c(1,2,3,4,5,1,2,3,4,5), ncol = 2)
# 
# geno_partial[matrix(c(1,2,3,4,5,1,2,3,4,5), ncol = 2)]
# as.matrix(geno_partial)[matrix(c(1,2,3,4,5,1,2,3,4,5), ncol = 2)]
# 
# gen_only <- midas_bimbam$Dat$geno %>% select(-site_id, -minor_allele, -major_allele) %>% as.matrix
# # gen_only <- gen_only[ rowSums(is.na(gen_only)) < 20, colSums(is.na(gen_only)) < 500 ]
# gen_only_mask <- !is.na(gen_only)
# samples_per_snp <- rowSums(gen_only_mask)
# snps_per_sample <- colSums(gen_only_mask)
# gen_only <- gen_only[ samples_per_snp > min_samples, snps_per_sample > min_snps ]
# # gen_only
# # image(gen_only)
# 
# gen_only_mask <- which(!is.na(gen_only))
# length(gen_only_mask)
# hidden_index <- sample(gen_only_mask,
#                        size = round(length(gen_only_mask) * args$hidden_proportion,
#                                     digits = 0), replace = FALSE)
# batch_name <- "batch0"
# res <- tibble(batch = batch_name, gen = gen_only[hidden_index])
# 
# # sum(!is.na(gen_only))
# gen_only[hidden_index] <- NA
# 
# # sum(!is.na(gen_only))
# gen_only <- cbind(midas_bimbam$Dat$geno %>%
#                     filter(samples_per_snp > min_samples) %>%
#                     select(site_id, minor_allele, major_allele),
#                   gen_only) %>%
#   as_tibble
# # gen_only <- cbind(midas_bimbam$Dat$geno %>% select(site_id, minor_allele, major_allele), gen_only) %>% as_tibble
# # gen_only
# # Write batch with hidden data
# hidden_gen_file <- file.path(Files$Dirs$data_hidden_dir, paste(c(batch_name, 'geno.bimbam'), collapse = "_"))
# write_tsv(gen_only, path = hidden_gen_file, col_names = FALSE)
# hidden_pheno_file <- file.path(Files$Dirs$data_hidden_dir, paste(c(batch_name, 'pheno.bimbam'), collapse = "_"))
# write_tsv(midas_bimbam$Dat$pheno %>% filter(snps_per_sample > min_snps) %>% select(phenotype),
#           path = hidden_pheno_file, col_names = FALSE)
# 
# # Impute batch
# Files$Dirs$imputed_dir <- file.path(args$outdir, "imputed_from_hidden")
# hidden_immputed_file <- bimbam_impute(geno_file = hidden_gen_file,
#                                       pheno_file = hidden_pheno_file,
#                                       pos_file = midas_bimbam$filenames$snp_file,
#                                       bimbam = args$bimbam,
#                                       outdir = Files$Dirs$imputed_dir,
#                                       em_runs = 10,
#                                       em_steps = 20,
#                                       em_clusters = 15,
#                                       prefix = batch_name)
# # Load imputed and compare
# hidden_imputed <- read_delim(hidden_immputed_file, col_names = FALSE, delim = " ",
#                            col_types = cols(X1 = col_character(),
#                                             X2 = col_character(),
#                                             X3 = col_character(),
#                                             .default = col_number()))
# # hidden_imputed
# # match(c('a','a', 'c'), letters[1:5] )
# # midas_bimbam$Dat$geno
# # Sorting
# hidden_imputed <- hidden_imputed[ match(midas_bimbam$Dat$geno$site_id, hidden_imputed$X1 ), ] %>% drop_na
# hidden_imputed <- hidden_imputed %>% select(-X1, -X2, -X3) %>% as.matrix
# res <- res %>% bind_cols(imputed = hidden_imputed[hidden_index]) %>% mutate(dir = args$midas_dir)
# write_tsv(res, file.path(args$outdir, "imputation_results.txt"))


cor(res$observed, res$imputed, use = "complete.obs")
hist(res$imputed)
hist(res$observed)
p1 <- ggplot(res, aes(x = observed, y = imputed)) +
  geom_point() +
  geom_smooth(method = "lm")
p1
summary(res$imputed)
# image(hidden_imputed)
# image(gen_only[,-(1:3)] %>% as.matrix())
# 
# gen_only
# x <- gen_only %>% select(-site_id, -minor_allele, -major_allele) %>% as.matrix %>% as.vector
# y <- hidden_imputed %>% as.vector
# 
# summary(x[!is.na(x)] - y[!is.na(x)])
