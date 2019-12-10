library(HMVAR)
library(tidyverse)


geno <- read_tsv("Neisseria_flavescens_61757/geno.bimbam",
                 col_names = FALSE,
                 col_types = cols(X1 = col_character(),
                                  X2 = col_character(),
                                  X3 = col_character(),
                                  .default = col_double()))
colnames(geno)[1:3] <- c("site_id", "minor_allele", "major_allele")
geno

snp <- read_tsv("Neisseria_flavescens_61757/snp.bimbam",
                col_names = FALSE,
                col_types = cols(X1 = col_character(),
                                 X2 = col_number(),
                                 X3 = col_character()))
colnames(snp) <- c("ID", "ps", "chr")
snp

seed <- 76543
verbose <- TRUE
block_size <- 100
m <- 5
p <- 0.1
outdir <- "test"

### benchmark
# dir.create(outdir)

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

# # Impute
# t <- system.time(imp <- mice_impute(geno = geno_hidden,
#                                     snp = snp,
#                                     outdir = outdir,
#                                     m = m,
#                                     verbose = verbose,
#                                     prefix = "imputed",
#                                     return_table = TRUE,
#                                     seed = seed,
#                                     block_size = block_size))

Res <- res
rm(res)


######### impute
imp <- geno_hidden %>%
  split(snp$chr) 
# %>%
#   purrr::map_dfr(~tidy_mice(.), m1 = m, verbose = verbose, seed = seed, block_size = block_size) %>%
#   dplyr::arrange(match(site_id, snp$ID))
d <- imp[[2]]
d

##### tidy
res <- d %>%
  dplyr::select(site_id, minor_allele, major_allele)
res

# Remove sites that are completeley missing
ii <- !(d %>% dplyr::select(-site_id, -minor_allele, -major_allele) %>% is.na %>% apply(1, all))
d <- d %>% dplyr::filter(ii)
d

if(sum(ii) == 0){
  return(res %>% dplyr::left_join(d, by = "site_id"))
}

if(block_size > 0){
  pred <- sapply(1:nrow(d), HMVAR:::create_blocks, block_size = block_size, total_pos = nrow(d))
  diag(pred) <- 0
}else{
  pred <- matrix(1, nrow = nrow(d), ncol = nrow(10))
  diag(pred) <- 0
}
# image(pred)

ids <- d$site_id
samples <- colnames(d)[-(1:3)]
d <- d %>%
  dplyr::select(-site_id, -minor_allele, -major_allele) %>%
  t %>% mice::mice(m = m, printFlag = verbose, seed = seed, method = "pmm", predictorMatrix = pred)

# plot(d)
# mice:::stripplot.mids(d)
d <- d %>%
  mice::complete() %>%
  t
colnames(d) <- samples
d <- bind_cols(tibble(site_id = ids), as_tibble(d))

res <- res %>% dplyr::left_join(d, by = "site_id")
#####



