library(HMVAR)
library(tidyverse)

indir <- "./"
spec <- "midas_output_small/"

# Eventually replace this with argparse
args <- list(midas_dir = file.path(indir, spec),
             map_file = "map.txt",
             outdir = "benchmark_imputation",
             hidden_proportion = 0.1,
             seed = 76543,
             m = 5)
args <- process_arguments()

# Main output directory
cat("Creating output directory...\n")
dir.create(args$outdir)

# Create list for filenames
Files <- list(Dirs = list(),
              Files = list())

# Read map
cat("Reading map...\n")
map <- read_tsv(args$map_file, col_types = 'cc')
map <- map %>% select(sample = ID, Group = Group)

# Convert to bimbam
cat("Converting MIDAS data to BIMBAM...\n")
Files$Dirs$bimbam_dir <- file.path(args$outdir, "original")
midas_bimbam <- midas_to_bimbam(midas_dir = args$midas_dir,
                                map = map,
                                outdir = Files$Dirs$bimbam_dir,
                                focal_group = "",
                                prefix = NULL)
Files$Files$midas_geno_file <- midas_bimbam$filenames$geno_file
Files$Files$pheno_file <- midas_bimbam$filenames$pheno_file
Files$Files$snp_file <- midas_bimbam$filenames$snp_file
rm(map)


### Benchmark imputation
cat("Benchmarking imputation...\n")
Files$Dirs$data_hidden_dir <- file.path(args$outdir, "data_hidden_geno_files/")
Res <- benchmark_imputation(geno = midas_bimbam$Dat$geno,
                            snp = midas_bimbam$Dat$snp,
                            outdir = Files$Dirs$data_hidden_dir,
                            p = args$hidden_proportion,
                            m = args$m,
                            verbose = TRUE,
                            seed = args$seed)
Files$Files$imputed_geno_file <- Res$imputed_geno_file


#' Imputation via mice
#' 
#' Imputes missing genotypes using mice
#'
#' @param geno A data table with genotype information. First three columns
#' must be site_id, minor_allele and major allele; folowed by one column per
#' sample with the sample ID as column name. Missing genotypes must be encoded
#' as NA values
#' @param snp SNP information data table. Must have columns ID, pos and chr.
#' ID mut correspond to the "site_id" column in \code{geno}, and the SNPs must
#' be in the same order in both tables.
#' @param outdir Output directory to create and write the imputed table in
#' BIMBAM format.
#' @param m Number of multiple imputations. See \link{mice} documentation.
#' @param verbose Whether to print progress during imputation.
#' @param prefix Prefix for imputed file.
#' @param return_table If TRUE, the function will return a list with both
#' the output file path and the imputed table. If FALSE, only the file path
#' will be returned.
#' @param seed Seed for imputation. See \link{mice} documentation.
#'
#' @return Either the output file path, or a list containing the file path
#' and the imputed table.
#' 
#' @export
#' 
#' @importFrom magrittr %>%
mice_impute <- function(geno, snp,
                        outdir = "imputed/",
                        m = 5,
                        verbose = FALSE,
                        prefix = "imputed",
                        return_table = FALSE,
                        seed = NA){
  
  geno <- midas_bimbam$Dat$geno
  snp <- midas_bimbam$Dat$snp
  m <- 5
  outdir <- 'imputed/'
  verbose <- FALSE
  prefix <- 'imputed'
  return_table <- TRUE
  seed <- 12345
  
  if(any(geno$site_id != snp$ID)){
    stop("ERROR: geno and snp tables do not match", call. = TRUE)
  }
  
  # imp <- geno %>%
  #   split(snp$chr)
  # d <- imp[[1]]
  # d
  
  imp <- geno %>%
    split(snp$chr) %>%
    purrr::map_dfr(~tidy_mice(.), m1 = m, verbose = verbose, seed = seed) %>%
    dplyr::arrange(match(site_id, snp$ID))
  
  # Write results
  dir.create(outdir)
  
  filename <- file.path(outdir,
                        paste(c(prefix, "mean.genotype.txt"), collapse = "."))
  readr::write_tsv(imp, path = filename, col_names = FALSE)
  
  if(return_table){
    return(list(imputed_file = filename, imp = imp))
  }else{
    return(filename)
  }
}

#' Tidy mice
#' 
#' Internal utility function
#' 
#' Calls mice on data table
#'
#' @param d 
#' @param m 
#' @param verbose 
#'
#' @return A tibble with imputed results
#' 
#' @importFrom magrittr %>%
tidy_mice <- function(d, m = 5, verbose = FALSE, seed = NA){
  res <- d %>%
    dplyr::select(site_id, minor_allele, major_allele)
  
  # Remove sites that are completeley missing
  ii <- !(d %>% dplyr::select(-site_id, -minor_allele, -major_allele) %>% is.na %>% apply(1, all))
  d <- d %>% dplyr::filter(ii)
  
  if(sum(ii) == 0){
    return(res %>% dplyr::left_join(d, by = "site_id"))
  }
  
  ids <- d$site_id
  samples <- colnames(d)[-(1:3)]
  d <- d %>%
    dplyr::select(-site_id, -minor_allele, -major_allele) %>%
    t %>% mice::mice(m = 1, printFlag = TRUE, seed = seed, method = "pmm", maxit = 1) %>%
    mice::complete() %>%
    t
  colnames(d) <- samples
  d <- bind_cols(tibble(site_id = ids), as_tibble(d))
  # d <- d %>%
  #   dplyr::select(-minor_allele, -major_allele) %>%
  #   mice::mice(m = m, printFlag = verbose, seed = seed) %>% 
  #   mice::complete() %>%
  #   dplyr::as_tibble()
  
  res %>% dplyr::left_join(d, by = "site_id")
}


create_blocks <- function(i, block_size=3, total_pos=10){
  start <- i - ceiling(block_size / 2) + 1
  end <- start + block_size - 1
  start <- max(start, 1)
  end <- min(total_pos, end)
  
  v <- rep(0, length.out = total_pos)
  v[start:end] <- 1
  
  v
}



# res <- d %>%
#   dplyr::select(site_id, minor_allele, major_allele)
# 
# # Remove sites that are completeley missing
# ii <- !(d %>% dplyr::select(-site_id, -minor_allele, -major_allele) %>% is.na %>% apply(1, all))
# d <- d %>% dplyr::filter(ii)
# 
# if(sum(ii) == 0){
#   return(res %>% dplyr::left_join(d, by = "site_id"))
# }
# 
# ids <- d$site_id
# samples <- colnames(d)[-(1:3)]
# d <- d %>%
#   dplyr::select(-site_id, -minor_allele, -major_allele) %>%
#   t
# # colnames(d) <- res$site_id
# d <- d %>% mice::mice(m = 1, printFlag = TRUE, seed = seed, method = "pmm", maxit = 1) 
# d <- t(mice::complete(d))
# colnames(d) <- samples
# d <- bind_cols(tibble(site_id = ids), as_tibble(d))
# 
# res <- res %>% dplyr::left_join(d, by = "site_id")
# res





d <- d %>%
  dplyr::select(-minor_allele, -major_allele)
mice::mice(m = m, printFlag = verbose, seed = seed) %>% 
  mice::complete() %>%
  dplyr::as_tibble()

res %>% dplyr::left_join(d, by = "site_id")
