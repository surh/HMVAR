# (C) Copyright 2018-2019 Sur Herrera Paredes

# This file is part of HMVAR.
# 
# HMVAR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# HMVAR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with HMVAR.  If not, see <http://www.gnu.org/licenses/>.

#' Benchmark mice imputation
#' 
#' Hide a percent of the observations and compare
#' the results of imputation via \link[mice]{mice}
#' with the original observations.
#'
#' @param geno A data table in BIMBAM format. Must contain
#' one row per SNP, the first three columns must be site_id,
#' minor_allele and major allele, followed by one column per
#' sample. Missing values must be encoded as NA.
#' @param snp A data table with snp information in BIMBAM
#' format. Must contain columns ID, pos and chr in that order.
#' Column ID must correspond to column site_id in geno.
#' @param outdir Output directory to write the results.
#' @param p Fraction of observations to hide.
#' @param m Number of imputations performed. See \link[mice]{mice}
#' documentation.
#' @param verbose Whether to print progress on imputation. 
#' @param seed Seed for random sambling of data and for imputation
#' if needed.
#' @param block_size Number of markers to use to impute each marker.
#'
#' @return A list with elements r, p.imputed, res, and imputed_geno_file.
#' @export
#' 
#' @importFrom magrittr %>%
benchmark_imputation <- function(geno, snp, outdir, p = 0.1 ,
                                 m = 5, verbose = FALSE, seed = NA,
                                 block_size = 100){
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
                                      seed = seed,
                                      block_size = block_size))
  
  # Check imputation output
  if( any(imp$imp$site_id != geno_hidden$site_id)){
    stop("ERROR")
  }
  if( any(colnames(imp$imp) != colnames(geno_hidden))){
    stop("ERROR")
  }
  
  # Collect imputed values
  res$imputed <- (imp$imp %>%
                    dplyr::select(-site_id, -minor_allele, -major_allele) %>%
                    as.matrix)[ii]
  res$path <- outdir
  
  # Calculate correlation
  r <- cor(res$observed, res$imputed, use = "complete.obs")
  p.imputed <- 1 - (sum(is.na(res$imputed)) / nrow(res))
  
  # Plot
  p1 <- ggplot(res, aes(x = observed, y = imputed)) +
    geom_point() +
    geom_smooth(method = "lm") +
    AMOR::theme_blackbox()
  filename <- file.path(outdir, "observed_vs_imputed.svg")
  ggsave(filename, p1, width = 5, height = 5)
  
  p1 <- res %>%
    gather(key = "Type", value = "allele_frequency", observed, imputed) %>%
    ggplot(aes(x=allele_frequency)) +
    facet_grid(Type ~ .) +
    geom_histogram(bins = 20) +
    AMOR::theme_blackbox()
  filename <- file.path(outdir, "alllele_freq_histograms.svg")
  ggsave(filename, p1, width = 12, height = 5)
  
  return(list(r = r, p.imputed = p.imputed, res = res, imputed_geno_file = imp$imputed_file, t = t))
}

#' Impute genotypes with BIMBAM
#' 
#' @param geno_file Path to mean genotype file. See BIMBAM docummentation
#' for details.
#' @param pheno_file Path to phenotype file. See BIMBAM docummentation
#' for details.
#' @param pos_file Path to SNP position file. See BIMBAM docummentation
#' for details.
#' @param bimbam BIMBAM executable.
#' @param outdir Output directory.
#' @param em_runs Number of EM algorithm runs.
#' @param em_steps Steps of each EM run.
#' @param em_clusters Number of clusters in EM algorithm.
#' @param prefix Prefix of all output files
#' 
#' @references
#' http://www.haplotype.org/download/bimbam-manual.pdf
#' 
#' @export
bimbam_impute <- function(geno_file, pheno_file, pos_file,
                          bimbam = 'bimbam',
                          outdir = "imputed/",
                          em_runs = 10,
                          em_steps = 20,
                          em_clusters = 15,
                          prefix = "imputed"){
  
  cmd <- paste(bimbam,
               "-g", geno_file,
               "-p", pheno_file,
               "-pos", pos_file,
               "-e", em_runs,
               "-s", em_steps,
               "-c", em_clusters,
               "--nobf",
               "-o", prefix,
               "-wmg",
               "-gmode", 1)
  out <- run_command(cmd)
  
  # Re-organize files
  dir.create(outdir)
  
  filename <- paste0("output/", paste(c(prefix, "log.txt"), collapse = "."))
  file.copy(filename, outdir)
  file.remove(filename)
  
  filename <- paste0("output/", paste(c(prefix, "snpinfo.txt"), collapse = "."))
  file.copy(filename, outdir)
  file.remove(filename)
  
  filename <- paste0("output/", paste(c(prefix, "mean.genotype.txt"), collapse = "."))
  file.copy(filename, outdir)
  file.remove(filename)
  imputed_file <- file.path(outdir, paste(c(prefix, "mean.genotype.txt"), collapse = "."))
  
  return(imputed_file)
}

#' Calculate kinship matrix with GEMMA
#' 
#' @param geno_file Mean genotype file in BIMBAM format. See GEMMA
#' manual for details.
#' @param pheno_file Phenotype file in BIMBAM format. See GEMMA
#' manual for details.
#' @param snp_file SNP annotation file in BIMBAM format. See GEMMA
#' manual for details.
#' @param gemma GEMMA executable.
#' @param outdir directory to write output.
#' @param prefix Prefix for output filenames.
#' 
#' @references 
#' http://www.xzlab.org/software/GEMMAmanual.pdf
#' 
#' @export
gemma_kinship <- function(geno_file, pheno_file, snp_file,
                          gemma = 'gemma',
                          outdir = "kinship/",
                          prefix = "kinship"){
  
  cmd <- paste(gemma,
               "-g", geno_file,
               "-p", pheno_file,
               "-a", snp_file,
               "-gk", 1,
               "-o", prefix)
  out <- run_command(cmd)
  
  # Re-organize files
  dir.create(outdir)
  
  filename <- paste0("output/", paste(c(prefix, "log.txt"), collapse = "."))
  file.copy(filename, outdir)
  file.remove(filename)
  
  filename <- paste0("output/", paste(c(prefix, "cXX.txt"), collapse = "."))
  file.copy(filename, outdir)
  file.remove(filename)
  kinship_file <- file.path(outdir, paste(c(prefix, "cXX.txt"), collapse = "."))
  
  return(kinship_file)
}

#' Run LMM via GEMMA
#' 
#' Tested with GEMMA 0.93b modified by bugwas package.
#' 
#' @param geno_file Mean genotype file in BIMBAM format. See GEMMA
#' manual for details.
#' @param pheno_file Phenotype file in BIMBAM format. See GEMMA
#' manual for details.
#' @param snp_file SNP annotation file in BIMBAM format. See GEMMA
#' manual for details.
#' @param kinship_file Kinship file in GEMMA format. See GEMMA
#' manual for details.
#' @param cov_file Covariates file in BIMBAM format. See GEMMA
#' manual for details.
#' @param gemma GEMMA executable.
#' @param outdir directory to write output.
#' @param maf Minor allele frequency threshold for SNPs to test.
#' @param prefix Prefix for output filenames.
#' 
#' @references 
#' http://www.xzlab.org/software/GEMMAmanual.pdf
#' 
#' @export
gemma_lmm <- function(geno_file, pheno_file, snp_file, kinship_file,
                      cov_file = NULL,
                      gemma = "gemma",
                      outdir = "lmm/",
                      maf = 0,
                      prefix = "lmm"){
  
  # Run lmm
  cmd <- c(gemma,
           "-g", geno_file,
           "-p", pheno_file,
           "-a", snp_file,
           "-k", kinship_file,
           "-lmm", 2,
           "-o", prefix,
           "-maf", maf)
  if(!is.null(cov_file)){
    cmd <- c(cmd,
             "-c", cov_file)
  }

  cmd <- paste(cmd, collapse = " ")
  out <- run_command(cmd)
  
  # Re-organize files
  dir.create(outdir)
  
  filename <- paste0("output/", paste(c(prefix, "log.txt"), collapse = "."))
  file.copy(filename, outdir)
  file.remove(filename)
  lmm_log_file <- file.path(outdir, paste(c(prefix, "log.txt"), collapse = "."))
  
  filename <- paste0("output/", paste(c(prefix, "assoc.txt"), collapse = "."))
  file.copy(filename, outdir)
  file.remove(filename)
  lmm_assoc_file <- file.path(outdir, paste(c(prefix, "assoc.txt"), collapse = "."))
  
  return(c(lmm_log_file, lmm_assoc_file))
}

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
#' @param block_size Number of markers to use to impute each marker.
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
                        seed = NA,
                        block_size = 100){
  
  if(any(geno$site_id != snp$ID)){
    stop("ERROR: geno and snp tables do not match", call. = TRUE)
  }
  
  imp <- geno %>%
    split(snp$chr) %>%
    purrr::map_dfr(~tidy_mice(.), m1 = m, verbose = verbose, seed = seed, block_size = block_size) %>%
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
#' @param d A tibble
#' @param m Number of imputations
#' @param verbose Print info while running
#' @param seed Seed for mice
#' @param block_size Number of markers to use to impute each marker.
#'
#' @return A tibble with imputed results
#' 
#' @importFrom magrittr %>%
tidy_mice <- function(d, m = 5, verbose = FALSE, seed = NA, block_size = 100){
  cat("Welcome to tidy...\n")
  res <- d %>%
    dplyr::select(site_id, minor_allele, major_allele)
  
  # Remove sites that are completeley missing
  ii <- !(d %>% dplyr::select(-site_id, -minor_allele, -major_allele) %>% is.na %>% apply(1, all))
  d <- d %>% dplyr::filter(ii)
  
  if(sum(ii) == 0){
    return(res %>% dplyr::left_join(d, by = "site_id"))
  }
  
  if(block_size > 0){
    pred <- sapply(1:nrow(d), create_blocks, block_size = block_size, total_pos = nrow(d))
    diag(pred) <- 0
  }else{
    pred <- matrix(1, nrow = nrow(d), ncol = nrow(10))
    diag(pred) <- 0
  }
  
  ids <- d$site_id
  samples <- colnames(d)[-(1:3)]
  d <- d %>%
    dplyr::select(-site_id, -minor_allele, -major_allele) %>%
    t %>% mice::mice(m = m, printFlag = verbose, seed = seed, method = "pmm", predictorMatrix = pred) %>%
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
