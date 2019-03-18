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
  
  if(any(geno$site_id != snp$ID)){
    stop("ERROR: geno and snp tables do not match", call. = TRUE)
  }
  
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
    return(NULL)
  }
  
  d <- d %>%
    dplyr::select(-minor_allele, -major_allele) %>%
    mice::mice(m = m, printFlag = verbose, seed = seed) %>% 
    mice::complete() %>%
    dplyr::as_tibble()
  
  res %>% dplyr::left_join(d, by = "site_id")
}
