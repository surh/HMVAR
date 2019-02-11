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
#' @References
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
                      gemma = "gemma",
                      outdir = "lmm/",
                      maf = 0,
                      prefix = "lmm"){
  
  # Run lmm
  cmd <- paste(gemma,
               "-g", geno_file,
               "-p", pheno_file,
               "-a", snp_file,
               "-k", kinship_file,
               "-lmm", 2,
               "-o", prefix,
               "-maf", maf)
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
