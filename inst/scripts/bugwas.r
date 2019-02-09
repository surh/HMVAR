library(HMVAR)
library(tidyverse)
library(bugwas)

#' Run system command
#' 
#' Internal
#' 
#' @param cmd System command
run_command <- function(cmd){
  cat("Running\n\t>", cmd, "\n")
  out <- system(cmd)
  
  return(out)
}

# Run bimbam
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

# The following will be based on bugwas:
# Steps
# 1. Impute genotypes with BIMBAM
# 2. Get kinship matrix
# 3. Run lmm

Sys.setenv(LD_LIBRARY_PATH="/opt/modules/pkgs/eqtlbma/git/lib/")
# indir <- commandArgs(trailingOnly = TRUE)[1]
# spec <- commandArgs(trailingOnly = TRUE)[2]
# indir <- "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/"
# spec <- "Actinomyces_odontolyticus_57475"
indir <- "/godot/users/sur/exp/fraserv/2019/today/"
spec <- "midas_output_small/"

args <- list(midas_dir = file.path(indir, spec),
             map_file = "map.txt",
             outdir = "metawas",
             prefix = spec,
             gemma = "~/bin/gemma.0.93b",
             bimbam = "~/bin/bimbam",
             gemma_version = 'bugwas')
rm(indir, spec)

# test_genes <- c("411466.7.peg.516", "411466.7.peg.602",
#                 "411466.7.peg.603", "411466.7.peg.604",
#                 "411466.7.peg.605", "411466.7.peg.606",
#                 "411466.7.peg.607", "411466.7.peg.608",
#                 "411466.7.peg.609", "411466.7.peg.965",
#                 "411466.7.peg.975", "411466.7.peg.1738")


# Main output directory
dir.create(args$outdir)

# Create list for filenames
Files <- list(Dirs = list(),
              Files = list())

# Read map
map <- read_tsv(args$map_file, col_types = 'cc')
map <- map %>% select(sample = ID, Group = Group)

# Convert to bimbam
Files$Dirs$bimbam_dir <- file.path(args$outdir, "bimbam")
midas_bimbam <- midas_to_bimbam(midas_dir = args$midas_dir,
                                map = map,
                                outdir = Files$Dirs$bimbam_dir,
                                prefix = NULL)
Files$Files$midas_geno_file <- midas_bimbam$filenames$geno_file
Files$Files$pheno_file <- midas_bimbam$filenames$pheno_file
Files$Files$snp_file <- midas_bimbam$filenames$snp_file
rm(map)

# Impute
Files$Dirs$imputed_dir <- file.path(args$outdir, "imputed")
Files$Files$imputed_geno_file <- bimbam_impute(geno_file = midas_bimbam$filenames$geno_file,
                                               pheno_file = midas_bimbam$filenames$pheno_file,
                                               pos_file = midas_bimbam$filenames$snp_file,
                                               bimbam = args$bimbam,
                                               outdir = Files$Dirs$imputed_dir,
                                               em_runs = 10,
                                               em_steps = 20,
                                               em_clusters = 15,
                                               prefix = "imputed")

# Get kinship matrix
# Works with both gemma v0.93b & v0.98.1
# I am ingoring patterns since genotypes are not fixed but frequencies instead
Files$Dirs$kinship_dir <- file.path(args$outdir, "kinship")
Files$Files$kinship_file <- gemma_kinship(geno_file = Files$Files$imputed_geno_file,
                                          pheno_file = Files$Files$pheno_file,
                                          snp_file = Files$Files$snp_file,
                                          gemma = args$gemma,
                                          outdir = Files$Dirs$kinship_dir,
                                          prefix = 'kinship')

# Run lmm
Files$Dirs$lmm_dir <- file.path(args$outdir, "lmm")
res <- gemma_lmm(geno_file = Files$Files$imputed_geno_file,
          pheno_file = Files$Files$pheno_file,
          snp_file = Files$Files$snp_file,
          kinship_file = Files$Files$kinship_file,
          gemma = args$gemma,
          outdir = Files$Dirs$lmm_dir,
          maf = 0,
          prefix = "lmm")
Files$Files$lmm_log_file <- res[1]
Files$Files$lmm_assoc_file <- res[2]
rm(res)


# Prepare data for gemma
Dat_gemma <- list(geno = read_table2(file.path(Files$Dirs$imputed_dir, "imputed.mean.genotype.txt"),
                                     col_names = colnames(midas_bimbam$Dat$geno),
                                     col_types = paste0('c', 'c', 'c',
                                                        paste(rep('n', ncol(midas_bimbam$Dat$geno) - 3),
                                                              collapse ="" ))) %>%
                    arrange(factor(site_id, levels = midas_bimbam$Dat$snp$ID)),
                  pheno = midas_bimbam$Dat$pheno,
                  snp = midas_bimbam$Dat$snp)
# Cleanup
rm(midas_bimbam, cmd, out)
gc()

# Run lmm
cmd <- paste(args$gemma,
             "-g", Files$Files$imputed_geno_file,
             "-p", Files$Files$pheno_file,
             "-a", Files$Files$snp_file,
             "-k", Files$Files$kinship_file,
             "-lmm", 2,
             "-o", "lmm",
             "-maf", 0)
cat("Running\n\t>", cmd, "\n")
out <- system(cmd)
# Re-organize files
Files$Dirs$lmm_dir <- file.path(args$outdir, "lmm")
dir.create(Files$Dirs$lmm_dir)
file.copy("output/lmm.assoc.txt", Files$Dirs$lmm_dir)
file.copy("output/lmm.log.txt", Files$Dirs$lmm_dir)
file.remove("output/lmm.log.txt")
file.remove("output/lmm.assoc.txt")
Files$Files$lmm_assoc_file <- file.path(Files$Dirs$lmm_dir, "lmm.assoc.txt")
Files$Files$lmm_log_file <- file.path(Files$Dirs$lmm_dir, "lmm.log.txt")

# Process lmm results
# Get lognull and lambda
# For newer GEMMA versions I need vg/ve from two diff lines.
lognull <- scan(Files$Files$lmm_log_file,
                what = character(0),
                sep = "\n")[17] %>%
  strsplit(" ") %>%
  unlist %>%
  last %>%
  as.numeric
lambda <- scan(Files$Files$lmm_log_file,
                what = character(0),
                sep = "\n")[13] %>%
  strsplit(" ") %>%
  unlist %>%
  last %>%
  as.numeric

# # # Read results
# # lmm_res <- read_tsv(Files$Files$lmm_assoc_file,
# #                     col_types = cols(rs = 'c')) %>% 
# #   select(chr, rs, ps, lambda_H1 = l_mle,
# #          loglik_H1 = logl_H1, p.value = p_lrt) %>%
# #   mutate(minus.log.p = -log10(p.value))
# # 
# # lmm <- list(lmm = lmm_res,
# #             lognull = as.numeric(lognull),
# #             lambda = as.numeric(lambda))
# lmm_res <- read_tsv(Files$Files$lmm_assoc_file,
#                     col_types = cols(rs = 'c')) %>%
#   mutate(negLog10 = -log10(p_lrt))
# 
# lmm.bi <- list(lmm = lmm_res, "lognull" = as.numeric(lognull),
#      lambda = as.numeric(lambda))
#      
# cor.geno <- bugwas:::get_correlations(XX = geno,
#                                       pca = geno.pca$x,
#                                       npcs = length(Dat_gemma$pheno$id),
#                                       id = Dat_gemma$pheno$id)
# lmm <- list(logreg.bi = NULL, lmm.bi = lmm.bi,
#             lognull = lognull,
#             lambda = lambda, cor.XX = cor.geno)


# # SVD & PCA
# # Since there are no patterns all genotypes have the same weighth
geno <- Dat_gemma$geno %>% select(-site_id, -minor_allele, -major_allele) %>%
  as.matrix %>% t
geno.svd <- svd(geno)
geno.pca <- prcomp(geno)
# rm(geno)

# Wald test
# wald <- wald_test(y = y,
#                   XX = XX,
#                   svd.XX = svd.XX,
#                   lambda = biallelic$lambda,
#                   XX.all = XX.all,
#                   prefix = prefix,
#                   npcs = npcs,
#                   pca = pca$pca)

# fit.lmm <- ridge_regression(y, XX, svdX=svd.XX,
#                             lambda_init=as.numeric(lambda)/sum(XX.all$bippat),
#                             maximize=FALSE, skip.var=TRUE)
m1.ridge <- bugwas:::ridge_regression(y = Dat_gemma$pheno$phenotype, x = geno,
                                      svdX = geno.svd, lambda_init = lambda / nrow(Dat_gemma$geno),
                                      maximize = FALSE, skip.var = TRUE)

# Fit the grand null model
# fit.0 <- lm(y~1)
m0 <- lm(phenotype ~ 1, data = Dat_gemma$pheno)

# LRT for the LMM null vs grand null
# LRTnullVgrand <- -log10(pchisq(2*(fit.lmm$ML - as.numeric(logLik(fit.0))), 1, low=F)/2)
# cat(paste0("## LRT for the LMM null vs grand null = ", LRTnullVgrand),
#     file = paste0(prefix, "_logfile.txt"), sep="\n", append = TRUE)
lrt.pval <- pchisq(2*(m1.ridge$ML - as.numeric(logLik(m0))), df = 1, lower.tail = FALSE) / 2

# Heritability
# fit.lmm.ypred <- XX %*% fit.lmm$Ebeta
# cat(paste0("## Heritability (R^2) = ", cor(fit.lmm.ypred,y)^2),
#     file=paste0(prefix, "_logfile.txt"), sep="\n", append=TRUE)
y.pred <- geno %*% m1.ridge$Ebeta
# plot(y.pred, Dat_gemma$pheno$phenotype)
h2 <- cor(y.pred, Dat_gemma$pheno$phenotype)^2
# h2

# Get full posterior covariance matrix for Bayesian Wald Test
# Need the full posterior covariance matrix for the Bayesian Wald test,
# to get the posterior uncertainty for each point estimate
# wald_input <- get_wald_input(fit.lmm = fit.lmm, pca = pca, svd.XX = svd.XX,
#                              y = y, npcs = npcs, XX = XX)
wald_input <- bugwas:::get_wald_input(fit.lmm = m1.ridge,
                                      pca = geno.pca,
                                      svd.XX = geno.svd,
                                      y = Dat_gemma$pheno$phenotype,
                                      npcs = length(Dat_gemma$pheno$phenotype),
                                      XX = geno)

# Bayesian Wald Test
# pca.bwt <- wald_input$Ebeta^2/diag(wald_input$Vbeta)
# p.pca.bwt <- -log10(exp(1))*pchisq(pca.bwt, 1, low=F, log=T)
# cat(paste0("## Bayesian Wald Test for PCs range = ", paste(range(p.pca.bwt), collapse=" ")),
#     file=paste0(prefix, "_logfile.txt"), sep="\n", append=TRUE)
# write.table(p.pca.bwt, file = paste0(prefix, "_Bayesian_Wald_Test_negLog10.txt"),
#             sep="\t", row=T, col = F, quote=F)
bwt.pvals <- -log10(exp(1)) *
  pchisq(wald_input$Ebeta^2 / diag(wald_input$Vbeta),
         df = 1, lower.tail = FALSE, log.p = TRUE)
# bwt.pvals

# Predict phenotype using effect sizes.
# There are small numerical differences from this prediction but it is the same
# as y.pred
# effect <- t(t(XX) * as.vector(fit.lmm$Ebeta))
# pred <- rowSums(effect)
# pheno.pred <- rowSums(t(t(geno) * as.vector(m1.ridge$Ebeta)))


# y.pred <- geno %*% m1.ridge$Ebeta
# pheno.pred <- t(geno) * as.vector(m1.ridge$Ebeta)
# summary(colSums(pheno.pred) - as.vector(y.pred))
# 
# mat <- matrix(1:10, ncol = 5 )
# beta <- matrix(1:5, ncol = 1)
# mat
# beta
# mat %*% beta
# 
# t(mat)
# t(mat) * as.vector(beta)
# colSums(t(mat) * as.vector(beta))

# return(list("pc_order" = pc_order, "p.pca.bwt" = p.pca.bwt, "pred" = pred,
#             "signif_cutoff" = signif_cutoff))









# # Get a dendrogram of samples
# # TEMPORARY. MOVE EARLIER
# tre <- hclust(dist(dist(geno)))
# tre <- ape::as.phylo(tre)
# Files$Files$phylo_file <- file.path(args$outdir, "bimbam/phylo.tre")
# ape::write.tree(phy = tre, file = Files$Files$phylo_file)
# 
# # Get list of all tree info
# treeInfo <- bugwas:::get_tree(phylo = Files$Files$phylo_file,
#                               prefix = 'actino',
#                               XX.ID = Dat_gemma$pheno$id,
#                               pca = geno.pca,
#                               npcs = length(Dat_gemma$pheno$id),
#                               allBranchAndPCCor = args$PC_branch_cors)
# 
# 
# 
# # Ridge regression
# wald <- bugwas:::wald_test(y = Dat_gemma$pheno$phenotype,
#                            XX = geno,
#                            svd.XX = geno.svd,
#                            lambda = lmm$lambda,
#                            XX.all = XX.all,
#                            prefix = 'actino',
#                            npcs = length(Dat_gemma$pheno$id),
#                            pca = geno.pca)
# 
# 
# 
# ggplot(lmm_res, aes(x = ps, y = negLog10)) +
#   facet_grid(~chr, space = "free_x", scales = "free_x") +
#   geom_hline(yintercept = 8, color = "red", size = 3) +
#   geom_point()
# 
# ggplot(lmm_res, aes(x = p_lrt)) +
#   geom_histogram(bins = 20) +
#   AMOR::theme_blackbox()
