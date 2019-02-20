library(HMVAR)
library(tidyverse)
library(bugwas)

#' Function to test PCs as composite of genetic effects
#' 
#' From bugwas
#' 
#' @export
bwt_pcs <- function(pheno, geno, x.svd, x.pca, lambda.init){
  # NOTE: move svd and PCA inside?
  
  # fit.lmm <- ridge_regression(y, XX, svdX=svd.XX,
  #                             lambda_init=as.numeric(lambda)/sum(XX.all$bippat),
  #                             maximize=FALSE, skip.var=TRUE)
  m1.ridge <- bugwas:::ridge_regression(y = pheno,
                                        x = geno,
                                        svdX = x.svd,
                                        lambda_init = lambda.init,
                                        maximize = FALSE,
                                        skip.var = TRUE)
  
  # Fit the grand null model, why?
  # fit.0 <- lm(y~1)
  m0 <- lm(pheno ~ 1)
  
  # LRT for the LMM null vs grand null
  # LRTnullVgrand <- -log10(pchisq(2*(fit.lmm$ML - as.numeric(logLik(fit.0))), 1, low=F)/2)
  # cat(paste0("## LRT for the LMM null vs grand null = ", LRTnullVgrand),
  #     file = paste0(prefix, "_logfile.txt"), sep="\n", append = TRUE)
  lrt.pval <- pchisq(2*(m1.ridge$ML - as.numeric(logLik(m0))),
                     df = 1, lower.tail = FALSE) / 2
  
  # Heritability
  # fit.lmm.ypred <- XX %*% fit.lmm$Ebeta
  # cat(paste0("## Heritability (R^2) = ", cor(fit.lmm.ypred,y)^2),
  #     file=paste0(prefix, "_logfile.txt"), sep="\n", append=TRUE)
  y.pred <- geno %*% m1.ridge$Ebeta
  # plot(y.pred, Dat_gemma$pheno$phenotype)
  h2 <- cor(y.pred, pheno)^2
  # h2
  
  # Get full posterior covariance matrix for Bayesian Wald Test
  # Need the full posterior covariance matrix for the Bayesian Wald test,
  # to get the posterior uncertainty for each point estimate
  # wald_input <- get_wald_input(fit.lmm = fit.lmm, pca = pca, svd.XX = svd.XX,
  #                              y = y, npcs = npcs, XX = XX)
  wald_input <- bugwas:::get_wald_input(fit.lmm = m1.ridge,
                                        pca = x.pca,
                                        svd.XX = x.svd,
                                        y = pheno,
                                        npcs = length(pheno),
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
  
  # NOTE: write output
  return(bwt.pvals)
}


# The following will be based on bugwas:
# Steps
# 1. Impute genotypes with BIMBAM
# 2. Get kinship matrix
# 3. Run lmm
# 4. Perform bayesian Wald test on principal components

# Required in fraserv for the bugwas modified gemma
Sys.setenv(LD_LIBRARY_PATH="/opt/modules/pkgs/eqtlbma/git/lib/")

# Setting up options for test
# indir <- commandArgs(trailingOnly = TRUE)[1]
# spec <- commandArgs(trailingOnly = TRUE)[2]
# indir <- "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/"
# spec <- "Actinomyces_odontolyticus_57475"
indir <- "/godot/users/sur/exp/fraserv/2019/2019-02-08.test_metawas/"
spec <- "midas_output_small/"

# Eventually replace this with argparse
args <- list(midas_dir = file.path(indir, spec),
             map_file = "map.txt",
             outdir = "metawas",
             prefix = spec,
             gemma = "~/bin/gemma.0.93b",
             bimbam = "~/bin/bimbam",
             gemma_version = 'bugwas',
             pcs = "hmp_SPvsTD_relabun_10pcs.txt")
rm(indir, spec)

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


# Get covariates
pcs <- read_tsv(args$pcs)
pcs <- pcs %>% slice(match(midas_bimbam$Dat$pheno$id, ID))
pcs$ID <- 1
Files$Files$pc_covariates <- file.path(Files$Dirs$bimbam_dir, "pcs.bimbam")
write_tsv(pcs, Files$Files$pc_covariates, col_names = FALSE)


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


# # Process lmm results
# # Get lognull and lambda
# # For newer GEMMA versions I need vg/ve from two diff lines.
# # lognull <- scan(Files$Files$lmm_log_file,
# #                 what = character(0),
# #                 sep = "\n")[17] %>%
# #   strsplit(" ") %>%
# #   unlist %>%
# #   last %>%
# #   as.numeric
# lambda <- scan(Files$Files$lmm_log_file,
#                what = character(0),
#                sep = "\n")[13] %>%
#   strsplit(" ") %>%
#   unlist %>%
#   last %>%
#   as.numeric
# 
# # Prepare data for gemma
# Dat_gemma <- list(geno = read_table2(file.path(Files$Dirs$imputed_dir, "imputed.mean.genotype.txt"),
#                                      col_names = colnames(midas_bimbam$Dat$geno),
#                                      col_types = paste0('c', 'c', 'c',
#                                                         paste(rep('n', ncol(midas_bimbam$Dat$geno) - 3),
#                                                               collapse ="" ))) %>%
#                     arrange(factor(site_id, levels = midas_bimbam$Dat$snp$ID)),
#                   pheno = midas_bimbam$Dat$pheno,
#                   snp = midas_bimbam$Dat$snp)
# # Cleanup
# rm(midas_bimbam, cmd, out)
# gc()
# 
# # SVD & PCA
# # Since there are no patterns all genotypes have the same weight
# # Need to recenter genotype
# geno <- Dat_gemma$geno %>% select(-site_id, -minor_allele, -major_allele) %>%
#   as.matrix %>% t
# geno <- t(t(geno) - colMeans(geno))
# geno.svd <- svd(geno)
# geno.pca <- prcomp(geno)
# 
# # Bayesian wald test on PCS
# bwt.pvals <- bwt_pcs(pheno = Dat_gemma$pheno$phenotype,
#                      geno = geno,
#                      x.svd = geno.svd,
#                      x.pca = geno.pca,
#                      lambda.init = lambda / ncol(geno))
# 
# 
# 
# 
# 
# # Get a dendrogram of samples
# # TEMPORARY. MOVE EARLIER
# tre <- hclust(dist(dist(geno)))
# tre <- ape::as.phylo(tre)
# Files$Files$phylo_file <- file.path(args$outdir, "bimbam/phylo.tre")
# ape::write.tree(phy = tre, file = Files$Files$phylo_file)
# 
# # Get list of all tree info
# tree_info <- bugwas:::get_tree(phylo = Files$Files$phylo_file,
#                                prefix = 'actino',
#                                XX.ID = Dat_gemma$pheno$id,
#                                pca = geno.pca,
#                                npcs = length(Dat_gemma$pheno$id),
#                                allBranchAndPCCor = FALSE)
# 
# 
# 
# 
# 
# 
# ### Plots
# pc_alpha <- 0.5 ## FOR TEST PURPOSES!!
# 
# # Get a random sample of colours to use for all plots, equal to the number of significant PCs
# # colourPalette <- getColourPalette(p.pca.bwt = p.pca.bwt, signifCutOff = signifCutOff, pc.lim = pc.lim)
# pc_colors <- rep("grey50", nrow(bwt.pvals))
# ii <- bwt.pvals[,1] > -log10(pc_alpha / nrow(bwt.pvals))
# pc_colors[ ii ] <- rep(c(RColorBrewer::brewer.pal(n=12, name = 'Paired'),
#                          "#000000"),
#                        length.out = sum(ii))
# 
# 
# # sampleCount = length(Dat_gemphenotype)
# # m = match(o[pc.lim], which.mtp.pc)
# n_samples <- length(Dat_gemma$pheno$phenotype)
# matched_lineages <- match(which(ii), tree_info$cor.tree$which.pc)
# 
# # pc.lim 1:n_significant_pcs or NULL if no significant
# 
# #Bayesian Wald test for genome-wide PCs
# # p.genomewidepc = .testGenomeWidePCs(prefix = prefix,
# #                                     pc.lim = pc.lim,
# #                                     pca = pca,
# #                                     bippat = bippat,
# #                                     ipat = ipat,
# #                                     o = o)
# # message("Bayesian Wald test for genome-wide PCs has been completed successfully.")
# # No patterns
# bugwas:::.testGenomeWidePCs(prefix = 'testpc',
#                             pc.lim = 1,
#                             pca = geno.pca,
#                             bippat = rep(1, ncol(geno)),
#                             ipat = 1:(ncol(geno)),
#                             o = order(bwt.pvals[,1], decreasing = TRUE))
# 
# 
# 
# 
# #The barplot for the Bayesian wald test for genome-wide PCs
# .BayesianWaldTestPCsBarplot(prefix = prefix,
#                             p.pca.bwt = p.pca.bwt,
#                             colourPalette = colourPalette,
#                             o = o,
#                             m = m,
#                             p.genomewidepc = p.genomewidepc,
#                             pc.lim = pc.lim)
# message("The barplot for the Bayesian wald test for genome-wide PCs has been completed successfully.")



















# ggplot(lmm_res, aes(x = ps, y = negLog10)) +
#   facet_grid(~chr, space = "free_x", scales = "free_x") +
#   geom_hline(yintercept = 8, color = "red", size = 3) +
#   geom_point()
# 
# ggplot(lmm_res, aes(x = p_lrt)) +
#   geom_histogram(bins = 20) +
#   AMOR::theme_blackbox()
