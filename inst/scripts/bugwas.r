library(HMVAR)
library(tidyverse)
library(bugwas)

# The following will be based on bugwas:
# Steps
# 1. Impute genotypes with BIMBAM
# 2. Get kinship matrix
# 3. Run lmm

Sys.setenv(LD_LIBRARY_PATH="/opt/modules/pkgs/eqtlbma/git/lib/")
args <- list(midas_dir = "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/Actinomyces_odontolyticus_57475/",
             map_file = "hmp_SPvsTD_map.txt",
             outdir = "testout/",
             prefix = "Actinomyces_odontolyticus_57475",
             gemma = "~/bin/gemma.0.93b",
             bimbam = "~/bin/bimbam",
             loglik_null_line = 17)

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
bimbam_dir <- file.path(args$outdir, "bimbam")
Files$Dirs$bimbam_dir <- bimbam_dir
midas_bimbam <- midas_to_bimbam(midas_dir = args$midas_dir, map = map, outdir = bimbam_dir, prefix = NULL)
Files$Files$midas_geno_file <- midas_bimbam$filenames$geno_file
Files$Files$pheno_file <- midas_bimbam$filenames$pheno_file
Files$Files$snp_file <- midas_bimbam$filenames$snp_file

# Run bimbam
cmd <- paste(args$bimbam,
             "-g", midas_bimbam$filenames$geno_file,
             "-p", midas_bimbam$filenames$pheno_file,
             "-pos", midas_bimbam$filenames$snp_file,
             "-e", 10,
             "-s", 20,
             "-c", 15,
             "--nobf",
             "-o", "imputed",
             "-wmg",
             "-gmode", 1)
cat("Running\n\t>", cmd, "\n")
out <- system(cmd)
# Re-organize files
imputed_dir <- file.path(args$outdir, "imputed")
dir.create(imputed_dir)
Files$Dirs$imputed_dir <- imputed_dir
file.copy("output/imputed.log.txt", imputed_dir)
file.copy("output/imputed.mean.genotype.txt", imputed_dir)
file.copy("output/imputed.snpinfo.txt", imputed_dir)
file.remove("output/imputed.log.txt")
file.remove("output/imputed.mean.genotype.txt")
file.remove("output/imputed.snpinfo.txt")
Files$Files$imputed_geno_file <- file.path(imputed_dir, "imputed.mean.genotype.txt")

# Prepare data for gemma
Dat_gemma <- list(geno = read_table2(file.path(imputed_dir, "imputed.mean.genotype.txt"),
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

# Get kinship matrix
# Works with both gemma v0.93b & v0.98.1
# I am ingoring patterns since genotypes are not fixed but frequencies instead
cmd <- paste(args$gemma,
             "-g", Files$Files$imputed_geno_file,
             "-p", Files$Files$pheno_file,
             "-a", Files$Files$snp_file,
             "-gk", 1,
             "-o", "kinship")
cat("Running\n\t>", cmd, "\n")
out <- system(cmd)
# Re-organize files
Files$Dirs$kinship_dir <- file.path(args$outdir, "kinship")
dir.create(Files$Dirs$kinship_dir)
file.copy("output/kinship.cXX.txt", Files$Dirs$kinship_dir)
file.copy("output/kinship.log.txt", Files$Dirs$kinship_dir)
file.remove("output/kinship.log.txt")
file.remove("output/kinship.cXX.txt")
Files$Files$kinship_file <- file.path(Files$Dirs$kinship_dir, "kinship.cXX.txt")


# SVD & PCA
# Since there are no patterns all genotypes have the same weighth
geno <- Dat_gemma$geno %>% select(-site_id, -minor_allele, -major_allele) %>%
  as.matrix %>% t
geno.svd <- svd(geno)
geno.pca <- prcomp(geno)

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
                sep = "\n")[args$loglik_null_line] %>%
  strsplit(" ") %>%
  unlist %>%
  last
lambda <- scan(Files$Files$lmm_log_file,
                what = character(0),
                sep = "\n")[13] %>%
  strsplit(" ") %>%
  unlist %>%
  last

# Read results
lmm_res <- read_tsv(Files$Files$lmm_assoc_file,
                    col_types = cols(rs = 'c')) %>% 
  select(chr, rs, ps, lambda_H1 = l_mle,
         loglik_H1 = logl_H1, p.value = p_lrt) %>%
  mutate(minus.log.p = -log10(p.value))

lmm <- list(lmm = lmm_res,
            lognull = as.numeric(lognull),
            lambda = as.numeric(lambda))





# Get a dendrogram of samples

# tre <- hclust(dist(t(Dat$freq[,-1])))
tre <- hclust(dist(t(apply(t(Dat$depth[,-1]),1, scale))))
# tre <- hclust(dist(t(Dat$depth[,-1])))
# plot(tre)
tre <- ape::as.phylo(tre)


phylo_file <- paste0(args$outdir, "/bimbam/phylo.tre")
ape::write.tree(phy = tre, file = phylo_file)


rm(Dat, dat, map, geno, pheno, snp, tre)
gc()

## Call bugwas
b1 <- linLocGEMMA(gemmaGenFile = gen_file,
                  gemmaSnpFile = snp_file,
                  pheno = phen_file,
                  phylo = phylo_file,
                  prefix = args$prefix,
                  gem.path = args$gemma)

