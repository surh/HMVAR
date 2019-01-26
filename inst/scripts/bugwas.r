library(HMVAR)
library(tidyverse)
library(bugwas)

# The following will be based on bugwas:
# Steps
# 1. Impute genotypes with BIMBAM

Sys.setenv(LD_LIBRARY_PATH="/opt/modules/pkgs/eqtlbma/git/lib/")
args <- list(midas_dir = "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/Actinomyces_odontolyticus_57475/",
             map_file = "hmp_SPvsTD_map.txt",
             outdir = "testout/",
             prefix = "Actinomyces_odontolyticus_57475",
             gemma = "~/bin/gemma.0.93b",
             bimbam = "~/bin/bimbam")

# test_genes <- c("411466.7.peg.516", "411466.7.peg.602",
#                 "411466.7.peg.603", "411466.7.peg.604",
#                 "411466.7.peg.605", "411466.7.peg.606",
#                 "411466.7.peg.607", "411466.7.peg.608",
#                 "411466.7.peg.609", "411466.7.peg.965",
#                 "411466.7.peg.975", "411466.7.peg.1738")


# Main output directory
dir.create(args$outdir)

# Read map
map <- read_tsv(args$map_file, col_types = 'cc')
map <- map %>% select(sample = ID, Group = Group)

# Convert to bimbam
bimbam_dir <- file.path(args$outdir, "bimbam")
midas_bimbam <- midas_to_bimbam(midas_dir = args$midas_dir, map = map, outdir = bimbam_dir, prefix = NULL)

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
file.copy("output/imputed.log.txt", imputed_dir)
file.copy("output/imputed.mean.genotype.txt", imputed_dir)
file.copy("output/imputed.snpinfo.txt", imputed_dir)
file.remove("output/imputed.log.txt")
file.remove("output/imputed.mean.genotype.txt")
file.remove("output/imputed.snpinfo.txt")






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

