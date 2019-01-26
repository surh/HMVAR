library(HMVAR)
library(tidyverse)
library(bugwas)

Sys.setenv(LD_LIBRARY_PATH="/opt/modules/pkgs/eqtlbma/git/lib/")


args <- list(midas_dir = "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/Actinomyces_odontolyticus_57475/",
             map_file = "hmp_SPvsTD_map.txt",
             outdir = "testout/",
             prefix = "Actinomyces_odontolyticus_57475",
             gemma = "~/bin/gemma.0.93b")



# test_genes <- c("411466.7.peg.516", "411466.7.peg.602",
#                 "411466.7.peg.603", "411466.7.peg.604",
#                 "411466.7.peg.605", "411466.7.peg.606",
#                 "411466.7.peg.607", "411466.7.peg.608",
#                 "411466.7.peg.609", "411466.7.peg.965",
#                 "411466.7.peg.975", "411466.7.peg.1738")




map <- read_tsv(args$map_file, col_types = 'cc')
map <- map %>% select(sample = ID, Group = Group)


res <- midas_to_bimbam(midas_dir = args$midas_dir, map = map, outdir = 'bimbam', prefix = 'Actino')



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

