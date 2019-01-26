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





gemmaGen.df = read.table(file=gemmaGenFile, header=F, as.is=T)
gemmaSNP.df = read.table(file=gemmaSnpFile, header=F, as.is=T)
# bipCount = nrow(gemmaGen.df)
# dataCount = ncol(gemmaGen.df) - 3

n_snvs <- nrow(Dat_gemma$geno)
n_samples <- ncol(Dat_gemma$geno) - 3

if(nrow(Dat_gemma$pheno) != n_samples){
  stop("ERROR: Number of samples does not match between genotype and phenotype")
}
if(nrow(Dat_gemma$snp != n_snvs)){
  stop("ERROR: Number of SNVs does not match between genotype and position file")
}


# m = matrix(nrow=2,ncol=bipCount)
# gen = gemmaGen.df[,-c(1:3)]
# m[2,] = rowSums(gemmaGen.df[,-c(1:3)])
# m[1,] = dataCount - m[2,]
m <- matrix(nrow = 2, ncol = n_snvs)


allele.id <- matrix(c(0, 1)[apply(m, 2, order, decreasing=TRUE)], nrow=2)

# Output filenames
# biallelic polymorphisms encoded -1 (missing) 0 (allele 0) 1 (allele 1)
bip_outfile <- paste0(prefix, ".gemma.bip.patterns.txt");		

# positional and allelic information for biallelic polymorphisms
bipinfo_outfile <- paste0(prefix, ".gemma.bipinfo.txt");		

# Allocate memory for bip and snp, because need to transform so cannot output on the fly
bip <- matrix(NA, bipCount, dataCount)

for(i in 1:ncol(gen)) {
  # Read the mapcall file
  fa <- gen[,i]
  bip[, i] <- -1
  bip[fa==allele.id[1, ], i] = 0
  bip[fa==allele.id[2, ], i] = 1
  
}

# Convert BIP and SNP patterns to factors to identify equivalencies
bip.pat <- factor(apply(bip, 1, paste, collapse=""))
# Record only unique patterns, and record the pattern equivalence in the bipinfo file
bip.pat1 <- match(levels(bip.pat), bip.pat)

# Output compacted bip and snp objects
if(is.null(id)){
  id = paste("id", c(1:ncol(bip)), sep="")
}
colnames(bip) = id
write.table(bip[bip.pat1, ], bip_outfile, row=FALSE, col=TRUE, quote=FALSE, sep="\t")

# Output info files
pos = as.numeric(gemmaSNP.df[,2])
bipinfo <- data.frame("Position"= pos,
                      "Allele0"=allele.id[1,],
                      "Allele1"=allele.id[2,],
                      "0"=m[1,],
                      "1"=m[2,],
                      "Pattern"=as.numeric(bip.pat));

write.table(bipinfo, bipinfo_outfile, row=FALSE, col=TRUE, quote=FALSE, sep="\t")
ps.bips <- bipinfo$Position
bippat <- sapply(1:max(bipinfo$Pattern), function(x)sum(bipinfo$Pattern==x))
ipat <- bipinfo$Pattern
rm(bipinfo)		

return(list("XX" = bip[bip.pat1, ],
            "XX.tritetra" = NULL,
            "bippat" = bippat,
            "snppat" = NULL,
            "pattern" = ipat,
            "pattern.snps" = NULL,
            "ps" = ps.bips,
            "ps.snps" = NULL,
            "n.triallelic" = 0,
            "n.tetraallelic" = 0))








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

