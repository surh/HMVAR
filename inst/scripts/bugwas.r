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


#' Convert MIDAS merge output to BIMBAM input
#' 
#' @param midas_dir Directory where midas merge output for one genome
#' is located. Must contain files snps_info.txt, snps_depth.txt and
#' snps_freq.txt
#' @param map Data frame or tibble that maps samples to groups. It
#' must have columns 'sample' and 'ID'.
#' @param outdir Directory where to store the results. It will be
#' created if it does not exists already. If it exists, and files
#' with the output file names exist, they will be overwriten.
#' @param prefix Prefix to append to all filenames.
#' 
#' 
#' 
#' @export
#' 
#' 
midas_to_bimbam <- function(midas_dir, map, outdir, prefix = NULL){
  Dat <- read_midas_data(midas_dir = midas_dir,
                         map = map,
                         genes = NULL,
                         cds_only = FALSE)
  
  # Keep only full covered
  # Dat$freq <- Dat$freq %>% filter(rowSums(Dat$depth[,-1] == 0) == 0)
  # Dat$info <- Dat$info %>% filter(rowSums(Dat$depth[,-1] == 0) == 0)
  # Dat$depth <- Dat$depth %>% filter(rowSums(Dat$depth[,-1] == 0) == 0)
  
  # Match freqs and depth
  Dat$depth <- Dat$depth %>% gather(key = "sample", value = 'depth', -site_id)
  Dat$freq <- Dat$freq %>% gather(key = "sample", value = 'freq', -site_id)
  Dat$info <- Dat$info %>% select(site_id, ref_id, ref_pos, major_allele, minor_allele)
  
  # Set sites without coverage to NA
  dat <- Dat$depth %>%
    inner_join(Dat$freq, by = c("site_id", "sample"))
  dat$freq[ dat$depth < 1 ] <- NA
  Dat$freq <- dat %>% select(-depth) %>% spread(sample, freq)
  
  # Create BIMBAM tables
  geno <- Dat$info %>%
    select(site_id, minor_allele, major_allele) %>%
    left_join(Dat$freq, by = "site_id")  
  
  pheno <- map %>%
    filter(sample %in% colnames(geno)) %>%
    arrange(factor(sample, levels = colnames(geno)[-(1:3)])) %>%
    mutate(phenotype = 1*(Group == "Supragingival.plaque")) %>%
    select(id = sample, phenotype)
  
  snp <- Dat$info %>% select(ID = site_id, pos = ref_pos, chr = ref_id)
  
  # Write bimbam tables
  dir.create(args$outdir)
  dir.create(paste0(args$outdir, "/bimbam/"))
  gen_file <- paste0(args$outdir, "/bimbam/geno.bimbam")
  write_tsv(geno, path = gen_file, col_names = FALSE)
  # write_csv(geno, path = gen_file, col_names = FALSE, na = '??')
  # write.table(geno, gen_file, sep = ", ", na = 'NA', col.names = FALSE, quote = FALSE, row.names = FALSE)
  
  phen_file <- paste0(args$outdir, "/bimbam/pheno.bimbam")
  # write_tsv(pheno, path = phen_file)
  write_tsv(pheno %>% select(phenotype),
            path = phen_file, col_names = FALSE)
  
  snp_file <- paste0(args$outdir, "/bimbam/snp.bimbam")
  write_tsv(snp, path = snp_file, col_names = FALSE)
  # write_csv(snp, path = snp_file, col_names = FALSE)
  
  
  
}





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

