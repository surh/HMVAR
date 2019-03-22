# (C) Copyright 2019 Sur Herrera Paredes
# 
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

# Generic midas management functions

#' Read MIDAS abundance file
#' 
#' Reads either the snps_depth.txt or snps_freq.txt
#' file produced by midas_merge.py
#'
#' @param file File path
#'
#' @return A tibble
#' @export
#'
#' @importFrom readr read_tsv
read_midas_abun <- function(file){
  abun <- readr::read_tsv(file,
                          na = 'NA', n_max = 10)
  abun <- readr::read_tsv(file,
                          na = 'NA',
                          col_types = paste(c('c',rep('n', ncol(abun) - 1)),
                                            collapse = ''))
  
  return(abun)
}

#' Read midas snp merge data
#' 
#' Reads the output of midas_merge.py
#'
#' @param midas_dir Directory where the output from
#' midas_merge.py is located. It must include files:
#' 'snps_info.txt', 'snps_depth.txt' and 'snps_freq.txt'.
#' @param map A data table associating samples
#' in the MIDAS otput to groups. It must have an 'sample'
#' and a 'Group' column.
#' @param genes The list of genes that are to be tested.
#' Must correspond to entries in the 'genes_id' column
#' of the 'snps_info.txt' file. If NULL, all genes will
#' be tested.
#' @param cds_only If TRUE, it will remove all sites that
#' do not have a gene_id. If a gene list is passed its effect
#' is redundant.
#'
#' @return A list with three data tables: info, freq
#' and depth, corresponding to the three files from 
#' MIDAS
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @importFrom dplyr select filter
#' @importFrom tidyselect starts_with
read_midas_data <- function(midas_dir, map, genes, cds_only = TRUE){
  # Read data
  info <- readr::read_tsv(paste0(midas_dir, "/snps_info.txt"),
                          col_types = 'ccncccnnnnnccccc',
                          na = 'NA')
  depth <- read_midas_abun(paste0(midas_dir, "/snps_depth.txt"))
  freq <- read_midas_abun(paste0(midas_dir, "/snps_freq.txt"))
  
  # Process data
  # Clean info
  info <- info %>% 
    dplyr::select(-locus_type, -tidyselect::starts_with("count_"))
  # Clean depth and freq
  depth <- select_samples_from_abun(depth, map)
  freq <- select_samples_from_abun(freq, map)
  # Clean map
  map <- map %>% 
    dplyr::filter(sample %in% colnames(depth))
  
  # Select gene data
  if(!is.null(genes)){
    info <- info %>%
      dplyr::filter(gene_id %in% genes)
  }else if(cds_only){
    info <- info %>% 
      dplyr::filter(!is.na(gene_id))
  }
  freq <- freq %>% 
    dplyr::filter(site_id %in% info$site_id)
  depth <- depth %>% 
    dplyr::filter(site_id %in% info$site_id)
  
  return(list(info = info, freq = freq, depth = depth))
}


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
#' @param focal_group A character string with the group that is going
#' to be compared against everything else. If NULL, then the function
#' assumes that column group contains only two levels.
#' @param prefix Prefix to append to all filenames.
#' 
#' @return A list with elements filenames and Dat. The first element
#' contains the relative paths to the three BIMBAM files, and the
#' second contains tibbles with the data written to those files
#' 
#' @export
#' 
#' 
midas_to_bimbam <- function(midas_dir, map, outdir, focal_group = NULL,
                            prefix = NULL){
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
    arrange(factor(sample, levels = colnames(geno)[-(1:3)]))
  if(is.null(focal_group)){
    pheno <- pheno %>%
      # mutate(phenotype = 1*(Group == "Supragingival.plaque")) %>%
      mutate(phenotype = as.numeric(factor(Group)) - 1) %>%
      select(id = sample, phenotype)
  }else{
    pheno <- pheno %>%
      mutate(phenotype = 1*(Group == focal_group)) %>%
      # mutate(phenotype = as.numeric(factor(Group)) - 1) %>%
      select(id = sample, phenotype)
  }
  
  snp <- Dat$info %>% select(ID = site_id, pos = ref_pos, chr = ref_id)  %>%
    mutate(chr = as.numeric(factor(chr)))
  
  # Write bimbam tables
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  
  gen_file <- file.path(outdir, paste(c(prefix, 'geno.bimbam'), collapse = "_"))
  write_tsv(geno, path = gen_file, col_names = FALSE)
  # write_csv(geno, path = gen_file, col_names = FALSE, na = '??')
  # write.table(geno, gen_file, sep = ", ", na = 'NA', col.names = FALSE, quote = FALSE, row.names = FALSE)
  
  phen_file <- file.path(outdir, paste(c(prefix, 'pheno.bimbam'), collapse = "_"))
  # write_tsv(pheno, path = phen_file)
  write_tsv(pheno %>% select(phenotype),
            path = phen_file, col_names = FALSE)
  
  snp_file <- file.path(outdir, paste(c(prefix, 'snp.bimbam'), collapse = "_"))
  write_tsv(snp, path = snp_file, col_names = FALSE)
  # write_csv(snp, path = snp_file, col_names = FALSE)
  
  
  return(list(filenames = list(geno_file = gen_file,
                               pheno_file = phen_file,
                               snp_file = snp_file),
              Dat = list(geno = geno,
                         pheno = pheno,
                         snp = snp)))
}
