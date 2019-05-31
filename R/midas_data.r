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

#' Checks loaded midas data
#' 
#' Internal
#' 
#' @param Dat A MIDAS data
#'
#' @return TRUE if all tests pass
check_midas_data <- function(Dat){
  
  # Check name of elements
  if(!all(c('info', 'freq', 'depth') %in% names(Dat))){
    stop("ERROR: The midas_data object must have elements c('info', 'freq', 'depth')", call. = TRUE)
  }
  
  # Check element types
  if(!tibble::is_tibble(Dat$info))
    stop("ERROR: info must be a tibble.", call. = TRUE)
  if(!tibble::is_tibble(Dat$freq))
    stop("ERROR: freq must be a tibble.", call. = TRUE)
  if(!tibble::is_tibble(Dat$depth))
    stop("ERROR: depth must be a tibble.", call. = TRUE)
  
  # Check dimensions of freq and depth
  if(any(dim(Dat$freq) != dim(Dat$depth)))
    stop("ERROR: freq and depth must have matching dimmensions.", call. = TRUE)
  
  # Check names of samples and site ids
  if(any(colnames(Dat$freq) != colnames(Dat$depth)))
    stop("ERROR: column names of depth and freq don't match")
  if(any(Dat$freq$site_id != Dat$depth$site_id))
    stop("ERROR: site_id differs between freq and depth.", call. = TRUE)
  if(any(Dat$info$site_id != Dat$freq$site_id))
    stop("ERROR: site_id differens between info and freq")
  
  return(TRUE)
}

#' Read MIDAS abundance file
#' 
#' Reads either the snps_depth.txt or snps_freq.txt
#' file produced by midas_merge.py
#'
#' @param file File path
#'
#' @return A tibble
#' @export
read_midas_abun <- function(file){
  abun <- readr::read_tsv(file,
                          na = 'NA', 
                          col_types = readr::cols(site_id = readr::col_character(),
                                                  .default = readr::col_double()))
  
  # abun <- readr::read_tsv(file,
  #                         na = 'NA', n_max = 10)
  # abun <- readr::read_tsv(file,
  #                         na = 'NA',
  #                         col_types = paste(c('c',rep('n', ncol(abun) - 1)),
  #                                           collapse = ''))
  
  return(abun)
}

#' Get all MIDAS data from a gene
#' 
#' Returns a list similar to the one by \link{read_midas_data} but
#' only with sites from a given gene.
#'
#' @param Dat A list. See output from \link{read_midas_data}
#' @param gene Name of the gene to be extracted, must exist in
#' 'gene_id' column of the 'info' element from Dat.
#' @param depth_thres Minimum depth trheshold.
#' @param freq_thres Freq threshold.
#'
#' @return A tibble that results of merging the info, dat and freq elements
#' of Dat. Plus allele assignments.
#' 
#' @export
#' @importFrom magrittr %>%
gene_midas_data <- function(Dat, gene, depth_thres = 1, freq_thres = 0.5){
  # gene <- genes[1]
  check_midas_data(Dat)
  
  # Get allele_freq and data from genes
  i <- Dat$info %>%
    dplyr::filter(gene_id == gene)
  f <- Dat$freq %>%
    dplyr::filter(site_id %in% i$site_id) %>%
    arrange(match(site_id, i$site_id))
  d <- Dat$depth %>%
    dplyr::filter(site_id %in% i$site_id) %>%
    arrange(match(site_id, i$site_id))
  all_sites <- match_freq_and_depth(freq = f, depth = d, info = i,
                                    map = NULL, depth_thres = depth_thres)
  all_sites <- all_sites %>%
    dplyr::select(-snp_type, -site_type, -amino_acids, -gene_id, - ref_id)
  
  
  all_sites <- assign_allele(all_sites,
                             depth_thres = depth_thres,
                             freq_thres = freq_thres,
                             sequence = TRUE, na_rm = TRUE)
  
  return(all_sites)
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
#' 
#' @examples 
#' library(HMVAR)
#' library(tidyverse)
#' 
#' # Get file paths
#' midas_dir <- system.file("toy_example/merged.snps/", package = "HMVAR")
#' map <- read_tsv(system.file("toy_example/map.txt", package = "HMVAR"),
#'                 col_types = cols(.default = col_character())) %>%
#'   select(sample = ID, Group)
#' 
#' # Read data
#' midas_data <- read_midas_data(midas_dir = midas_dir, map = map, cds_only = TRUE)
#' midas_data
read_midas_data <- function(midas_dir, map, genes = NULL, cds_only = TRUE){
  # Read data
  info <- readr::read_tsv(paste0(midas_dir, "/snps_info.txt"),
                          col_types = readr::cols(.default = readr::col_character(),
                                                  ref_pos = readr::col_number(),
                                                  count_samples = readr::col_number(),
                                                  count_a = readr::col_number(),
                                                  count_c = readr::col_number(),
                                                  count_g = readr::col_number(),
                                                  count_t = readr::col_number()),
                          na = 'NA')
  depth <- read_midas_abun(paste0(midas_dir, "/snps_depth.txt"))
  freq <- read_midas_abun(paste0(midas_dir, "/snps_freq.txt"))
  
  # Process data
  # Clean info
  info <- info %>% 
    dplyr::select(-tidyselect::starts_with("count_"))
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
      dplyr::filter(locus_type == 'CDS')
  }
  info <- info %>% dplyr::select(-locus_type)
  
  freq <- freq %>% 
    dplyr::filter(site_id %in% info$site_id)
  depth <- depth %>% 
    dplyr::filter(site_id %in% info$site_id)
  
  return(list(info = info, freq = freq, depth = depth))
}

#' Match freq and depth
#' 
#' Takes two data.frames or tibbles representing site x sample allele
#' frequencies (freq) and sequence coverage (depth), and used \link[tidyr]{gather}
#' to convert them into a tibble with one line per site per sample. It combines
#' both sources of info into one table, it filters out positions below a
#' sequencing depth threshold, and if another data table with metadata (info)
#' about the sites is passed, it also joins that information into the sites
#' that passed the threshold.
#'
#' @param freq A site by sample data frame or tibble. Must contain a column
#' named "site_id".
#' @param depth A site by sample data frame or tibble. Must containa a column
#' named "site_id".
#' @param info A site by variable data frame or tibble. Must contain a column
#' named "site_id".
#' @param map A sample by variable data frame or tibble. Must conatin a colum
#' named "sample".
#' @param depth_thres Minimum sequence coverage (depth) for a site to be kept
#' in the final output.
#'
#' @return A data frame or tibble with columns site_id, freq, depth, and one
#' column per column in info and map.
#' 
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' freq <- tibble::tibble(site_id = paste('snv' , 1:4, sep = ""),
#'                        sample1 = c(1,1,0,1), sample2 = c(1,1,1,1),
#'                        sample3=c(0,0,0,1))
#' depth <- tibble::tibble(site_id = paste('snv' , 1:4, sep = ""),
#'                         sample1 = c(1,0,1,1), sample2 = c(4,1,1,0),
#'                         sample3=c(0,0,0,1))
#' match_freq_and_depth(freq, depth, depth_thres = 1)
match_freq_and_depth <- function(freq, depth, info = NULL, map = NULL, depth_thres = 1){
  if(!("site_id" %in% colnames(freq))){
    stop("ERROR: freq mut contain a site_id column", call. = TRUE)
  }
  if(!("site_id" %in% colnames(depth))){
    stop("ERROR: depth mut contain a site_id column", call. = TRUE)
  }
  
  cat("\tReformatting data...\n")
  depth <- depth %>% tidyr::gather(key = "sample", value = 'depth', -site_id)
  freq <- freq %>% tidyr::gather(key = "sample", value = 'freq', -site_id)
  
  cat("\tMatching freq and depth...\n")
  dat <- depth %>%
    dplyr::inner_join(freq, by = c("site_id", "sample")) %>%
    dplyr::filter(depth >= depth_thres)
  
  if(!is.null(map)){
    cat("\tMatching map...\n")
    dat <- dat %>%
      dplyr::left_join(map, by = "sample")
  }
  
  if(!is.null(info)){
    cat("\tMatching info...\n")
    dat <- dat %>%
      dplyr::left_join(info, by = "site_id")
  }
  
  return(dat)
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
#' @importFrom magrittr %>%
midas_to_bimbam <- function(midas_dir, map, outdir, focal_group = NULL,
                            prefix = NULL){
  Dat <- read_midas_data(midas_dir = midas_dir,
                         map = map,
                         genes = NULL,
                         cds_only = FALSE)
  
  # Match freqs and depth
  Dat$depth <- Dat$depth %>% tidyr::gather(key = "sample", value = 'depth', -site_id)
  Dat$freq <- Dat$freq %>% tidyr::gather(key = "sample", value = 'freq', -site_id)
  Dat$info <- Dat$info %>% dplyr::select(site_id, ref_id, ref_pos, major_allele, minor_allele)
  
  # Set sites without coverage to NA
  dat <- Dat$depth %>%
    dplyr::inner_join(Dat$freq, by = c("site_id", "sample"))
  dat$freq[ dat$depth < 1 ] <- NA
  Dat$freq <- dat %>% dplyr::select(-depth) %>% tidyr::spread(sample, freq)
  
  # Create BIMBAM tables
  geno <- Dat$info %>%
    dplyr::select(site_id, minor_allele, major_allele) %>%
    dplyr::left_join(Dat$freq, by = "site_id")  
  
  pheno <- map %>%
    dplyr::filter(sample %in% colnames(geno)) %>%
    dplyr::arrange(factor(sample, levels = colnames(geno)[-(1:3)]))
  if(is.null(focal_group)){
    pheno <- pheno %>%
      dplyr::mutate(phenotype = as.numeric(factor(Group)) - 1) %>%
      dplyr::select(id = sample, phenotype)
  }else{
    pheno <- pheno %>%
      dplyr::mutate(phenotype = 1*(Group == focal_group)) %>%
      dplyr::select(id = sample, phenotype)
  }
  
  snp <- Dat$info %>%
    dplyr::select(ID = site_id, pos = ref_pos, chr = ref_id)
  # snp <- Dat$info %>% select(ID = site_id, pos = ref_pos, chr = ref_id)  %>%
  #   mutate(chr = as.numeric(factor(chr)))

  # Write bimbam tables
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  
  gen_file <- file.path(outdir, paste(c(prefix, 'geno.bimbam'), collapse = "_"))
  readr::write_tsv(geno, path = gen_file, col_names = FALSE)
  
  phen_file <- file.path(outdir, paste(c(prefix, 'pheno.bimbam'), collapse = "_"))
  readr::write_tsv(pheno %>%
                     dplyr::select(phenotype),
                   path = phen_file, col_names = FALSE)
  
  snp_file <- file.path(outdir, paste(c(prefix, 'snp.bimbam'), collapse = "_"))
  readr::write_tsv(snp, path = snp_file, col_names = FALSE)

  return(list(filenames = list(geno_file = gen_file,
                               pheno_file = phen_file,
                               snp_file = snp_file),
              Dat = list(geno = geno,
                         pheno = pheno,
                         snp = snp)))
}
