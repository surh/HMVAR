# (C) Copyright 2018-2019 Sur Herrera Paredes
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

#' Get files for MK results
#' 
#' Takes a directory and finds all files that
#' correspond to results of MKtests.py, and returns
#' a list of full paths for those files
#' 
#' @param d Directory
#' @param pattern Pattern to match for files
#' containing results from MKtests.py
#' 
#' @return Vector of files
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
get_mk_results_files <- function(d, pattern = "^mk_results"){
  
  files <- list.files(d)
  chosen <- files[ grep(pattern = pattern, x = files) ]
  
  if(length(chosen) > 0)
    chosen <- paste(d, chosen, sep = "/")
  
  return(chosen)
}

#' Check p-values in MKtest.py outfile
#' 
#' Opens a file from MKtest.py and uses the qvalue methods
#' of Storey and Tibshirani to estimate the proportion
#' of True Negatives (Pi0) in the data
#' 
#' @param file file name
#' @param which either 'all' or the name of the test to
#' be analyzed
#' @param plot whether to return a histogram of the p-values
#' 
#' @author Sur Herrera Paredes
#' 
#' @return  A named list where each element corresponds to
#' a test and each element is alist with pi0 and p1 entries.
#' See \link{check_pvalues} for more info
#' 
#' @export
check_pvals_in_file <- function(file, which, plot=TRUE){
  Tab <- read.table(file, header = TRUE, sep = "\t")
  default <- c("gene", "contig", "start", "end", "Dn", "Ds", "Pn", "Ps")
  
  if(nrow(Tab) == 0)
    return(list())
  
  # Get pval columns
  if(which == "all"){
    tests <- setdiff(colnames(Tab), default)
    tests <- tests[ grep(pattern = ".pval", x = tests, invert = TRUE) ]
    tests <- tests[ grep(pattern = ".perm", x = tests, invert = TRUE) ]
  }else{
    tests <- which
  }
  
  Res <- list()
  for(t in tests){
    # t <- tests[8]
    res <- check_pvalues(Tab[,t], Tab[,paste(t,".pval",sep="")], plot = plot)
    Res[[t]] <- res
  }
  
  return(Res)
}

#' Check p-values
#' 
#' Takes a pair of vectors of estimates and their p-values
#' for some hypothesis test, and uses qvalue methods from
#' storey and Tibshirani to estimate the number of True
#' Negatives (pi0) and creat a p-value histogram. It doesn't
#' consider cases where the estimate is undefined.
#' 
#' @param estimates numeric verctor of estimates
#' @param pvals numeric vector of p-values
#' @param plot plot or not
#' 
#' @return a list with elements pi0 with the point estimate
#' of pi0 and p1 which is a ggplot2 plot or NULL in case parameter
#' plot was FALSE
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
check_pvalues <- function(estimates, pvals, plot = TRUE){
  pvals <- pvals[ !is.na(estimates) ]
  qvals <- qvalue::qvalue(pvals)
  # qvals.sum <- summary(qvals)
  pi0 <- qvals$pi0
  
  p1 <- NULL
  if(plot){
    p1 <- ggplot(data.frame(pvals),aes(x = pvals)) +
      geom_histogram(bins = 20) +
      AMOR::theme_blackbox
  }
  
  res <- list(pi0 = pi0, p1 = p1)
  
  return(res)
}


############# MIDAS MKTEST ###################
# Code for obtaining MKtest from midas merge output

#' Determine the effect of a coding variant on the aminoacid sequence
#' 
#' Takes an table corresponding to the contents of the snp_info.txt
#' file from midas_merge.py and adds a column indicating whether
#' the variant is synonymous or non-synonymous
#'
#' @param info A data table corresponding to the contents of
#' the snp_info.txt file produced by midas_merge. It must have columns:
#' 'site_id', 'major_allele', 'minor_allele' and 'amino_acids'. The aminoacid
#' column must contain a string of four comma-separated values indicating
#' the aminoacid encoded for each variant (eg. 'V,L,V,L').
#' @param nucleotides Named vector indicating the position of each
#' nucleotide that corresponds to the amino_acids column in info. The
#' default corresponds to MIDAS v1.3.1 default.
#'
#' @return The same data table passed as info with a factor column
#' 'snp_effect' added.
#' @export
#'
#' @importFrom magrittr %>%
#' 
#' @examples
#' library(HMVAR)
#' 
#' # Get file paths
#' midas_dir <- system.file("toy_example/merged.snps/", package = "HMVAR")
#' map <- readr::read_tsv(system.file("toy_example/map.txt", package = "HMVAR"),
#'                        col_types = readr::cols(.default = readr::col_character())) %>%
#'   dplyr::select(sample = ID, Group)
#' 
#' # Read data
#' midas_data <- read_midas_data(midas_dir = midas_dir, map = map, cds_only = TRUE)
#' 
#' info <- determine_snp_effect(midas_data$info) %>%
#'   determine_snp_dist(freq = midas_data$freq,
#'                      depth = midas_data$depth, map = map,
#'                      depth_thres = 1, freq_thres = 0.5)
#' info
#' 
#' mktable <- info %>%
#'   split(.$gene_id) %>%
#'   purrr::map_dfr(mkvalues,
#'                  .id = "gene_id")
#' mktable
determine_snp_effect <- function(info, nucleotides=c(A = 1, C = 2, G = 3, T = 4)){
  snp_effect <- info %>%
    dplyr::select(site_id, major_allele, minor_allele, amino_acids) %>%
    purrr::pmap_chr(function(site_id, major_allele, minor_allele,
                             amino_acids, nucleotides){
      aa = stringr::str_split(string = amino_acids,
                              pattern = ",",
                              simplify = TRUE)
      if(all(is.na(aa))){
        type <- NA
      }else if( aa[nucleotides[major_allele]] == aa[nucleotides[minor_allele]] ){
        type <- "synonymous"
      }else{
        type <- "non-synonymous"
      }
      return(type)},
      nucleotides = nucleotides)
  
  # Add column
  info <- info %>%
    tibble::add_column(snp_effect = factor(snp_effect,
                                           levels = c('synonymous',
                                                      'non-synonymous')))
  return(info)
}

#' Determine whether a variant site is fixed or polymorphic
#' 
#' This is an internal utility function
#' 
#' Looks at the distribution of alleles at one site and
#' determines how it distributes with respect to a grouping factor.
#' If the site perfectly separates both groups it is 'Fixed'
#' and 'Polymorphic' otherwise.
#' 
#' Assumes grouping factor with two levels (NO CHECKS!!)
#' 
#' Assumes bi-allelic site indicated by character string (NO CHECKS!!)
#' 
#' Assumes only data from one site is give (NO CHECKS!!)
#' 
#' USE WITH CARE!!
#'
#' @param d A data table that matches each site at each
#' sample with the group of that sample. Must have column
#' 'allele' which must have two values only; column 'Group'
#' which is the grouping factor (can be character) with two
#' levels.
#' @param group_thres The minimum number of samples per
#' group to consider the site.
#'
#' @return A character string indicating the type of site:
#' 'Fixed', 'Polymorphic', 'Invariant' or NA.
get_site_dist <- function(d, group_thres = 2){
  groups <- split(d$allele, d$group)
  n1 <- NA
  n2 <- NA
  if(length(groups) == 1){
    dist <- NA
  }else{
    if(length(groups[[1]]) < group_thres || length(groups[[2]]) < group_thres){
      dist <- NA
    }else{
      g1 <- unique(groups[[1]])
      g2 <- unique(groups[[2]])
      
      n1 <- length(groups[[1]])
      n2 <- length(groups[[2]])
      
      if(length(g1) > 1 || length(g2) > 1){
        dist <- 'Polymorphic'
      }else if(g1 == g2){
        dist <- 'Invariant'
      }else if(g1 != g2){
        dist <- 'Fixed'
      }else{
        dist <- 'What!!!'
      }
    }
  }
  # return(dist)
  tibble::tibble(snp_dist = dist,
                 n1 = n1,
                 n2 = n2)
}


# Alternative function to determine site distribution
# determine_site_dist <- function(d){
#   # dat <- subset(dat, site_id == "703112")
#   
#   counts <- ftable(Group ~ allele, d)
#   # counts
#   
#   d1 <- diag(counts)
#   d2 <- counts[row(counts) - col(counts) != 0]
#   if(all(d1 > 0) && all(d2 == 0)){
#     type <- "fixed"
#   }else if(all(d1 == 0) && all(d2 > 1)){
#     type <- "fixed"
#   }else if(any(colSums(counts) == 0)){
#     type <- NA
# 
#   }else{
#     type <- "polymorphic"
#   }
#   # return(type)
#   return(tibble(site_id = d$site_id[1],
#                 type = type))
# }
# 
# system.time(determine_site_dist(d))
# T1 <- NULL
# T2 <- NULL
# for(i in 1:100){
#   t1 <- system.time(fun2(d))
#   t2 <- system.time(determine_site_dist(d))
#   
#   T1 <- rbind(T1, t1)
#   T2 <- rbind(T2, t2)
# }

#' SNP distribution between sites
#' 
#' Determines how snps distribute between sites. Requires
#' output from midas_merge.py and a mapping file mapping 
#' samples to sites.
#' 
#' Only samples in both the map and the depth and freq tables
#' are considered. Everything else is removed (inner_join)
#'
#' @param info Data table corresponding to the 'snps_info.txt'
#' file from MIDAS. Must have columns 'site_id' and 'sample'
#' @param freq A data table corresponding to the 'snps_freq.txt'
#' file from MIDAS. Must have a 'site_id' column, and one more
#' column per sample. Each row is the frequency of the minor
#' allele for the corresponding site in the corresponding sample.
#' @param depth A data table corresponding to the 'snps_depth.txt'
#' file from MIDAS. Must have a 'site_id' column, and one more
#' column per sample. Each row is the sequencing depth for the
#' corresponding site in the corresponding sample.
#' @param map A data table associating samples with groups (sites).
#' must have columns 'sample' and 'Group'.
#' @param depth_thres Minimum number of reads (depth) at a site at
#' a sample to be considered.
#' @param freq_thres Frequency cuttoff for minor vs major allele.
#' The value represents the distance from 0 or 1, for a site to be
#' assigned to the major or minor allele respectively. It must be
#' a value in [0,1].
#' @param clean Whether to remove sites that had no valid
#' distribution.
#'
#' @return A data table which is the same and info bnut with
#' a 'distribution' column indicating the allele distribution
#' between sites in the  given samples.
#' @export
#' 
#' @importFrom magrittr %>%
#' 
#' @examples
#' library(HMVAR)
#' 
#' # Get file paths
#' midas_dir <- system.file("toy_example/merged.snps/", package = "HMVAR")
#' map <- readr::read_tsv(system.file("toy_example/map.txt", package = "HMVAR"),
#'                        col_types = readr::cols(.default = readr::col_character())) %>%
#'   dplyr::select(sample = ID, Group)
#' 
#' # Read data
#' midas_data <- read_midas_data(midas_dir = midas_dir, map = map, cds_only = TRUE)
#' 
#' info <- determine_snp_effect(midas_data$info) %>%
#'   determine_snp_dist(freq = midas_data$freq,
#'                      depth = midas_data$depth, map = map,
#'                      depth_thres = 1, freq_thres = 0.5)
#' info
#' 
#' mktable <- info %>%
#'   split(.$gene_id) %>%
#'   purrr::map_dfr(mkvalues,
#'                  .id = "gene_id")
#' mktable
determine_snp_dist <- function(info, freq, depth, map,
                               depth_thres = 1,
                               freq_thres = 0.5,
                               clean = TRUE){
  
  # Process freq_thres
  if(freq_thres < 0 || freq_thres > 1)
    stop("ERROR: freq_thres must have values in [0, 1]", call. = TRUE)
  if(ncol(freq) <= 1 || ncol(depth) <= 1){
    rlang::warn("\tWARNING:There are no samples in freq or depth. Returning empty tibble")
    info <- info %>% dplyr::filter(is.na(site_id) & !is.na(site_id)) %>%
      bind_cols(distribution = factor(levels = c('Fixed', 'Invariant', 'Polymorphic')))
    return(info)
  }
  
  freq_thres <- min(freq_thres, 1 - freq_thres)
  
  # Reformat
  depth <- depth %>% tidyr::gather(key = "sample", value = 'depth', -site_id)
  freq <- freq %>% tidyr::gather(key = "sample", value = 'freq', -site_id)
  
  # Last lines can be re-written for speed!!
  dat <- depth %>%
    dplyr::inner_join(freq, by = c("site_id", "sample")) %>%
    dplyr::left_join(map, by = "sample") %>%
    dplyr::filter(depth >= depth_thres) %>%
    dplyr::filter(freq != 0.5) %>%
    dplyr::mutate(allele = replace(freq, freq <= freq_thres, 'major')) %>%
    dplyr::mutate(allele = replace(allele, freq >= (1 - freq_thres), 'minor')) %>%
    dplyr::mutate(allele = replace(allele, (freq > freq_thres) & (freq < (1 - freq_thres)), NA)) %>%
    dplyr::filter(!is.na(allele))
  
  site_dist <- dat %>%
    split(.$site_id) %>%
    purrr::map_chr(get_site_dist)
  site_dist <- tibble::tibble(site_id = names(site_dist),
                              distribution = factor(site_dist,
                                                    levels = c('Fixed', 'Invariant', 'Polymorphic')))
  info <- info %>%
    dplyr::left_join(site_dist, by = "site_id")
  
  if(clean){
    info <- info %>%
      dplyr::filter(!is.na(distribution))
  }
  
  return(info)
}

#' Entries in McDonald-Kreitman table.
#' 
#' Calculates the four entries (Dn, Ds, Pn, Ps)
#' in the McDonald-Kreitman contingency table.
#'
#' @param info A data table corresponding to the
#' MIDAS snps_info.txt output file. Must have been
#' processed with \link{determine_snp_effect} and
#' \link{determine_snp_dist}, and so it must contain
#' columns 'snp_effect' and 'distribution'
#'
#' @return A tibble with one column with the count
#' of each substitution type.
#' @export 
#' 
#' @importFrom magrittr %>%
#' 
#' @examples
#' library(HMVAR)
#' 
#' # Get file paths
#' midas_dir <- system.file("toy_example/merged.snps/", package = "HMVAR")
#' map <- readr::read_tsv(system.file("toy_example/map.txt", package = "HMVAR"),
#'                        col_types = readr::cols(.default = readr::col_character())) %>%
#'   dplyr::select(sample = ID, Group)
#' 
#' # Read data
#' midas_data <- read_midas_data(midas_dir = midas_dir, map = map, cds_only = TRUE)
#' 
#' info <- determine_snp_effect(midas_data$info) %>%
#'   determine_snp_dist(freq = midas_data$freq,
#'                      depth = midas_data$depth, map = map,
#'                      depth_thres = 1, freq_thres = 0.5)
#' info
#' 
#' mktable <- info %>%
#'   split(.$gene_id) %>%
#'   purrr::map_dfr(mkvalues,
#'                  .id = "gene_id")
#' mktable
mkvalues <- function(info){
  if(nrow(info) < 1){
    rlang::warn("\tWARNING: No SNPs for mkvalues. Returning empty tibble.\n")
    return( tibble::tibble(Dn = integer(0), Ds = integer(0),
                           Pn = integer(0), Ps = integer(0)))
  }else{
    tab <- info %>% 
      dplyr::select(snp_effect, distribution) %>%
      table(exclude = NULL, useNA = 'always')
    
    return(dplyr::tibble(Dn = tab['non-synonymous', 'Fixed'],
                         Ds = tab['synonymous', 'Fixed'],
                         Pn = tab['non-synonymous', 'Polymorphic'],
                         Ps = tab['synonymous', 'Polymorphic']))
  }
}

#' Perform McDonald-Kreitman test on MIDAS SNPs
#' 
#' Take output from MIDAS script midas_merge.py
#' and perform a McDonald-Kreitman test between
#' two groups.
#'
#' @param midas_dir Directory where the output from
#' midas_merge.py is located. It must include files:
#' 'snps_info.txt', 'snps_depth.txt' and 'snps_freq.txt'.
#' @param map Either  a file path or a tibble. If a path,
#' it mos point to a mapping file associating samples
#' in the MIDAS otput to groups. It must have an 'ID'
#' and a 'Group' column. If a tibble. It must have columns
#' 'sample' and 'Group'.
#' @param map_file Same as map. Present for backwards compatibility.
#' @param genes The list of genes that are to be tested.
#' Must correspond to entries in the 'genes_id' column
#' of the 'snps_info.txt' file. If NULL, all genes will
#' be tested.
#' @param depth_thres The minimum number of reads at
#' a position in a given sample for that position in that
#' sample to be included.
#' @param freq_thres Frequency cuttoff for minor vs major allele.
#' The value represents the distance from 0 or 1, for a site to be
#' assigned to the major or minor allele respectively. It must be
#' a value in [0,1].
#' @param focal_group A string indicating the group that should
#' be compared to everything else. If different from NA, all values
#' in the 'Group' column of the mapping file will be converted to
#' "non.<focal_group>", ensuring a dichotomous grouping.
#'
#' @return A data table containing the McDonald-Kreitman
#' contingency table per gene.
#' 
#' @export
#' 
#' @importFrom magrittr %>%
#' 
#' @examples 
#' library(HMVAR)
#' 
#' # Get paths
#' midas_dir <- system.file("toy_example/merged.snps/", package = "HMVAR")
#' map_file <- system.file("toy_example/map.txt", package = "HMVAR")
#' 
#' # Process map yourself
#' map <- readr::read_tsv(map_file,
#'                        col_types = readr::cols(.default = readr::col_character())) %>%
#'   dplyr::select(sample = ID, Group)
#' midas_mktest(midas_dir = midas_dir,
#'              map = map)
#' 
#' # Give a map file path to midas_mktest
#' midas_mktest(midas_dir = midas_dir,
#'              map = map_file)
midas_mktest <- function(midas_dir, map,
                         map_file,
                         genes = NULL,
                         depth_thres = 1,
                         freq_thres = 0.5,
                         focal_group = NA){
  
  map <- check_map(map = map,
                   map_file = map_file,
                   focal_group = focal_group)
  
  # Read data
  Dat <- read_midas_data(midas_dir = midas_dir,
                         map = map,
                         genes = genes)
  
  # ML table
  Res <- calculate_mktable(info = Dat$info,
                           freq = Dat$freq,
                           depth = Dat$depth,
                           map = map,
                           depth_thres = depth_thres,
                           freq_thres = freq_thres)
  
  return(Res)
}

#' Calculate McDonald-Kreitman contingency table
#' 
#' Takes SNVs in MIDAS format and produces a McDonald-Kreitman
#' contingency table.
#'
#' @param info A snv x variable data frame or table. Must have columns
#' 'site_id', 'minor_allele', 'major_allele', 'gene_id' and 'amino_acids'
#' as defined in the midas_merge.py snps from MIDAS.
#' @param freq A snv x sample data frame or tibble with minor allele
#' frequenncies. Must have a 'site_id' column.
#' @param depth A snv x sample data frame or tibble with sequence coverage.
#' Must have a 'site_id' column.
#' @param map A sample x group data frame or tibble. Must have columns
#' 'sample' and 'Group'.
#' @param depth_thres The minimum number of reads at
#' a position in a given sample for that position in that
#' sample to be included.
#' @param freq_thres Frequency cuttoff for minor vs major allele.
#' The value represents the distance from 0 or 1, for a site to be
#' assigned to the major or minor allele respectively. It must be
#' a value in [0,1].
#'
#' @return A tibble with columns gene_id, Dn, Ds, Pn, and Ps, which
#' correspond to the McDonald-Kreitman contingency table.
#' 
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' library(tidyverse)
#' 
#' map <- read_tsv(system.file("toy_example/map.txt", package = "HMVAR"),
#'                 col_types = cols(.default = col_character())) %>%
#'   select(sample = ID, Group)
#' 
#' midas_data <- read_midas_data(system.file("toy_example/merged.snps/", package = "HMVAR"),
#'                               map = map,
#'                               genes = NULL,
#'                               cds_only = TRUE)
#' 
#' calculate_mktable(info = midas_data$info, freq = midas_data$freq, depth = midas_data$depth, map = map)
calculate_mktable <- function(info, freq, depth, map, depth_thres = 1, freq_thres = 0.5){
  # Calcualate snp effect
  info <- determine_snp_effect(info)
  
  # Calculate snp dist
  info <- determine_snp_dist(info = info,
                             freq = freq,
                             depth = depth,
                             map = map,
                             depth_thres = depth_thres,
                             freq_thres = freq_thres)
  
  Res <- info %>%
    split(.$gene_id) %>%
    purrr::map_dfr(mkvalues,
                   .id = "gene_id") %>%
    dplyr::bind_rows(tibble::tibble(gene_id = character(0),
                                    Dn = integer(0), Ds = integer(0),
                                    Pn = integer(0), Ps = integer(0)))
  
  return(Res)
}

#' McDonald-Kreitman test
#'
#' @param alleles site x sample alleles table. First column must be
#' "site_id"
#' @param info Table with information about SNVs. Must include columns
#' "site_id", "locus_type", "gene_id" and "snp_effect".
#' @param map Mapping file. must have columns "sample" and "group"
#'
#' @return A tibble
#' @export
#' @importFrom magrittr %>%
#' 
mktest <- function(alleles, info, map){
  if(!("site_id" %in% colnames(alleles))){
    stop("ERROR: alleles table needs column site_id", call. = TRUE)
  }
  if(!all(c("site_id", "locus_type", "gene_id", "snp_effect") %in% colnames(info))){
    stop("ERROR: info table needs columns 'site_id', 'locus_type', 'gene_id', and 'snp_effect'",
         call. = TRUE)
  }
  if(!all(c("sample", "group") %in% colnames(map))){
    stop("ERROR: map must have columns sample and group.", call. = TRUE)
  }
  if(any(alleles$site_id != info$site_id)){
    stop("ERROR: there are inconsistent entries between alleles and info tables",
         call. = TRUE)
  }
  
  # Get group variable per sample
  sample_ids <- colnames(alleles)[-1]
  meta <- map %>%
    dplyr::filter(sample %in% sample_ids)
  
  if(length(setdiff(sample_ids, meta$sample)) > 0){
    stop("ERROR: Thera are samples that do not appear in the map", call. = TRUE)
  }
  
  Res <- info %>%
    dplyr::filter(locus_type == "CDS") %>%
    dplyr::select(site_id, gene_id, snp_effect) %>%
    dplyr::inner_join(alleles, by = "site_id") %>%
    split(.$gene_id) %>%
    purrr::map_dfr(function(d, meta){
      # For each gene, determine if each SNV is fixed or
      # polymorphic
      snv_dist <- d %>%
        tidyr::pivot_longer(c(-site_id, -gene_id, -snp_effect),
                            names_to = "sample",
                            values_to = "allele") %>%
        dplyr::filter(!is.na(allele)) %>%
        dplyr::left_join(meta %>%
                           dplyr::select(sample, group),
                         by = "sample") %>%
        split(.$site_id) %>%
        # purrr::map_chr(HMVAR:::get_site_dist)
        purrr::map_dfr(HMVAR:::get_site_dist,
                       .id = "site_id")
        
      
      # Combine results
      d <- d %>%
        dplyr::select(site_id, snp_effect) %>%
        # dplyr::left_join(tibble::tibble(site_id = names(snv_dist),
        #                                 snp_dist = as.character(snv_dist)),
        #           by = "site_id") %>%
        dplyr::left_join(snv_dist,
                  by = "site_id") %>%
        dplyr::filter(!is.na(snp_dist))
      
      # Fill MK table
      tibble::tibble(med.n1 = median(d$n1),
                     med.n2 = median(d$n2),
                     Dn = sum(d$snp_effect == "non-synonymous" & d$snp_dist == "Fixed"),
                     Ds = sum(d$snp_effect == "synonymous" & d$snp_dist == "Fixed"),
                     Pn = sum(d$snp_effect == "non-synonymous" & d$snp_dist == "Polymorphic"),
                     Ps = sum(d$snp_effect == "synonymous" & d$snp_dist == "Polymorphic"))
    }, meta = meta, .id = "gene_id") %>%
    dplyr::filter(Ds > 0 & Pn > 0) %>%
    purrr::pmap_dfr(function(gene_id, med.n1, med.n2, Dn, Ds, Pn, Ps){
      # Make test
      test <- fisher.test(matrix(c(Dn, Ds, Pn, Ps), ncol = 2))
      tibble::tibble(gene_id = gene_id,
                     med.n1 = med.n1,
                     med.n2 = med.n2,
                     Dn = Dn,
                     Ds = Ds,
                     Pn = Pn,
                     Ps = Ps,
                     OR = test$estimate,
                     p.value = test$p.value)
    })
  
  Res
}



