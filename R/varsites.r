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


#' Determine DNA substitution type
#' 
#' Determines if a substitution was made of transitions or
#' transversions.
#'
#' @param info A data.frame or tibble. It should have columns
#' 'major_allele' and 'minor_allele'. These columns must have
#' single DNA characters in upper case.
#' @param clean Indicate if sites should be removed when
#' substitution type could not be correctly inferred.
#'
#' @return The same data.frame or tibble as info, with an extra
#' column called 'substitution'.
#' 
#' @export
#' 
#' @importFrom magrittr %>%
#'
#' @examples
#' library(HMVAR)
#' 
#' i <- dplyr::tibble(major_allele = c('A', 'C', 'G', 'T', 'a', NA),
#'                    minor_allele = c('T', 'T', 'A', 'A', 'T', 'G'))
#' determine_substitution_type(i, clean = FALSE)
determine_substitution_type <- function(info, clean = TRUE){
  info <- info %>%
    dplyr::bind_cols(substitution = info %>%
                       purrr:::pmap_chr(function(major_allele, minor_allele, ...){
                         base_type <- c(A = "purine", C = "pyrimidine", G = "purine", T = "pyrimidine")
                         major <- base_type[ major_allele ]
                         minor <- base_type[ minor_allele ]
                         if(any(is.na(c(major, minor)))){
                           substitution <- NA
                         }else if(major == minor){
                           substitution <- "transition" 
                         }else{
                           substitution <- "transversion"
                         }
                         return(substitution)
                       }))
  if(clean){
    info <- info %>%
      dplyr::filter(!is.na(substitution))
  }
  
  return(info)
}


#' Calculate prob of error
#' 
#' Internal function used to estimate prob that minor allele
#' at site is the product of error only.
#'
#' @param dat A data.frame or tibble, must have columns freq
#' and depth
#' @param error_prob Baseline err probability
#'
#' @return
#' 
#' @importFrom magrittr %>%
calculate_p_error <- function(dat, error_prob = 0.01){
  pvals <- dat %>%
    dplyr::mutate(x = purrr::pmap_dbl(list(freq, b = 1-freq), ~min(...))) %>%
    dplyr::transmute(x = round(x * depth),
                     n = depth) %>%
    purrr::pmap(binom.test, p = error_prob, alternative = "greater") %>%
    purrr::map_dbl(~ .x$p.value)
  
  return(pvals)
}

#' Determine distributiion within sample
#' 
#' For each site within each sample, it determines if
#' the dominant allele is homogeneous (likeley the only
#' allele present), or if there is evidence of a
#' second allele present (heterogeneous).
#'
#' @param dat A data.frame or tibble. One site per line,
#' must have columns depth and freq.
#' @param thres p-value threshold to consider that the evidence
#' against error only model is too large. Sites within samples
#' with p-value below threshold will be considered heterogeneous.
#' @param error_prob Baseline error probability in sequences.
#'
#' @return A data frame or tibble with columns for depth, freq
#' and sample_dist, as well as any other columns previously present.
#' @export
#' 
#' @importFrom magrittr %>%
determine_sample_dist <- function(dat, thres = 0.05, error_prob = 0.01){
  if(thres <= 0 || thres >= 1){
    stop("ERROR: thres must be in open interval (0,1)", call. = TRUE)
  }
  
  sample_dist <- calculate_p_error(dat, error_prob = error_prob) %>%
    cut(breaks = c(0, thres, 1), right = FALSE, include.lowest = TRUE)
  levels(sample_dist) <- c("heterogeneous", "homogeneous")
  
  dat <- dat %>% bind_cols(sample_dist = sample_dist)
  
  return(dat)
}


#' Size of groups
#' 
#' Internal function
#'
#' @param d A data frame or tibble
#' @param columns Names of columns to group
#'
#' @return A tibble
#' @importFrom magrittr %>%
group_size <- function(d, columns){
  d %>% split(.[,columns]) %>%
    purrr:::map_int(~nrow(.)) %>%
    t %>%
    tibble::as_tibble()
}

#' Calculated distribution of variable per site
#' 
#' Calculates the distribution of values of a categorical variable
#' per site from a table that contains one row per site per sample.
#'
#' @param dat A data frame or tibble containing columns "site_id",
#' "ref_id" and "ref_pos". Each row must correspond to a site per
#' sample.
#' @param variable Column name of variable to evaluate. It must be
#' a categorical variable.
#' @param group If passed, it must correspond to a column name in
#' dat. That column must be a grouping factor and the distribution
#' will be calculated independently for each group.
#'
#' @return A tibble with columns "site_id", "ref_id", and "ref_pos".
#' There will also be one column per level in `variable`, and,
#' optionally, one column for `group`.
#' 
#' @export
#' @importFrom magrittr %>%
#' 
#' @examples
#' library(magrittr)
#' map <- readr::read_tsv(system.file("toy_example/map.txt",
#'                                    package = "HMVAR"),
#'                        col_types = readr::cols(ID = readr::col_character(),
#'                                                Group = readr::col_character())) %>%
#'   dplyr::select(sample = ID,
#'                 tidyselect::everything())
#' Dat <- read_midas_data(midas_dir = system.file("toy_example/merged.snps/",
#'                                                package = "HMVAR"),
#'                        map = map,
#'                        cds_only = FALSE)
#' 
#' dat <- match_freq_and_depth(freq = Dat$freq,
#'                             depth = Dat$depth,
#'                             info = Dat$info,
#'                             map = map) %>%
#'   determine_sample_dist()
#' dat
variable_dist_per_site <- function(dat, variable, group = NULL){
  if(!all(c("site_id", "ref_id", "ref_pos", group, variable) %in% colnames(dat))){
    stop("ERROR: missing columns in dat")
  }
  
  if(is.null(group)){
    res <- dat %>%
      split(.$site_id) %>%
      purrr::map_dfr(group_size, columns = 'sample_dist', .id = 'site_id')
  }else{
    res <- dat %>%
      split(.$site_id) %>%
      purrr::map_dfr(function(d, column){
        d %>% split(.[,column]) %>%
          purrr::map_dfr(group_size, columns = 'sample_dist', .id = column)},
        column = 'Group',
        .id = 'site_id')
  }
  
  dat %>%
    dplyr::select(site_id, ref_id, ref_pos) %>%
    dplyr::filter(!duplicated(.)) %>%
    dplyr::full_join(res, by = "site_id")
}