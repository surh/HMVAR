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

