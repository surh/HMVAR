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

#' Assign alleles bases on allele frequency and depth
#' 
#' Internal
#'
#' @param dat A data.frame or tibble
#' @param depth_thres depth_threshold
#' @param freq_thres freq_threshold
#' @param sequence Should the actual nucleodites be returned? If FALSE
#' the allele column will contain 'major' or 'minor'. If TRUE, the
#' allele column will contain the actual nucleotide
#' @param na_rm Should sites from samples that could not be assigned be returned?
#'
#' @return Same tibble as dat but with an allele column
#' 
#' @importFrom magrittr %>%
assign_allele <- function(dat, depth_thres = 1, freq_thres = 0.5, sequence = FALSE, na_rm = TRUE){
  if(!is.data.frame(dat))
    stop("ERROR: dat must be a data.frame.", call. = TRUE)
  if(!all(c('depth', 'freq', 'major_allele', 'minor_allele') %in% colnames(dat)))
    stop("ERROR: dat must have columns c('depth', 'freq', 'major_allele', 'minor_allele').", call. = TRUE)
  
  # Assign allele
  dat <- dat %>%
    dplyr::filter(depth >= depth_thres) %>%
    dplyr::filter(freq != 0.5) %>%
    dplyr::mutate(allele = replace(freq, freq <= freq_thres, 'major')) %>%
    dplyr::mutate(allele = replace(allele, freq >= (1 - freq_thres), 'minor')) %>%
    dplyr::mutate(allele = replace(allele, (freq > freq_thres) & (freq < (1 - freq_thres)), NA))
  
  if(na_rm){
    dat <- dat %>%
      dplyr::filter(!is.na(allele))
  }
  
  if(sequence){
    dat$allele[ dat$allele == 'major' ] <- dat$major_allele[ dat$allele == 'major' ]
    dat$allele[ dat$allele == 'minor' ] <- dat$minor_allele[ dat$allele == 'minor' ]
  }
  
  return(dat)
}
