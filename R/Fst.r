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
# along with HMVAR.  If not, see <http://www.gnu.org/licenses/>


#' Fixation Index (Fst) for one site
#' 
#' Takes allele frequency data from one site only
#' and calculates Wright's fixation index (Fst).
#' 
#' It works only for haploid individuals (e.g proportion
#' of heterozygous is 0).
#'
#' @param freq A 1 x (n+1) tibble. One of the columns must correspond
#' to 'site_id', and the rest of the columns correspond to one sample
#' each. Entries under each sample correspond to allele frequencies.
#' @param support A 1 x (n+1) tibble. Columns must match columns in freq.
#' Values must be numeric and indicate support of each allele frequencuy
#' @param info A tiblle with one row. One column must correspond to
#' 'site_id'. The rest of the columns contain information about the site
#' @param map A tibble with columns 'sample' (corresponding to column names)
#' in freq and support, and 'Group' corresponding to the populations used
#' for Fst calculations.
#' @param support_thres Minimum value in support tibble to keep a sample
#' for this site
#' @param method The method to calculate Fst. Currently only the
#' Weir-Cockerham method from 1984 is implemented.
#'
#' @return A tibble
#' @export
#' 
#' @importFrom magrittr %>%
site_fst <- function(freq, support, info,
                     map, support_thres = 1,
                     method = "Weir-Cockerham"){
  dat <- match_freq_and_depth(freq = freq,
                              depth = support,
                              info = info,
                              map = map,
                              depth_thres = support_thres,
                              verbose = FALSE)
  
  # Using Weir-Cockerham 1984 method
  # Get basic quantities
  n_i <- table(dat$Group)
  p_i <- dat %>% split(.$Group) %>% purrr::map_dbl(~mean(.$freq))
  p_i <- p_i[names(n_i)]
  r <- length(n_i)
  n_mean <- nrow(dat) / r
  n_c <- ((r * n_mean) - (sum(n_i ^ 2) / (r * n_mean))) / (r - 1)
  # p_mean <- sum(n_i * p_i) / (r * n_mean)
  p_mean <- mean(dat$freq)
  S_sqrd <- sum(n_i * ((p_i - p_mean) ^ 2)) / ((r-1) * n_mean)
  h_mean <- 0
  
  # Calculate parts
  a <- (n_mean / n_c) * (S_sqrd - (1 / (n_mean - 1)) * ((p_mean * (1 - p_mean)) - ((r - 1) * S_sqrd / r) - (h_mean / 4)))
  b <- (n_mean / (n_mean - 1)) * (p_mean * (1 - p_mean) - ((r - 1) * S_sqrd / r) - ((2 * n_mean - 1) * h_mean / (4 * n_mean)) )
  c <- h_mean / 2
  
  tibble::tibble(r = r, n_mean = n_mean, n_c = n_c,
                 p_mean = p_mean, S_sqrd = S_sqrd,
                 h_mean = h_mean,
                 a = a, b = b,  c = c)
  
}

# Tidyverse based function removed because it takes too much time.
# Probably if I join first it would be faster
# calculate_fst <- function(sites, Dat, map, depth_thres = 1, method = "Weir-Cockerham"){
#   res <- sites %>%
#     purrr::map_dfr(function(site, Dat, map, depth_thres = 1, method = "Weir-Cockerham"){
#       dat <- Dat %>%
#         purrr::map(function(d, site){d %>% filter(site_id == site)}, site = site)
#       site_fst(freq = dat$freq, depth = dat$depth,
#                info = dat$info, map = map, depth_thres = depth_thres)
#     }, Dat = Dat, map = map, depth_thres = 1, method = method, .id = "site_id")
#   
#   return(res)
# }
