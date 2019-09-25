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

#' Calculate Fixation Index (Fst)
#' 
#' Calculates Fst from imported midas data
#'
#' @param Dat A MIDAS dataset. Imported into R with \link{read_midas_data}
#' @param map A tibble with columns sample and Group corresponding
#' to sample names in MIDAS dataset and sample population respectiveley.
#' @param method Weir-Cockerham or Fstpool.
#' @param support_thres calculate_fst uses the \code{depth} slot of
#' the MIDAS data as a 'support' value to decide if a site in a given
#' sample will be considered. If it corresponds to depth, this is the
#' minimum coverage at a site at a given sample.
#' @param w_size Window size in bp. If both w_size and s_size
#' are numeric. Fst/Fst^{pool} values are calculated for a sliding
#' window.
#' @param s_size Window step size in bp. If both w_size and s_size
#' are numeric. Fst/Fst^{pool} values are calculated for a sliding
#' window.
#' @param sorted Logical indicating if sites are sorted in ascending
#' orderby ref_pos within each ref_id. Only used if sliding windows are calculated
#' @param verbose Logical indicating if progress during site level calculations
#' must be printed.
#'
#' @return A list with tibbles fst and w_fst corresponding to site
#' and window level Fst calculations. w_fst is NULL if no window
#' calculations are performed.
#' 
#' @export
#' @importFrom magrittr %>%
calculate_fst <- function(Dat, map,
                          method = "Weir-Cockerham",
                          support_thres = 1,
                          w_size = NULL,
                          s_size = NULL,
                          sorted = FALSE,
                          verbose = TRUE){
  fst <- NULL
  for(i in 1:nrow(Dat$info)){
    f <- site_fst(freq = Dat$freq[i,],
                  support = Dat$depth[i,],
                  info = Dat$info[i,],
                  map = map,
                  support_thres = support_thres,
                  method = method)
    
    f$site_id <- Dat$info$site_id[i]
    f$ref_id <- Dat$info$ref_id[i]
    f$ref_pos <- Dat$info$ref_pos[i]
    fst <- fst %>%
      dplyr::bind_rows(f)
    
    if(verbose && ((i %% 1000) == 0))
      cat(i, "\n")
  }
  
  if(method == "Weir-Cockerham"){
    if(verbose)
      cat("Calculating Fst.\n")
    fst <- fst %>%
      dplyr::mutate(Fst = a / (a + b + c)) %>%
      dplyr::mutate(Fst = replace(Fst, Fst < 0, 0))
  }else if(method == "Fstpool"){
    if(verbose)
      cat("Calculating Fst_pool.\n")
    fst <- fst %>%
      dplyr::mutate(Fst_pool = (MSP - MSI) / (MSP + ((n_c - 1) * MSI)) ) %>%
      dplyr::mutate(Fst_pool = replace(Fst_pool, Fst_pool < 0, 0))
  }
  
  w_fst <- NULL
  if(is.numeric(w_size) && is.numeric(s_size)){
    if(verbose)
      cat("Calculating multilocus Fst across sliding window.\n")
    w_fst <- fst %>%
      split(.$ref_id) %>%
      purrr::map_dfr(ref_window_fst, w_size = 1000, s_size = 300, sorted = sorted) 
  }
  
  return(list(fst = fst, w_fst = w_fst))
}

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
#' 'Weir-Cockerham' method from 1984, and the 'Fstpool' methods are implemented.
#'
#' @return A tibble
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
  
  if(method == "Weir-Cockerham"){
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
    
    res <- tibble::tibble(r = r, n_mean = n_mean, n_c = n_c,
                          p_mean = p_mean, S_sqrd = S_sqrd,
                          h_mean = h_mean,
                          a = a, b = b,  c = c)
    # res <- res %>%
    #   dplyr::mutate(Fst = a / (a + b + c))
  }else if(method == "Fstpool"){
    # Make pools
    dat <- dat %>%
      split(.$Group) %>%
      purrr::map_dfr(function(d){
        tibble::tibble(#site_id = unique(d$site_id),
                       depth = sum(d$depth),
                       freq = sum(d$depth * d$freq) / sum(d$depth),
                       count = sum(round(d$freq * d$depth)),
                       n_ind = nrow(d))
      }, .id = "Group")
    
    # Parameters
    C_1 <- sum(dat$depth)
    C_2 <- sum(dat$depth ^ 2)
    D_2 <- sum((dat$depth / dat$n_ind) + ((dat$n_ind - 1) / dat$n_ind))
    D_2ast <- sum(dat$depth * (dat$depth + dat$n_ind + 1) / dat$n_ind) / C_1 
    pi_k <- sum(dat$depth * dat$freq) / sum(dat$depth)
    n_c <- (C_1 - C_2/ C_1) / (D_2 - D_2ast)
    
    # MSI assume only bi-allelic sites
    MSI <- (1 / (C_1 - D_2)) * 2 * sum(dat$depth * dat$freq * (1 - dat$freq))
    
    # MSP assume bi-allelic sites
    MSP <- (1 / (D_2 - D_2ast)) * 2 * sum(dat$depth * (dat$freq - pi_k) ^ 2)
    
    
    res <- tibble::tibble(r = nrow(dat),
                          n_c = n_c,
                          p_mean = pi_k,
                          MSI = MSI,
                          MSP = MSP)
    # res <- res %>%
    #   dplyr::mutate(Fst_pool = (MSP - MSI) / (MSP + ((n_c - 1) * MSI)) )
  }else{
    stop("ERROR: Invalid method", call. = TRUE)
  }
  
  return(res)
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

#' Sliding window Fst calculation
#' 
#' Takes a tibble that has the a, b, and c parameters of
#' for Fst calculation (Weir & Cockerham 1984) and calculates
#' Fst over a sliding window. It only worts for a single
#' contig.
#'
#' @param dat A tibble with columns 'ref_id', 
#' 'ref_pos', and at leat one complete set of 
#' c('a', 'b','c') and/or c("n_c", "MSI", "MSP"). The first
#' set is used to calculate traditiona Fst according to
#' Weirr & Cockerham 1984 and the second calculates
#' Fst^{pool}.
#' @param w_size Window size in bp.
#' @param s_size Step size in bp.
#' @param sorted logical indicating if the sites in dat.
#' are already in sorted ascending order by ref_pos within each ref_id.
#'
#' @return A tibble with columns start, end, ref_id, nsites, and
#' Fst and/or Fst_pool columns depending on which estimators were calculated.
#' @export
#' 
#' @importFrom magrittr %>%
window_fst <- function(dat, w_size = 1000, s_size = 300, sorted = FALSE){
  # dat <- fst_pool$fst
  if(length(unique(dat$ref_id)) > 1){
    # If more than one contig/chromosome, recursively process each one
    Res <- dat %>%
      split(.$ref_id) %>%
      purrr::map_dfr(window_fst,
                     w_size = w_size,
                     s_size = s_size,
                     sorted = sorted) 
    # return(Res)
  }else{
    # If only one contig/chr
    if(!sorted){
      dat <- dat %>%
        dplyr::arrange(ref_pos)
    }
    
    Res <- NULL
    left_ii <- 1
    right_ii <- 2
    for(start in seq(from = 1,
                     to = max(dat$ref_pos) - s_size + 1,
                     by = s_size)){
      end <- start + w_size
      
      for(curr_left in left_ii:nrow(dat)){
        if(dat$ref_pos[curr_left] >= start)
          break
      }
      
      for(curr_right in right_ii:nrow(dat)){
        if(dat$ref_pos[curr_right] > end - 1)
          break
      }
      
      window <- dat[(curr_left):(curr_right - 1), ] %>%
        tidyr::drop_na()
      # dplyr::filter(!is.na(Fst))
      # window %>% print(n = w_size)
      
      res <- tibble::tibble(start = start, end = end,
                            ref_id = unique(window$ref_id),
                            n_sites = nrow(window))
      
      if(all(c("a", "b", "c") %in% colnames(window))){
        
        fst <- sum(window$a) / sum(window$a + window$b + window$c)
        res$Fst <- max(0, fst)
      }
      
      if(all(c("n_c", "MSI", "MSP") %in% colnames(window))){
        fst_pool <- sum(window$MSP - window$MSI) / sum(window$MSP + ((window$n_c - 1) * window$MSI))
        res$Fst_pool <- max(0, fst_pool)
      }
      
      Res <- Res %>%
        dplyr::bind_rows(res)
      left_ii <- curr_left
      right_ii <- curr_right - 1
    }
  }
  
  return(Res)
}
