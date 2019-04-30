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

#' Varsites pipeline
#' 
#' Analyzes variable sites per sample and produces counts of
#' different types of variants per sample. It also analyzes variable
#' sites within sample to determine if they are homogeneous (fixed within sample)
#' or heterogeneous (variable within sample). It produces several plots
#'
#' @param freq A site x sample allele frequency table as a data frame or
#' tibble. It must have a 'site_id' column.
#' @param depth A site x sample sequence coverage table as a data frame
#' or tibble. It must have a 'site_id' column.
#' @param info A site x variable table. Based on the MIDAS snp_info.txt
#' files. It must have columns 'site_id', 'ref_id', 'ref_pos', 'amino_acids',
#' 'major_allele', and 'minor_allele'.
#' @param map A data frame or tibble with sample metadata. It must have a
#' 'sample' column that matches the sample names in `freq` and `depth`,
#' as well as a 'Group' column with a categorical variable to group the
#' samples.
#' @param depth_thres Minimum sequence coverage to to keep a site in
#' a given sample.
#' @param freq_thres Frequency trheshold for allele assignemnt per
#' sample. See \link{determine_snp_dist}.
#' @param plot If TRUE several plot objects will be included in the ouptut.
#'
#' @return A list with elements varsites and varsites.pos. Optionally,
#' ggplot2 objects are also included.
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
#' Res <- varsites_pipeline(freq = Dat$freq,
#'                          depth = Dat$depth,
#'                          info = Dat$info,
#'                          map = map)
varsites_pipeline <- function(freq, depth, info, map,
                              depth_thres = 1, freq_thres = 0.5, plot = TRUE){
  
  if(!("Group" %in% colnames(map))){
    stop("ERROR: map must have a Group column", call. = TRUE)
  }
  
  cat("Determining snp effect...\n")
  info <- determine_snp_effect(info)
  
  cat("Determining snp distribution...\n")
  info <- determine_snp_dist(info = info,
                             freq = freq,
                             depth = depth,
                             map = map,
                             depth_thres = depth_thres,
                             freq_thres = freq_thres,
                             clean = FALSE)
  
  # Subsitution type
  cat("Determining substitution type...\n")
  info <- determine_substitution_type(info, clean = FALSE)
  
  # Reformat data
  cat("Matching freq and depth...\n")
  dat <- match_freq_and_depth(freq = freq,
                              depth = depth,
                              info = info,
                              map = map,
                              depth_thres = depth_thres)
  
  # For each site determine if it is homogeneous or heterogenous
  cat("Determining within sample distribution...\n")
  dat <- determine_sample_dist(dat)
  
  # Count number of sites and number of variable per sample
  cat("Calcuating variable site types per sample...\n")
  varsites <- dat %>%
    dplyr::filter(depth >= 1) %>%
    split(.$sample) %>%
    purrr::map_dfr(function(d){
      nsites <- nrow(d)
      # sites <- d$freq
      nhomogeneous <- sum(d$sample_dist == "homogeneous", na.rm = TRUE)
      nheterogeneous <- sum(d$sample_dist == "heterogeneous", na.rm = TRUE)
      ntransitions <- sum(d$substitution == "transition", na.rm = TRUE)
      ntransvertions <- sum(d$substitution == "transversion", na.rm = TRUE)
      nsynonymous <- sum(d$snp_effect == "synonymous", na.rm = TRUE)
      nnonsynonymous <- sum(d$snp_effect == "non-synonymous", na.rm = TRUE)
      return(tibble::tibble(nsites = nsites,
                            n_homogeneous = nhomogeneous,
                            n_heterogeneous = nheterogeneous,
                            n_transitions = ntransitions,
                            n_transversions = ntransvertions,
                            n_synonymous = nsynonymous,
                            n_non.synonymous = nnonsynonymous))
    }, .id = "sample") %>%
    dplyr::inner_join(map, by = "sample")
  
  # Varsites per position
  cat("Calcuating variable site types per site...\n")
  varsites.pos <- variable_dist_per_site(dat = dat,
                                         variable = "sample_dist",
                                         group = "Group")
  
  # Prepare ouptut
  Res <- list(varsites = varsites, varsites.pos = varsites.pos,
              dnds.plot = NULL,
              variability.plot = NULL,
              pos.plot = NULL)
  
  if(plot){
    cat("Plotting...\n")
    # Plot synonymous & non-synonymous
    p1 <- plotgg_stacked_columns(dat = varsites,
                                 x = "sample",
                                 columns = c("n_synonymous", "n_non.synonymous"),
                                 facet_formula = ~ Group,
                                 gather_key = "type",
                                 gather_value = "nloci")
    Res$dnds.plot <- p1
    
    # Plot homogeneous & heterogeneous
    p1 <- plotgg_stacked_columns(dat = varsites,
                                 x = "sample",
                                 columns = c("n_homogeneous", "n_heterogeneous"),
                                 facet_formula = ~ Group,
                                 gather_key = "type",
                                 gather_value = "nloci")
    Res$variability.plot <- p1
    
    # Plot homogeneous & heterogeneous
    p1 <- ggplot2::ggplot(varsites.pos,
                          ggplot2::aes(x = ref_pos,
                                       y  = homogeneous / (homogeneous + heterogeneous))) +
      ggplot2::facet_grid(~ ref_id, space = "free_x", scales = "free_x") +
      ggplot2::geom_line(ggplot2::aes(color = Group)) +
      ggplot2::geom_point(ggplot2::aes(color = Group), size = 1, alpha = 0.2) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::theme(panel.background = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     axis.text = ggplot2::element_text(color = "black"),
                     axis.text.x = ggplot2::element_text(angle = 90),
                     axis.line.x.bottom = ggplot2::element_line(),
                     axis.line.y.left = ggplot2::element_line())
    Res$pos.plot <- p1
  }
  
  return(Res)
}
