#' Calculate selection coefficient
#' 
#' Selection coefficient for alleles that change frequency between two
#' time points.
#' 
#' According to https://www.genetics.org/content/196/2/509
#'
#' @param Dat There are 3 options:
#' 
#' 1. A list like the one produced by \link{HMVAR::read_midas_data}. It
#' must have freq and depth elements.
#' 2. A data frame or tibble corresponding to the freq element of option 1.
#' 3. A data frame or tibble already including the initila (v0) and final (v1)
#' frequencies of every allele. It must have columns "site_id", "v0",and "v1".
#' For option 1 and 2, and additonal map is required with time information
#' for the samples. For option 3 an additional time parameter is needed.
#' @param aggregate Whether to calculate aggregated statistics from
#' multiple populations
#' @param time Only used (and required) for option 3 in Dat. Either a single
#' numeric value indicating the time between the first and last allele
#' frequency measurement, or a single character value indicating the name of
#' the column containing the time between first and last allele frequency
#' measurement.
#' 
#' @param map 
#' 
#' @return A tibble with selection coefficient calculations
#' @export
#' @importFrom magrittr %>%
s_coefficient <- function(Dat, ...){
  UseMethod("s_coefficient")
}

#' @export
s_coefficient.data.frame <- function(Dat, aggregate = TRUE, time, map){

  maf_cols <- c("site_id", "v0", "v1", "pop_id")
  
  if( all(maf_cols %in% colnames(Dat)) ){
    
    # Check that Dat and time are correctly passed
    if(missing(time)){
      stop("ERROR: Time must be passed when pasing a table of maf changes", 
           call. = TRUE)
    }
    if(length(time) != 1){
      stop("ERROR: time must be one value exactly", call. = TRUE)
    }
    
    # Extract needed data, add time and fix names
    if(class(time) == "character"){
      time_col <- time
      maf_cols <- c(maf_cols, time_col)
      Dat <- Dat %>%
        dplyr::select(dplyr::all_of(maf_cols)) %>%
        dplyr::rename(time = dplyr::all_of(time_col))
      
    }else if(class(time) == "numeric"){
      time_col <- "time"
      Dat <- Dat %>%
        dplyr::select(dplyr::all_of(maf_cols)) %>%
        dplyr::mutate(!!time_col := time)
      
    }
    
    # Calculate selection coefficient
    Res <- Dat %>%
      dplyr::filter(!(v0 == 0 & v1 == 0)) %>%
      dplyr::mutate(s = (1 / time) * log( (v1 / (1 - v1)) * ((1 - v0) / v0) ))
    
    # Aggregate multiple populations
    if(aggregate){
      Res <- Res %>%
        dplyr::filter(!is.na(s)) %>%
        dplyr::filter(!is.infinite(s)) %>%
        dplyr::group_by(site_id) %>%
        dplyr::summarise(s.mean = mean(s),
                         s.sd = sd(s),
                         n_pops = length(s),
                         .groups = 'drop') %>%
        dplyr::mutate(t.s = s.mean / (s.sd / sqrt(n_pops))) %>%
        dplyr::mutate(pval = 2 * pt(abs(t.s), df = n_pops - 1, lower.tail = FALSE)) %>%
        dplyr::arrange(pval) 
    }
  }else{
    stop("ERROR: only MAF changes method implemented")
  }
  
  return(Res)
}
  