#' Frequency Increment Test (FIT)
#' 
#' From https://www.genetics.org/content/196/2/509
#'
#' @param Dat One of three options:
#' 1. A list from \link{read_midas_data}
#' 2. A tibble corresponding to the freq element of 1.
#' 3. A tibble with all the in long format must have site_id, freq, pop_id
#' & time columns
#' @param map For 1 and 2 a map indicating sample population and timepoint
#' @param min_depth For 1, min depth to consider a site.
#' @param aggregate 
#'
#' @return A tibble. Empty if no site x pop_id combination has at least three
#' timepoints
#' @export
#' @importFrom magrittr %>%
FIT <- function(Dat, aggregate = TRUE, ...){ UseMethod("FIT") }

#' @export
FIT.data.frame <- function(Dat, aggregate = TRUE, map, min_depth = 1){
  
  expanded_cols <- c("site_id", "freq", "pop_id", "time")
  if( all(expanded_cols %in% colnames(Dat)) ){
    Res <- Dat %>%
      dplyr::filter(freq != 1 & freq != 0) %>%
      split(list(.$site_id, .$pop_id)) %>%
      purrr::map_dfr(function(d){
        
        if(nrow(d) < 3)
          return(NULL)
        
        pop_id <- unique(d$pop_id)
        site_id <- unique(d$site_id)
        # ref_id <- unique(d$ref_id)
        # ref_pos <- unique(d$ref_pos)
        
        d %>%
          dplyr::arrange(time) %>%
          dplyr::mutate(vprev = dplyr::lag(freq, n = 1, default = NA),
                        t = as.numeric(time - min(time))) %>%
          dplyr::mutate(tprev = dplyr::lag(t, n = 1, default = NA)) %>%
          dplyr::mutate(Y = (freq - vprev) / sqrt(2 * vprev * (1 - vprev) * (t - tprev))) %>%
          dplyr::summarise(Y.mean = mean(Y, na.rm = TRUE),
                           Y.s2 = var(Y, na.rm = TRUE),
                           L = length(Y) - 1) %>%
          dplyr::mutate(t.fi = Y.mean / sqrt(Y.s2 / L)) %>%
          dplyr::mutate(pval = 2 * (1 - pt(q = abs(t.fi), df = L - 1)),
                        pop_id = pop_id,
                        site_id = site_id)
        
      })
    
    if(aggregate){
      Res <- Res %>%
        dplyr::filter(!is.na(t.fi)) %>%
        dplyr::filter(!is.infinite(t.fi)) %>%
        dplyr::group_by(site_id) %>%
        dplyr::summarise(y.mean = sum(Y.mean/Y.s2) / sum(1/Y.s2),
                         y.var = 1 / sum(1/Y.s2),
                         n_pops = length(Y.mean)) %>%
        dplyr::mutate(z.score = y.mean / sqrt(y.var)) %>%
        dplyr::mutate(pval = 2 * pnorm(abs(z.score), mean = 0, sd = 1, lower.tail = FALSE)) %>%
        dplyr::arrange(desc(abs(z.score)))
    }
  }else{
    stop("ERROR: only option 3 for Dat is available now")
  }

  return(Res)
}
