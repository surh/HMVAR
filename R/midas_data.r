# Generic midas management functions

#' Read MIDAS abundance file
#' 
#' Reads either the snps_depth.txt or snps_freq.txt
#' file produced by midas_merge.py
#'
#' @param file File path
#'
#' @return A tibble
#' @export
#'
#' @examples
#' @importFrom readr read_tsv
read_midas_abun <- function(file){
  abun <- readr::read_tsv(file,
                          na = 'NA', n_max = 10)
  abun <- readr::read_tsv(file,
                          na = 'NA',
                          col_types = paste(c('c',rep('n', ncol(abun) - 1)),
                                            collapse = ''))
  
  return(abun)
}