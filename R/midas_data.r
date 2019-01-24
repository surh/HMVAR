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

#' Read midas snp merge data
#' 
#' Reads the output of midas_merge.py
#'
#' @param midas_dir Directory where the output from
#' midas_merge.py is located. It must include files:
#' 'snps_info.txt', 'snps_depth.txt' and 'snps_freq.txt'.
#' @param map A data table associating samples
#' in the MIDAS otput to groups. It must have an 'sample'
#' and a 'Group' column.
#' @param genes The list of genes that are to be tested.
#' Must correspond to entries in the 'genes_id' column
#' of the 'snps_info.txt' file. If NULL, all genes will
#' be tested.
#'
#' @return A list with three data tables: info, freq
#' and depth, corresponding to the three files from 
#' MIDAS
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @importFrom dplyr select filter
read_midas_data <- function(midas_dir, map, genes){
  # Read data
  info <- readr::read_tsv(paste0(midas_dir, "/snps_info.txt"),
                          col_types = 'ccncccnnnnnccccc',
                          na = 'NA')
  depth <- read_midas_abun(paste0(midas_dir, "/snps_depth.txt"))
  freq <- read_midas_abun(paste0(midas_dir, "/snps_freq.txt"))
  
  # Process data
  # Clean info
  info <- info %>% 
    dplyr::select(-locus_type, -tidyselect::starts_with("count_"))
  # Clean depth and freq
  depth <- select_samples_from_abun(depth, map)
  freq <- select_samples_from_abun(freq, map)
  # Clean map
  map <- map %>% 
    dplyr::filter(sample %in% colnames(depth))
  
  # Select gene data
  if(is.null(genes)){
    info <- info %>% 
      dplyr::filter(!is.na(gene_id))
  }else{
    info <- info %>%
      dplyr::filter(gene_id %in% genes)
  }
  freq <- freq %>% 
    dplyr::filter(site_id %in% info$site_id)
  depth <- depth %>% 
    dplyr::filter(site_id %in% info$site_id)
  
  return(list(info = info, freq = freq, depth = depth))
}