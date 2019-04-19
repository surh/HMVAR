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
#' # Get file paths
#' midas_dir <- system.file("toy_example/merged.snps/", package = "HMVAR")
#' map <- readr::read_tsv(system.file("toy_example/map.txt", package = "HMVAR"),
#'                        col_types = readr::cols(.default = readr::col_character())) %>%
#'   dplyr::select(sample = ID, Group)
#' 
#' # Read data
#' midas_data <- read_midas_data(midas_dir = midas_dir, map = map, cds_only = TRUE)
#' 
#' determine_substitution_type(midas_data$info, clean = FALSE)
determine_substitution_type <- function(info, clean = TRUE){
  info <- info %>%
    add_column(substitution = info %>%
                 pmap_chr(function(major_allele, minor_allele, ...){
                   base_type <- c(A = "purine", C = "pyrimidine", G = "purine", T = "pyrimidine")
                   if(base_type[major_allele] == base_type[minor_allele]){
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
