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

#' Process map
#' 
#' Internal
#'
#' @param map Either a data frame or tibble or the path to
#' a file containing a table. Must have columns ID and Group.
#' @param map_file Same as map, for backwars compatibility
#' @param focal_group Level from Group column to use as
#' focal grou, i.e. compare everything else to this group.
#' 
#' @importFrom magrittr %>%
#' 
#' @return A tibble with columns sample and Group.
check_map <- function(map, map_file, focal_group = NA){
  # Backwards compatibility with map_file.
  if(missing(map) && !missing(map_file)){
    map <- map_file
  }
  
  # Process map
  if(is.character(map) && length(map) == 1){
    # Read and process map
    map <- readr::read_tsv(map,
                           col_types = readr::cols(.default = readr::col_character()))
    # Rename map columns
    map <- map %>% 
      dplyr::select(sample = ID, Group) 
  }else if(is.data.frame(map)){
    if(!all(c("sample", "Group")  %in% colnames(map))){
      stop("ERROR: map must contain sample and Group columns", call. = TRUE)
    }
  }else{
    stop("ERROR: map must be a file path or a data.frame/tibble")
  }
  
  # Compare everything to focal_group
  if(!is.na(focal_group)){
    map <- map %>%
      dplyr::mutate(Group = replace(Group,
                                    Group != focal_group,
                                    paste0("non.", focal_group)))
  }
  
  return(map)
}

#' Run system command
#' 
#' Internal
#' 
#' @param cmd System command
run_command <- function(cmd){
  cat("Running\n\t>", cmd, "\n")
  out <- system(cmd)
  
  return(out)
}

#' Select samples that are present in mapping file
#'
#' @param abun A data table where the first column is called 'site_id', and all
#' the other columns correspond to sample names
#' @param samples A vector with sample IDs.
#'
#' @return A data table
#'
#' @importFrom magrittr %>%
select_samples_from_abun <- function(abun, samples){
  abun <- abun %>% dplyr::select(site_id,
                                 dplyr::intersect(samples, colnames(abun)) )
  
  return(abun)
}

