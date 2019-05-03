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

#' Calculate Direction of Selection
#' 
#' Calculates the Direction of Selection (DoS) statistic from
#' McDonald-Kreitman contingency tables.
#'
#' @param dat A gene by substitution type data frame or tibble.
#' Each row must correspond to a gene and it must have columns
#' Dn, Ds, Pn and Ps, corresponding to the McDonald-Kreitman
#' contingency table of each gene.
#' @param test It TRUE, a Z-test will be performed on the
#' DoS statistic.
#' @param clean If TRUE, genes where DoS is undefined (either
#' Dn and Ds are both zero, or Pn and Ps are both zero), will
#' be removed from the result.
#'
#' @return A tibble with the same columns as Dat plus a column
#' named DoS. If test is passed, a column named p.value will also
#' be included unless such a column already exists. In that case
#' a column with the name DoS.p.value will be included
#' 
#' @export
#' @importFrom magrittr %>%
calculate_dos <- function(dat, test = TRUE, clean = TRUE){
  dat <- dat %>%
    dplyr::mutate(DoS = (Dn / (Dn + Ds)) - (Pn / (Pn + Ps)))
  if(clean){
    dat <- dat %>%
      dplyr::filter(!is.na(DoS))
  }
  if(test){
    if('p.value' %in% colnames(dat)){
      dat <- dat %>%
        dplyr::mutate(DoS.p.value = 2*(1 - pnorm(q = abs(DoS) / sd(DoS, na.rm = TRUE))))
    }else{
      dat <- dat %>%
        dplyr::mutate(p.value = 2*(1 - pnorm(q = abs(DoS) / sd(DoS, na.rm = TRUE))))
    }
    
  }
  
  return(dat)
}

#' Calculate Direction of Selection from SNV data
#' 
#' Takes SNV data in MIDAS format and produces Direction
#' of Selection (DoS) statistics.
#' 
#' @param info Data table corresponding to the 'snps_info.txt'
#' file from MIDAS. Must have columns 'site_id' and 'sample'
#' @param freq A data table corresponding to the 'snps_freq.txt'
#' file from MIDAS. Must have a 'site_id' column, and one more
#' column per sample. Each row is the frequency of the minor
#' allele for the corresponding site in the corresponding sample.
#' @param depth A data table corresponding to the 'snps_depth.txt'
#' file from MIDAS. Must have a 'site_id' column, and one more
#' column per sample. Each row is the sequencing depth for the
#' corresponding site in the corresponding sample.
#' @param map A data table associating samples with groups (sites).
#' must have columns 'sample' and 'Group'.
#' @param depth_thres Minimum number of reads (depth) at a site at
#' a sample to be considered.
#' @param freq_thres Frequency cuttoff for minor vs major allele.
#' The value represents the distance from 0 or 1, for a site to be
#' assigned to the major or minor allele respectively. It must be
#' a value in [0,1].
#' @param test It TRUE, a Z-test will be performed on the
#' DoS statistic.
#' @param clean If TRUE, genes where DoS is undefined (either
#' Dn and Ds are both zero, or Pn and Ps are both zero), will
#' be removed from the result.
#'
#' @return A tibble with the same columns as Dat plus a column
#' named DoS. If test is passed, a column named p.value will also
#' be included unless such a column already exists. In that case
#' a column with the name DoS.p.value will be included
#' 
#' @export
#' @importFrom magrittr %>%
dos <- function(info, freq, depth, map,
                depth_thres = 1,
                freq_thres = 0.5,
                test = TRUE,
                clean = TRUE){
  
  Res <- calculate_mktable(info = info,
                           freq = freq,
                           depth = depth,
                           map = map,
                           depth_thres = depth_thres,
                           freq_thres = freq_thres) %>%
    calculate_dos(test = test, clean = clean)
  
  return(Res)
}

#' Direction of Selection statistics on MIDAS output
#' 
#' Take output from MIDAS script midas_merge.py
#' and obtain and test Direction of Selection (DoS)
#' statistics.
#'
#' @param midas_dir Directory where the output from
#' midas_merge.py is located. It must include files:
#' 'snps_info.txt', 'snps_depth.txt' and 'snps_freq.txt'.
#' @param map Either  a file path or a tibble. If a path,
#' it mos point to a mapping file associating samples
#' in the MIDAS otput to groups. It must have an 'ID'
#' and a 'Group' column. If a tibble. It must have columns
#' 'sample' and 'Group'.
#' @param map_file Same as map. Present for backwards compatibility.
#' @param genes The list of genes that are to be tested.
#' Must correspond to entries in the 'genes_id' column
#' of the 'snps_info.txt' file. If NULL, all genes will
#' be tested.
#' @param depth_thres The minimum number of reads at
#' a position in a given sample for that position in that
#' sample to be included.
#' @param freq_thres Frequency cuttoff for minor vs major allele.
#' The value represents the distance from 0 or 1, for a site to be
#' assigned to the major or minor allele respectively. It must be
#' a value in [0,1].
#' @param focal_group A string indicating the group that should
#' be compared to everything else. If different from NA, all values
#' in the 'Group' column of the mapping file will be converted to
#' "non.<focal_group>", ensuring a dichotomous grouping.
#'
#' @return A tibble with columns gene_id, Dn, Ds, Pn, Ps, DoS and p.value.
#' 
#' @export
#'
#' @examples
#' midas_dos(midas_dir = system.file("toy_example/merged.snps/", package = "HMVAR"),
#'           map = system.file("toy_example/map.txt", package = "HMVAR"))
midas_dos <- function(midas_dir, map,
                      map_file,
                      genes = NULL,
                      depth_thres = 1,
                      freq_thres = 0.5,
                      focal_group = NA){
  
  map <- check_map(map = map,
                   map_file = map_file,
                   focal_group = focal_group)
  
  # Read data
  Dat <- read_midas_data(midas_dir = midas_dir,
                         map = map,
                         genes = genes,
                         cds_only = TRUE)
  
  # DoS
  Res <- dos(info = Dat$info,
             freq = Dat$freq, 
             depth = Dat$depth,
             map = map,
             depth_thres = depth_thres,
             freq_thres = freq_thres,
             test = TRUE,
             clean = FALSE)
  
  return(Res)
}