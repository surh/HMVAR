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

#' Calculate genome-wide NI and alpha
#' 
#' Reads a file from MKtest.py and uses the TG
#' estimator to calculate a genome-wide neutrality
#' index as well as alpha
#' 
#' Uses estimatro NI_TG as fefined by Stoletzki & Eyre-Walker 2011.
#' Based on work from Tarone 1981 and Greenland 1982.
#' 
#' \deqn{NI_{TG} = \frac{\sum D_{si}P_{ni}/(P_{si} + D_{si})}{\sum P_{si}D_{ni}/(P_{si} + D_{si})}}
#' 
#' where i is the i-th gene.
#' 
#' @param file File to a tab-delimited file. Should have Dn,
#' Ds, Pn, and Ps columns, and one row per gene.
#' @param col_types Column definitions for \link{read_tsv}.
#' @param na Vector corresponding to NaN values.
#' 
#' @return vector with NI and alpha
#' 
#' @author Sur Herrera Paredes
#' 
#' @importFrom magrittr %>%
#' 
#' @export
genome_wide_ni <- function(file, na = c('nan',"NA","")){
  mkres <- readr::read_tsv(file,
                           col_types = readr::cols(.default = readr::col_double(),
                                                   gene_id = readr::col_character()),
                           na = na)

  # Cases where there were no genes that could be tested
  if(nrow(mkres) == 0)
    return(c(NA, NA))
  
  # Calculate genome-wide NI using the TG (Tarone 1981; Greenland 1982) estimator
  mkres <- mkres %>%
    dplyr::mutate(num = (Ds * Pn) / (Ps + Ds),
                  denom = (Ps * Dn) / (Ps + Ds))
  
  num <- sum(mkres$num, na.rm = TRUE)
  denom <- sum(mkres$denom, na.rm = TRUE)
  NI_TG <- num / denom
  alpha <- 1 - NI_TG
  
  return(c(NI_TG, alpha))
}
