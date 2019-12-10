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

#' Calculate taxon prevalence
#' 
#' Calculates the prevalence of each
#' taxon overall or by some grouping factor.
#' 
#' @param Dat a dataset object
#' @param thres Minimum number of reads for a taxon
#' in a sample to be counted as present.
#' @param group A grouping variable
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
calculate_prevalence <- function(Dat, thres = 1, group = NULL){
  if(class(Dat) != "Dataset")
    stop("ERROR: A Dataset object must be passed", call. = TRUE )
  
  if(is.null(group)){
    group <- rep("All", length.out = ncol(Dat.bin$Tab))
    group_n <- table(group)
    varnames <- c("Group", "Taxon")
  }else{
    group <- Dat$Map[ , group ]
    group_n <- table(group)
    varnames <- c("Taxon", "Group")
  }
  
  Dat.bin <- create_dataset(Tab = 1*(Dat$Tab >= thres), Map = Dat$Map, Tax = Dat$Tax)
  Dat.bin <- pool_samples.default(Tab = Dat.bin$Tab, groups = group, FUN = sum)
  
  # res <- data.frame(Taxon = row.names(Dat$Tab), Count = rowSums(Dat$Tab >= thres))
  # res$Proportion <- res$Count / ncol(Dat$Tab)
  
  Res <- reshape2::melt(Dat.bin$Tab, varnames = varnames, value.name = "Count")
  Res <- Res[ , c("Taxon", "Group", "Count") ]
  
  if(length(group_n) == 1)
    Res$Group <- "All"
  
  Res$Proportion <- as.vector(Res$Count / group_n[ as.character(Res$Group) ])
  
  return(Res)
}
