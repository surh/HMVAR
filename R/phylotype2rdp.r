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

#' Format HMMCP phylotype taxonomy for AMOR
#' 
#' @param x A character vector, one element per taxa
#' @param split.char A character string indicating the field
#' delimiter character
#' 
#' @return A character vector
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
phylotype2rdp <- function(x, split.char = ';'){
  for(i in 1:length(x)){
    x[i] <- gsub(pattern = "\\(\\d+\\)", replacement = "", x = x[i], perl = TRUE)
  }
  
  return(x)
}

#' Phylotype to RDP version 2
#' 
#' Slow version
#' 
#' @param x A character vector, one element per taxa
#' @param split.char A character string indicating the field
#' delimiter character
phylotype2rdp2 <- function(x, split.char = ';'){
  sapply(strsplit(x, split = split.char), function(x){
    
    # Get last quality
    # last <- x[ length(x) ]
    # pat <- "\\((\\d+)\\)"
    # match <- regexpr(pattern = pat, text = last, perl = TRUE)
    # start.match <- attr(x = match, which = "capture.start")
    # length.match <- attr(x = match, which = "capture.length")
    # q <- substr(x = last, start = start.match, stop = start.match + length.match - 1)
    #return(q)
    
    # Clean
    x <- paste(gsub(pattern = pat, replacement = "", x = x), collapse = ";")
    return(x)
  })
  
}
