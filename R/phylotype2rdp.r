#' Format HMMCP phylotype taxonomy for AMOR
#' 
#' @param x A character vector, one element per taxa
#' @param split.char A character string indicating the field
#' delimiter character
#' 
#' @return A character vector
#' 
#' @author Sur Herrera Paredes
phylotype2rdp <- function(x, split.char = ';'){
  for(i in 1:length(x)){
    x[i] <- gsub(pattern = "\\(\\d+\\)", replacement = "", x = x[i], perl = TRUE)
  }
  
  return(x)
}

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