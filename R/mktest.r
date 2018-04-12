#' Get files for MK results
#' 
#' Takes a directory and finds all files that
#' correspond to results of MKtests.py, and returns
#' a list of full paths for those files
#' 
#' @param d Directory
#' @param pattern Pattern to match for files
#' containing results from MKtests.py
#' 
#' @return Vector of files
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
get_mk_results_files <- function(d, pattern = "^mk_results"){
  
  files <- list.files(d)
  chosen <- files[ grep(pattern = pattern, x = files) ]
  
  if(length(chosen) > 0)
    chosen <- paste(d, chosen, sep = "/")

  return(chosen)
}


#' Calculate genome-wide NI and alpha
#' 
#' Reads a file from MKtest.py and uses the TG
#' estimator to calculate a genome-wide neutrality
#' index as well as alpha
#' 
#' @param file File
#' 
#' @return vector with NI and alpha
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
calculate_genome_wide_ni <- function(file){
  mkres <- read.table(file, sep = "\t", header = TRUE)
  
  # Cases where there were no genes that could be tested
  if(nrow(mkres) == 0)
    return(c(NA, NA))
  
  # Calculate genome-wide NI using the TG estimator
  num <- mkres$Ds * mkres$Pn / (mkres$Ps + mkres$Ds)
  denom <- mkres$Ps * mkres$Dn / (mkres$Ps + mkres$Ds)
  if(any(is.na(denom) != is.na(num))){
    stop("ERROR: There are not matching undefined values for NI_TG")
  }
  NI_TG <- sum(num, na.rm = TRUE) / sum(denom, na.rm = TRUE)
  alpha <- 1 - NI_TG
  
  return(c(NI_TG, alpha))
}
