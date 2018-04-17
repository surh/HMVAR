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

#' Check p-values in MKtest.py outfile
#' 
#' Opens a file from MKtest.py and uses the qvalue methods
#' of Storey and Tibshirani to estimate the proportion
#' of True Negatives (Pi0) in the data
#' 
#' @para file file name
#' @param which either 'all' or the name of the test to
#' be analyzed
#' @param plot whether to return a histogram of the p-values
#' 
#' @author Sur Herrera Paredes
#' 
#' @return  A named list where each element corresponds to
#' a test and each element is alist with pi0 and p1 entries.
#' See \link{check_pvalues} for more info
#' 
#' @export
check_pvals_in_file <- function(file, which, plot=TRUE){
  Tab <- read.table(file, header = TRUE, sep = "\t")
  default <- c("gene", "contig", "start", "end", "Dn", "Ds", "Pn", "Ps")
  
  if(nrow(Tab) == 0)
    return(list())
  
  # Get pval columns
  if(which == "all"){
    tests <- setdiff(colnames(Tab), default)
    tests <- tests[ grep(pattern = ".pval", x = tests, invert = TRUE) ]
    tests <- tests[ grep(pattern = ".perm", x = tests, invert = TRUE) ]
  }else{
    tests <- which
  }
  
  Res <- list()
  for(t in tests){
    # t <- tests[8]
    res <- check_pvalues(Tab[,t], Tab[,paste(t,".pval",sep="")], plot = plot)
    Res[[t]] <- res
  }
  
  return(Res)
}

#' Check p-values
#' 
#' Takes a pair of vectors of estimates and their p-values
#' for some hypothesis test, and uses qvalue methods from
#' storey and Tibshirani to estimate the number of True
#' Negatives (pi0) and creat a p-value histogram. It doesn't
#' consider cases where the estimate is undefined.
#' 
#' @param estimates numeric verctor of estimates
#' @param pvals numeric vector of p-values
#' @param plot plot or not
#' 
#' @return a list with elements pi0 with the point estimate
#' of pi0 and p1 which is a ggplot2 plot or NULL in case parameter
#' plot was FALSE
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
check_pvalues <- function(estimates, pvals, plot = TRUE){
  pvals <- pvals[ !is.na(estimates) ]
  qvals <- qvalue::qvalue(pvals)
  # qvals.sum <- summary(qvals)
  pi0 <- qvals$pi0
  
  p1 <- NULL
  if(plot){
    p1 <- ggplot(data.frame(pvals),aes(x = pvals)) +
      geom_histogram(bins = 20) +
      AMOR::theme_blackbox
  }
  
  res <- list(pi0 = pi0, p1 = p1)
  
  return(res)
}
