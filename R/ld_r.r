#' Calculate linkage disequilibrium (LD) R^2
#' 
#' Calculates r^2 = D ^2 / (pA(1-PA)pB(1-pB))
#'
#' @param x A SNP numerical matrix with SNPs as rows & samples
#' as columns 
#'
#' @return A distance object
#' @export
#' @importFrom magrittr %>%
#' 
ld_r <- function(x) {
  apply(x,1,function(A,B){
    
    pA <- mean(A, na.rm = TRUE)
    pB <- rowMeans(B, na.rm = TRUE)
    pAB <- colMeans(t(B) * A, na.rm = TRUE)
    
    D <- pAB - (pA * pB)
    
    # Calculate r.squared
    return((D^2)/(pA*(1-pA) * pB * (1-pB)))
  }, B = x) %>% as.dist()
  
}



