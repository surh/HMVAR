#' p-value qqplot
#' 
#' quantile-quantile plot of expected vs observed
#' p-values. With base graphics.
#' 
#' Converts zeroes to the smalles non-zero p-value which
#' is required for the log-transformation. It uses base
#' graphics becasue ggplot2 can be very inefficient when
#' there are millions of p-values.
#' 
#' @param pvector Vector of p-values
#' @param main Title
#' @param ... Extra parameters for plot
#'
#' @return A base plot
#' @export
#' 
#' @references 
#' \url{http://GettingGeneticsDone.blogspot.com/}
#' \url{http://gettinggeneticsdone.blogspot.com/p/copyright.html}
#' 
#' @examples
pval_qqplot <- function(pvector, main=NULL, ...) {
  
  # removing zeroes before log-transformation
  pvector[ pvector == 0 ] <- min(pvector[ pvector > 0 ])
  
  o <- -log10(sort(pvector,decreasing=F))
  e <- -log10( 1:length(o)/length(o) )
  plot(e, o, pch = 19, cex = 1, main = main, ...,
       xlab = expression(Expected~~-log[10](italic(p))),
       ylab = expression(Observed~~-log[10](italic(p))),
       xlim = c(0, max(e)), ylim=c(0, max(o)))
  lines(e, e, col="red")
}