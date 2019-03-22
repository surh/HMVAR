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
#' 
#' @export
#' 
#' @references 
#' \url{http://GettingGeneticsDone.blogspot.com/}
#' \url{http://gettinggeneticsdone.blogspot.com/p/copyright.html}
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
