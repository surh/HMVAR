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

#' Manhattan plot with ggplot2
#' 
#' Creates manhattan plot
#'
#' @param dat A data.frame or tibble with columns 'ref_id', 'ref_pos'
#' and one column matching `pval_column`.
#' @param refs A data.frame or tibble with columns 'ref_id',
#' 'start' and 'end'.
#' @param pval_thres p-value threshold value to draw a line.
#' @param pval_column Name of the column with the p-values.
#' @param log_transform If TRUE, values in `pval_column` will be
#' transformed via -log10(p-value).
#' @param color_refs If TRUE, points from different refs will have different
#' colors.
#' @param alpha Transparency parameters for the points
#' @param custom_theme A \link[ggplot2]{theme} function or call to format the
#' plot.
#'
#' @return A ggplot2 plot
#' 
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' # Create some fake data
#' set.seed(12345)
#' dat <- tibble::tibble(ref_id = rep(LETTERS[1:5], each = 20),
#'                       ref_pos = sample(1:200,size = 100, replace = FALSE),
#'                       p.value = runif(100)^4)
#' 
#' # Create some fake ref sizes
#' refs <- tibble::tibble(ref_id = LETTERS[1:5],
#'                        start = 0,
#'                        end = c(600,500,450,400,200))
#' 
#' # Make some plots
#' plotgg_manhattan(dat)
#' plotgg_manhattan(dat, refs = refs)
#' plotgg_manhattan(dat, refs = refs, color_refs = FALSE)
plotgg_manhattan <- function(dat, refs = NULL,
                             pval_thres = 1e-6,
                             pval_column = "p.value",
                             log_transform = TRUE,
                             color_refs = TRUE,
                             alpha = 0.2,
                             custom_theme = ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                           panel.grid = ggplot2::element_blank(),
                                                           axis.text = ggplot2::element_text(color = "black"),
                                                           axis.title = ggplot2::element_text(color = "black", face = "bold"),
                                                           axis.line.x.bottom = ggplot2::element_line(),
                                                           axis.line.y.left = ggplot2::element_line())){
  
  # Check columns in dat
  if(!all(c("ref_id", "ref_pos", pval_column) %in% colnames(dat))){
    stop("ERROR: missing columns in dat", call. = TRUE)
  }
  
  # Match sites with start and end of its ref
  if(!is.null(refs)){
    if(!all(c("ref_id", "start", "end") %in% colnames(refs))){
      stop("ERROR: missing columns in refs", call. = TRUE)
    }
    dat <- dat %>% dplyr::full_join(refs, by = "ref_id")
  }
  
  # Create base plot
  p1 <- ggplot2::ggplot(dat,ggplot2::aes_string(x = "ref_pos", y = pval_column)) +
    ggplot2::facet_grid(~ ref_id, scales = "free_x", space = "free_x") +
    ggplot2::geom_hline(yintercept = pval_thres, color = "red", size = 2)
  
  # Add points with or without color per ref
  if(color_refs){
    p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = ref_id), alpha = alpha)
  }else{
    p1 <- p1 + ggplot2::geom_point(alpha = alpha)
  }
  
  # Set axis of each ref to its defined start and end
  if(!is.null(refs)){
    p1 <- p1 +
      ggplot2::geom_blank(ggplot2::aes(x = start)) +
      ggplot2::geom_blank(ggplot2::aes(x = end))
  }
  
  # Log10-transform the y-axis
  if(log_transform){
    minus.log10 <- scales::trans_new("minus.log10",
                                     transform = function(x) -log10(x),
                                     inverse = function(x) 10^(-x),
                                     breaks = scales::log_breaks(base = 10),
                                     domain = c(0,Inf))
    p1 <- p1 +
      ggplot2::scale_y_continuous(trans = minus.log10)
    
  }
  
  # Format
  p1 <- p1 +
    ggplot2::guides(color = FALSE) +
    custom_theme
  
  return(p1)
}