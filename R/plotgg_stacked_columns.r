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


#' Plot columns as stacked bars
#' 
#' Takes a data frame or tibble with one row per observation and multiple
#' variables (columns) per observation. It then plots the selected variables
#' as stacked bars per observation, possibly faceting the plot.
#'
#' @param dat A data frame or tibble with one row per observation and multiple
#' variables (i.e. columns) per observations.
#' @param x The name of the column to be used as an observation name. It will
#' be used to as a label in the x axis of the plot.
#' @param columns A character vector containing the names of the columns to
#' stack. If NULL, all columns except for the column indicated by `x` and
#' the columns in `facet_formula` will be included.
#' @param facet_formula A formula to be used by \link[ggplot2]{facet_grid}. It
#' will split the stacked bars into panels.
#' @param gather_key The key name to be used by \link[tidyr]{gather}. It will indicate
#' the name for the gathered columns. The text itself will be used as the
#' legend title.
#' @param gather_value The value name to be used by \link[tidyr]{gather}. It will
#' indicate the sum of values in the stacked bars. The text itself will
#' be used as the y-axis title.
#' @param custom_theme The result of a call to \link[ggplot2]{theme}. It allows
#' for control over the plot format
#'
#' @return A ggplot2 plot.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#'
#' @examples
#' d <- tibble::tibble(sample = letters[1:4],
#'                     g1 = 4:1,
#'                     g2 = 1:4,
#'                     group = rep(LETTERS[1:2], 2))
#' plotgg_stacked_columns(d,
#'                        x = 'sample',
#'                        facet_formula = ~ group)
plotgg_stacked_columns <- function(dat, x, columns = NULL,
                                   facet_formula = NULL,
                                   gather_key = "key",
                                   gather_value = "value",
                                   custom_theme = ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                                 panel.grid = ggplot2::element_blank(),
                                                                 axis.text = ggplot2::element_text(color = "black"),
                                                                 axis.title = ggplot2::element_text(color = "black", face = "bold"),
                                                                 axis.line.x.bottom = ggplot2::element_line(),
                                                                 axis.line.y.left = ggplot2::element_line())){
  if(!(x %in% colnames(dat))){
    stop("ERROR: x must be a column of dat", call. = TRUE)
  }
  
  # Get column names if not specified
  if(is.null(columns)){
    columns <- setdiff(colnames(dat), c(all.names(facet_formula), x))
  }
  
  
  p1 <- dat %>%
    tidyr::gather(key = !!gather_key, value = !!gather_value, columns) %>%
    ggplot2::ggplot(ggplot2::aes_string(x = x, y = gather_value, fill = gather_key)) +
    ggplot2::geom_bar(stat = "identity") +
    custom_theme
  
  # Facet
  if(!is.null(facet_formula)){
    facet_formula <- formula(facet_formula)
    p1 <- p1 + ggplot2::facet_grid(facet_formula, scales = "free_x",
                                   space = "free_x")
  }
  
  return(p1)
}