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

#' Read eggNOG mapper annotations
#'
#' @param path File path
#' @param format 
#'
#' @return A tibble
#' 
#' @export
read_eggnog <- function(path, format = "emapper1"){
  if(format == "emapper1"){
    res <- readr::read_tsv(path,
                           comment = "#",
                           col_names = c("query_name", "seed_eggNOG_ortholog", "seed_ortholog_evalue",
                                         "seed_ortholog_score",	"predicted_gene_name", "GO_terms",
                                         "KEGG_KOs", "BiGG_reactions", "Annotation_tax_scope", "OGs",
                                         "bestOG|evalue|score", "COG_cat", "eggNOG_annot"),
                           col_types = readr::cols(.default = readr::col_character(),
                                                   seed_ortholog_score = readr::col_double(),
                                                   seed_ortholog_evalue = readr::col_double()))
  }else if(format == "uhgg"){
    # Read weird UHGG eggnog annots
    res <- read_lines(path, n_max = Inf, skip = 1) %>%
      str_split_fixed(pattern = "\t", n = 22) %>%
      as_tibble(.name_repair = function(x) paste0("X", 1:length(x)))
    res[ res == "" ] <- NA
    colnames(res) <- c("query_name", "seed_eggNOG_ortholog", "seed_ortholog_evalue",
                         "seed_ortholog_score", "best_tax_level",	"predicted_gene_name",
                         "GO_terms", "EC_terms", "KEGG_KOs", "KEGG_pathways","KEGG_modules",
                         "KEGG_reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
                         "BiGG_reactions", "tax_scope", "OGs",
                         "bestOG|evalue|score", "COG_cat", "eggNOG_annot")
  }else{
    stop("ERROR: format must be 'emapper1' or 'uhgg'", call. = TRUE)
  }
  

      
  return(res)
}
