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

expand_annot <- function(gene_id, terms){
  terms <- stringr::str_split(string = terms, pattern = ",") %>%
    unlist
  tibble::tibble(gene_id = gene_id,
                 term = terms)
}

annots_to_geneGO <- function(annots, direction = "geneID2GO"){
  if(!(c("gene_id", "terms") %in% colnames(annots))){
    stop("ERROR: missing columns in annots")
  }
  
  if(direction == "geneID2G"){
    annots %>%
      purrr::pmap_dfr(expand_annot) %>%
      dplyr::filter(!is.na(term)) %>%
      split(.$gene_id) %>%
      purrr::map(~ .x$term)
    
  }else if(direction == "GO2geneID"){
    annots <- annots %>%
      purrr::pmap_dfr(expand_annot) %>%
      dplyr::filter(!is.na(term)) %>%
      split(.$term) %>%
      purrr::map(~ .x$gene_id)
  }else{
    stop("ERROR: direction must be geneID2G or GO2geneID", call. = TRUE)
  }
  
  return(annots)
}


library(topGO)
library(ALL)
data(ALL)
data(geneList)

affyLib <- paste(annotation(ALL), "db", sep = ".")
affyLib
library(package = affyLib, character.only = TRUE)

sampleGOdata <- new("topGOdata",
                    description = "Simple session",
                    ontology = "BP",
                    allGenes = geneList,
                    geneSelectionFun = topDiffGenes,
                    nodeSize = 10,
                    annotationFun = annFUN.db,
                    affyLib = affyLib)



dat <- read_eggnog("~/micropopgen/exp/2019/2019-04-01.hmp_subsite_annotations/annotations/Catonella_morbi_61904.emapper.annotations")

gene_column <- "query_name"
annot_colum <- "GO_terms"

annots <- dat %>%
  dplyr::select(gene_id = gene_column, terms = annot_colum)
annots









