#' Formats input form HMMCP 16S files
#' 
#' Only to be used to process 16S files from MOTHUR \
#' pipeline of HMP (HMMCP). Count files must be pre-edited
#' to remove trailing tab
#' 
#' @return A Dataset object
#' 
#' @author Sur Herrera Paredes
#' 
format_input <- function(name, counts_file, taxonomy_file, Map = Map,
                         collapse_level = NULL){
  # name <- files$Name[1]
  # counts_file <- files$counts[1]
  # taxonomy_file <- files$taxonomy[1]
  
  # Count table
  Tab <- read.table(file = counts_file,
                    sep = "\t", row.names = 1, header = TRUE)
  Tab <- t(Tab)
  # dim(Tab)
  # Tab[1:5,1:5]
  
  # Taxonomy
  Tax <- read.table(file = taxonomy_file,
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  colnames(Tax) <- c("ID", "Taxonomy")
  Tax$Taxonomy <- phylotype2rdp(Tax$Taxonomy)
  Tax$ID <- paste("X", Tax$ID, sep = "")
  row.names(Tax) <- Tax$ID
  # head(Tax)
  
  # Create dataset
  to_remove <- c("positive_control.PPS", "positive_control.may1",
                 "positive_gd.PPS", "positive_mock.PPS", 
                 "water_blank.PPS", "water_blank.may1")
  Tab <- Tab[ , setdiff(colnames(Tab), to_remove) ]
  
  if(length(setdiff(colnames(Tab), row.names(Map))) > 0)
    stop("ERROR1")
  if(length( setdiff(row.names(Tab), row.names(Tax))) > 0)
    stop("ERROR2")
  
  Map <- Map[ colnames(Tab), ]
  Tax <- Tax[ row.names(Tab), ]
  Dat <- create_dataset(Tab = Tab, Map = Map, Tax = Tax)
  
  if(!is.null(collapse_level)){
    Dat <- collapse_by_taxonomy(Dat = Dat, Group = NULL,
                                level = collapse_level, sepchar = ";", FUN = sum)
  }
  
  return(Dat)
}