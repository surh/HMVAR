#!/usr/bin/env Rscript
library(HMVAR)
library(argparser)
library(tidyverse)

process_arguments <- function(){
  p <- argparser::arg_parser(paste0("Calculates genome-wide neutrality index from ",
                                    "Anchurus MKtest.py results"))
  
  # Positional arguments
  p <- argparser::add_argument(p, "input", help = paste0("Results from MKtest.py. ",
                                                         "must have Dn, Ds, Pn, and Ps columns, ",
                                                         "as well as one gene per line. Must ",
                                                         "correspond to a single genome file, or a ",
                                                         "directory of one file per genome."),
                               type = "character")
  
  # Optional arguments
  p <- argparser::add_argument(p, "--outfile", help = "File with results",
                               default = "genomewide_ni.txt", type = "character")
  p <- argparser::add_argument(p, "--type", help = "Type of input (file or dir).",
                               default = "file", type = "character")
  
  # Read arguments
  args <- argparser::parse_args(p)
  
  if(!(args$type %in% c("file", "dir"))){
    stop("ERROR: type must be either file or dir.", call. = TRUE)
  }
  
  return(args)
}

# Get arguments
args <- process_arguments()
print(args)

if(args$type == 'file'){
  Res <- genome_wide_ni(args$input)
  Res <- tibble(File = basename(args$input),
                NI = Res[1],
                alpha = Res[2])
}else if (args$type == "dir"){
  file_list <- list.files(args$input)
  Res <- NULL
  for(f in file_list){
    species <- basename(f)
    cat(species, "\n")
    
    res <- genome_wide_ni(paste0(args$input, "/", f))
    res <- tibble(File = species,
                  NI = res[1],
                  alpha = res[2])
    Res <- rbind(Res, res)
  }
}

# Write output
Res %>% arrange(desc(NI))
write_tsv(Res, args$outfile)
