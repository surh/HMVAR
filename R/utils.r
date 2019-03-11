#' Run system command
#' 
#' Internal
#' 
#' @param cmd System command
run_command <- function(cmd){
  cat("Running\n\t>", cmd, "\n")
  out <- system(cmd)
  
  return(out)
}

