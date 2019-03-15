library(HMVAR)
library(tidyverse)
library(mice)

#' Tidy mice
#' 
#' Internal utility function
#' 
#' Calls mice on data table
#'
#' @param d 
#' @param m 
#' @param verbose 
#'
#' @return A tibble with imputed results
#' 
#' @importFrom maggritr %>%
tidy_mice <- function(d, m = 5, verbose = FALSE){
  res <- d %>%
    dplyr::select(site_id, minor_allele, major_allele)
  
  d <- d %>%
    dplyr::select(-site_id, -minor_allele, -major_allele) %>%
    dplyr::filter_all(dplyr::any_vars(!is.na(.))) %>%
    mice::mice(m = 5, printFlag = verbose) %>% 
    mice::complete() %>%
    dplyr::as_tibble()
  
  res %>% dplyr::bind_cols(d)
}



mice_impute <- function(geno, snp,
                        outdir = "imputed/",
                        m = 5,
                        verbose = FALSE,
                        prefix = "imputed"){
  
  if(any(geno$site_id != snp$ID)){
    stop("ERROR: geno and snp tables do not match", call. = TRUE)
  }
  
  imp <- geno %>% split(snp$chr) %>%
    map_dfr(~tidy_mice(.), .id = "site_id", m1 = 1)
  imp
  
  
  # Re-organize files
  dir.create(outdir)
  
  filename <- paste0("output/", paste(c(prefix, "log.txt"), collapse = "."))
  file.copy(filename, outdir)
  file.remove(filename)
  
  filename <- paste0("output/", paste(c(prefix, "snpinfo.txt"), collapse = "."))
  file.copy(filename, outdir)
  file.remove(filename)
  
  filename <- paste0("output/", paste(c(prefix, "mean.genotype.txt"), collapse = "."))
  file.copy(filename, outdir)
  file.remove(filename)
  imputed_file <- file.path(outdir, paste(c(prefix, "mean.genotype.txt"), collapse = "."))
  
  return(imputed_file)
}




# data("airquality")
# airquality
# imp <- mice(airquality, m = 1)
# airquality_impute <- complete(imp)
# imp

map <- read_tsv("map.txt", col_types = cols(.default = col_character())) %>%
  select(sample = ID, Group = Group)
map
Dat <- midas_to_bimbam(midas_dir = "midas_output_small/", outdir = "bimbam", map = map, focal_group = "Supragingival.plaque", prefix = "test")


rm(imp)
imp <- Dat$Dat$geno %>% split(Dat$Dat$snp$chr) %>%
  map_dfr(~tidy_mice(.), m1 = 1)
imp




Dat$Dat$geno 

colnames(imp)[-1] == colnames(Dat$Dat$geno)[-(1:3)]
imp
Dat$Dat$geno %>% select(site_id, minor_allele)


table(is.na(Dat$Dat$geno[,-(1:3)]))
table(is.na(imp[,-1]))

summary(as.matrix(imp[,-1])[ is.na(Dat$Dat$geno[,-(1:3)]) & !is.na(imp[,-1]) ])
as.matrix(imp[,-1])[ is.na(Dat$Dat$geno[,-(1:3)]) & !is.na(imp[,-1]) ] %>% tibble(x = .) %>%
  ggplot(aes(x = x)) +
  geom_histogram(bins = 20) +
  AMOR::theme_blackbox()
