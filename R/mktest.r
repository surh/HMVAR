#' Get files for MK results
#' 
#' Takes a directory and finds all files that
#' correspond to results of MKtests.py, and returns
#' a list of full paths for those files
#' 
#' @param d Directory
#' @param pattern Pattern to match for files
#' containing results from MKtests.py
#' 
#' @return Vector of files
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
get_mk_results_files <- function(d, pattern = "^mk_results"){
  
  files <- list.files(d)
  chosen <- files[ grep(pattern = pattern, x = files) ]
  
  if(length(chosen) > 0)
    chosen <- paste(d, chosen, sep = "/")

  return(chosen)
}

#' Check p-values in MKtest.py outfile
#' 
#' Opens a file from MKtest.py and uses the qvalue methods
#' of Storey and Tibshirani to estimate the proportion
#' of True Negatives (Pi0) in the data
#' 
#' @param file file name
#' @param which either 'all' or the name of the test to
#' be analyzed
#' @param plot whether to return a histogram of the p-values
#' 
#' @author Sur Herrera Paredes
#' 
#' @return  A named list where each element corresponds to
#' a test and each element is alist with pi0 and p1 entries.
#' See \link{check_pvalues} for more info
#' 
#' @export
check_pvals_in_file <- function(file, which, plot=TRUE){
  Tab <- read.table(file, header = TRUE, sep = "\t")
  default <- c("gene", "contig", "start", "end", "Dn", "Ds", "Pn", "Ps")
  
  if(nrow(Tab) == 0)
    return(list())
  
  # Get pval columns
  if(which == "all"){
    tests <- setdiff(colnames(Tab), default)
    tests <- tests[ grep(pattern = ".pval", x = tests, invert = TRUE) ]
    tests <- tests[ grep(pattern = ".perm", x = tests, invert = TRUE) ]
  }else{
    tests <- which
  }
  
  Res <- list()
  for(t in tests){
    # t <- tests[8]
    res <- check_pvalues(Tab[,t], Tab[,paste(t,".pval",sep="")], plot = plot)
    Res[[t]] <- res
  }
  
  return(Res)
}

#' Check p-values
#' 
#' Takes a pair of vectors of estimates and their p-values
#' for some hypothesis test, and uses qvalue methods from
#' storey and Tibshirani to estimate the number of True
#' Negatives (pi0) and creat a p-value histogram. It doesn't
#' consider cases where the estimate is undefined.
#' 
#' @param estimates numeric verctor of estimates
#' @param pvals numeric vector of p-values
#' @param plot plot or not
#' 
#' @return a list with elements pi0 with the point estimate
#' of pi0 and p1 which is a ggplot2 plot or NULL in case parameter
#' plot was FALSE
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
check_pvalues <- function(estimates, pvals, plot = TRUE){
  pvals <- pvals[ !is.na(estimates) ]
  qvals <- qvalue::qvalue(pvals)
  # qvals.sum <- summary(qvals)
  pi0 <- qvals$pi0
  
  p1 <- NULL
  if(plot){
    p1 <- ggplot(data.frame(pvals),aes(x = pvals)) +
      geom_histogram(bins = 20) +
      AMOR::theme_blackbox
  }
  
  res <- list(pi0 = pi0, p1 = p1)
  
  return(res)
}


############# MIDAS MKTEST ###################
# Code for obtaining MKtest from midas merge output

#' Select samples that are present in mapping file
#'
#' @param abun A data table where the first column is called 'site_id', and all
#' the other columns correspond to sample names
#' @param map A data table where there is a column called 'sample' which
#' corresponds to the column names of 'abun'.
#'
#' @return A data table
#' @export
#'
#' @examples
#' @importFrom magrittr %>%
#' @importFrom dplyr select intersect
select_samples_from_abun <- function(abun, map){
  abun <- abun %>% dplyr::select(site_id, dplyr::intersect(map$sample, colnames(abun)) )
  
  return(abun)
}

#' Determine the effect of a coding variant on the aminoacid sequence
#' 
#' Takes an table corresponding to the contents of the snp_info.txt
#' file from midas_merge.py and adds a column indicating whether
#' the variant is synonymous or non-synonymous
#'
#' @param info A data table corresponding to the contents of
#' the snp_info.txt file produced by midas_merge. It must have columns:
#' 'site_id', 'major_allele', 'minor_allele' and 'amino_acids'. The aminoacid
#' column must contain a string of four comma-separated values indicating
#' the aminoacid encoded for each variant (eg. 'V,L,V,L').
#' @param nucleotides Named vector indicating the position of each
#' nucleotide that corresponds to the amino_acids column in info. The
#' default corresponds to MIDAS v1.3.1 default.
#'
#' @return The same data table passed as info with a factor column
#' 'snp_effect' added.
#' @export
#'
#' @examples
#' @importFrom dplyr select
#' @importFrom purrr pmap_chr
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#' @importFrom tibble add_column
find_snp_effect <- function(info, nucleotides=c(A = 1, C = 2, G = 3, T = 4)){
  snp_effect <- info %>%
    dplyr::select(site_id, major_allele, minor_allele, amino_acids) %>%
    purrr::pmap_chr(function(site_id, major_allele, minor_allele,
                             amino_acids, nucleotides){
      aa = stringr::str_split(string = amino_acids,
                              pattern = ",",
                              simplify = TRUE)
      if( aa[nucleotides[major_allele]] == aa[nucleotides[minor_allele]] ){
        type <- "synonymous"
      }else{
        type <- "non-synonymous"
      }
      return(type)},
      nucleotides = nucleotides)
  
  # Add column
  info <- info %>%
    tibble::add_column(snp_effect = factor(snp_effect,
                                           levels = c('synonymous',
                                                      'non-synonymous')))
  return(info)
}

determine_site_dist <- function(d, group_thres = 2){
  groups <- split(d$allele, d$Group)
  if(length(groups) == 1){
    dist <- NA
  }else{
    if(length(groups[[1]]) < group_thres || length(groups[[2]]) < group_thres){
      dist <- NA
    }else{
      g1 <- unique(groups[[1]])
      g2 <- unique(groups[[2]])
      
      if(length(g1) > 1 || length(g2) > 1){
        dist <- 'Polymorphic'
      }else if(g1 == g2){
        dist <- 'Invariant'
      }else if(g1 != g2){
        dist <- 'Fixed'
      }else{
        dist <- 'What!!!'
      }
    }
  }
  return(dist)
}

# determine_site_dist <- function(d){
#   # dat <- subset(dat, site_id == "703112")
#   
#   counts <- ftable(Group ~ allele, d)
#   # counts
#   
#   d1 <- diag(counts)
#   d2 <- counts[row(counts) - col(counts) != 0]
#   if(all(d1 > 0) && all(d2 == 0)){
#     type <- "fixed"
#   }else if(all(d1 == 0) && all(d2 > 1)){
#     type <- "fixed"
#   }else if(any(colSums(counts) == 0)){
#     type <- NA
# 
#   }else{
#     type <- "polymorphic"
#   }
#   # return(type)
#   return(tibble(site_id = d$site_id[1],
#                 type = type))
# }
# 
# system.time(determine_site_dist(d))
# T1 <- NULL
# T2 <- NULL
# for(i in 1:100){
#   t1 <- system.time(fun2(d))
#   t2 <- system.time(determine_site_dist(d))
#   
#   T1 <- rbind(T1, t1)
#   T2 <- rbind(T2, t2)
# }

calculate_snp_dist <- function(info, freq, depth, map, depth_thres = 1){
  # Reformat
  depth <- depth %>% gather(key = "sample", value = 'depth', -site_id)
  freq <- freq %>% gather(key = "sample", value = 'freq', -site_id)
  # meta <- info %>% select(site_id, ref_pos, snp_effect)
  
  dat <- depth %>%
    inner_join(freq, by = c("site_id", "sample")) %>%
    left_join(map, by = "sample") %>%
    filter(depth >= depth_thres) %>%
    mutate(allele = replace(freq, freq < 0.5, 'major')) %>%
    mutate(allele = replace(allele, freq >= 0.5, 'minor'))
  
  site_dist <- dat %>% split(.$site_id) %>% map_chr(determine_site_dist)
  site_dist <- tibble(site_id = names(site_dist), distribution = factor(site_dist,
                                                                        levels = c('Fixed', 'Invariant', 'Polymorphic')))
  info <- info %>% inner_join(site_dist, by = "site_id")
  
  return(info)
}

mkvalues <- function(info, depth_thres = 1){
  
  tab <- info %>% 
    select(snp_effect, distribution) %>%
    table(exclude = NULL, useNA = 'always')
  
  return(tibble(Dn = tab['non-synonymous', 'Fixed'],
                Ds = tab['synonymous', 'Fixed'],
                Pn = tab['non-synonymous', 'Polymorphic'],
                Ps = tab['synonymous', 'Polymorphic']))
}

midas_mktest <- function(midas_dir, map_file, genes, depth_thres = 1){
  # Read data
  map <- read_tsv(map_file)
  info <- read_tsv(paste0(midas_dir, "/snps_info.txt"),col_types = 'ccncccnnnnnccccc',
                   na = 'NA')
  depth <- read_midas_abun(paste0(midas_dir, "/snps_depth.txt"))
  freq <- read_midas_abun(paste0(midas_dir, "/snps_freq.txt"))
  
  # Process data
  # Rename map columns
  map <- map %>% select(sample = ID, Group) 
  # Clean info
  info <- info %>% select(-locus_type, -starts_with("count_"))
  # Clean depth and freq
  depth <- select_samples_from_abun(depth, map)
  freq <- select_samples_from_abun(freq, map)
  # Clean map
  map <- map %>% filter(sample %in% colnames(depth))
  
  # Select gene data
  info <- info %>% filter(gene_id %in% genes)
  freq <- freq %>% filter(site_id %in% info$site_id)
  depth <- depth %>% filter(site_id %in% info$site_id)
  
  # Calculate MK parameters
  # Calcualate snp effect
  info <- find_snp_effect(info)
  # Calculate snp dist
  info <- calculate_snp_dist(info = info,
                             freq = freq,
                             depth = depth,
                             map = map,
                             depth_thres = depth_thres)
  Res <- info %>%
    split(.$gene_id) %>%
    map_dfr(mkvalues,
            depth_thres = depth_thres,
            .id = "gene_id")
  
  return(Res)
}