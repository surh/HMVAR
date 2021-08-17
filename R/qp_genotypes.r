#' Quasi-phased genotypes
#' 
#' From allele frequency data in the format of midas_merge.py snps it determines
#' which samples and which positions in each sample are quasi-phaseable and
#' makes the corresponding genotype calls.
#'
#' @param midas_dir Directory with outptu from mdas_merge.py snps. Either
#' this or midas_dat needs to be provided.
#' @param midas_dat A list with info, freq, and depth elements corresponding
#' to the result of read_midas_data. Either this or midas_dir needs to be
#' provided
#' @param min_depth Minimum sequence depth at a given sample and position
#' to be abele to call a gnotype 
#' @param min_snv_prop Minimum proportion of SNVs in a sample, among SNVs that
#' pass the min_depth filter, that pass the maf_thres filter. Samples that
#' have a smaller proportion of SNVs passing this condition are considered
#' non-quasi-phaseable and genotype calls are not made.
#' @param maf_thres A value in the range [0, 0.5]. SNVs are only considered 
#' quasi-phaseable if their frequency (f) fullfills f <= maf_thres | f >= 1 -
#' maf_thres, and if they are in a quasi-phaseable sample. SNVs in
#' quasi-phaseable samples, that do not pass the maf_thres or min_depth
#' filters are considered missing data.
#' @param map A data table associating samples in the MIDAS otput to groups.
#' It must have an 'sample' and a 'Group' column. Only used if midas_dat was
#' passed and only samples in map are processed.
#'
#' @return A tiblle with columns site_id, and one column per quasi-phaseable
#' sample, which contain the genotype calls. 
#' @export
#' @importFrom magrittr %>%
qp_genotypes <- function(midas_dir, midas_dat,
                         map = NULL,
                         min_depth = 5,
                         min_snv_prop = 0.8,
                         maf_thres = 0.2){
  
  # Check parameters
  if(missing(midas_dir) & missing(midas_dat)){
    stop("ERROR: either midas_dir or midas_dat needs to be provided",
         call. = TRUE)
  }
  if(maf_thres < 0 || maf_thres > 0.5){
    stop("ERROR: maf_thres must be in the range [0, 0.5]", call. = TRUE)
  }
  if(min_snv_prop < 0 || min_snv_prop > 1){
    stop("ERROR: min_snv_prop must be in the range [0, 1]", call. = TRUE)
  }
  
  # Get or check midas_dat
  if(!missing(midas_dir)){
    midas_dat <- read_midas_data(midas_dir = midas_dir, map = map)
  }else if(!all(c("info", "freq", "depth") %in% names(midas_dat))){
    stop("ERROR: midas_dat must have info, freq, and depth elements",
         call. = TRUE)
  }
  
  
  Res <- HMVAR:::match_freq_and_depth(freq = midas_dat$freq,
                               depth = midas_dat$depth,
                               depth_thres = min_depth) %>%
    split(.$sample) %>%
    purrr::map_dfr(function(d, maf_thres = 0.2, min_snv_prop = 0.8){
      n_snvs <- nrow(d
      )
      # Identify SNVs in QP range
      d <- d[ d$freq <= maf_thres | d$freq >= 1 - maf_thres, ]
      
      # Keep only samples with enough QP snvs
      if(nrow(d) >= min_snv_prop * n_snvs){
        # assign allele. 1 for minor allele, 2 for major allele
        d$allele <- 1*(d$freq >= 1 - maf_thres) + 1
        d <- d %>%
          dplyr::select(site_id, allele)
      }else{
        d <- NULL
      }
      
      return(d)
    }, maf_thres = min(maf_thres, 1 - maf_thres),
    min_snv_prop = min_snv_prop, .id = "sample")
  
  if(nrow(Res) > 0){
    Res <- Res %>%
      dplyr::left_join(midas_dat$info %>%
                         dplyr::select(site_id, major_allele, minor_allele),
                       by = "site_id") %>%
      dplyr::mutate(allele = replace(allele, allele == 1, minor_allele[allele == 1])) %>%
      dplyr::mutate(allele = replace(allele, allele == "2", major_allele[allele == "2"])) %>%
      dplyr::select(-major_allele, -minor_allele) %>%
      tidyr::pivot_wider(names_from = "sample", values_from = "allele")
  }
    
  Res
}
