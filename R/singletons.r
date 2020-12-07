# (C) Copyright 2020 Sur Herrera Paredes
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

#' Find singletons
#' 
#' Finds variants that have singleton alleles and returns the
#' positions and the samples where they ocurr.
#'
#' @param alleles site x sample alleles table. Must have column site_id plus
#' one column per sample
#' @param midas_dir Instead of alleles table. A path with output files
#' in the format of midas merge can be provided. Must have freqs, depth and
#' info files.
#' @param info Optional site info table to append metadata to singleton
#' sites. Must have a site_id column.
#' @param meta Optional sample metadata to append to samplew with singleton
#' sites. Must have column sample.
#'
#' @return a tibble with one row per singleton
#' @export
#' 
#' @importFrom magrittr %>%
find_singletons <- function(alleles, midas_dir, info,
                            meta){
  if(missing(alleles) & missing(midas_dir)){
    stop("ERROR: either alleles or midas_dir needs to be provided",
         call. = TRUE)
  }else if(missing(alleles) & !missing(midas_dir)){
    stop("ERROR: finding singletons from midas_dir no implemented",
         call. = TRUE)
    # Here create alleles and info objects
  }
  
  cat("\tFinding singletons...\n")
  Sing <- alleles %>%
    tidyr::pivot_longer(-site_id,
                        names_to = "sample",
                        values_to = "allele") %>%
    split(.$site_id) %>%
    purrr::map_dfr(function(d){
      tab <- table(d$allele, useNA = "no")
      
      # Fail if site not bi-allelic
      if(length(tab) != 2){
        stop("ERROR")
      }
      
      if( sum(tab == 1) == 1 ){
        allele <-  names(tab)[which(tab == 1)]
        sample_id <- d$sample[d$allele == allele]
        sample_id <- sample_id[!is.na(sample_id)]
        
        # Return singleton
        return(tibble::tibble(n_samples = sum(!is.na(d$allele)),
                              allele = allele,
                              sample = sample_id))
        
      }else{
        # ignore non-singletons and cases of two samples with a different
        # allele each
        return(NULL)
      }
    }, .id = "site_id")
  
  # combine site and sample metadata
  cat("\tMerging metadata...\n")
  if(!missing(site)){
    Sing <- Sing %>%
      left_join(info,
                by = 'site_id')
  }
  if(!missing(meta)){
    Sing <- Sing  %>%
      left_join(meta, by = "sample")
  }
  
  return(Sing)
}


#' Test for singleton enrichments
#' 
#' Takes table of singletons and overall SNVs and tests for an enrichment
#' of overall singletons of enrichments per gene, as well as an erichment
#' of non-synonymous over synonymous mutations among singletons versus
#' among non-singletons
#'
#' @param Dat A table with one singelton site per row. Must have columns
#' site_id, gene_id and snp_effect.
#' @param info A table with information about all SNV sites (singleton and
#' non-singleton) in the genes
#' being tested. Must have columns site_id, gene_id, and snp_effect.
#'
#' @return A tibble
#' @export
#' 
#' @importFrom magrittr %>%
test_singleton_enrichment <- function(Dat, info){
  if(missing(Dat) || missng(info)){
    stop("ERROR: both Dat and info must be provided")
  }
  
  Dat %>%
    # Calculate # singletons by effect
    dplyr::group_by (gene_id) %>%
    dplyr::summarise(n_genesing = length(site_id),
                     n_genesing_ns  = sum(snp_effect == "non-synonymous"),
                     n_genesing_s = sum(snp_effect == "synonymous"),
                     .groups = "drop") %>%
    # Join with # snvs per gene by effect
    dplyr::left_join(info %>%
                       dplyr::group_by(gene_id) %>%
                       dplyr::summarise(n_genesnvs = length(site_id),
                                        n_genesnvs_ns = sum(snp_effect == "non-synonymous"),
                                        n_genesnvs_s = sum(snp_effect == "synonymous"),
                                        .groups = "drop"),
                     by = "gene_id") %>%
    # Perform tests
    purrr::pmap_dfr(function(gene_id, n_genesing, n_genesing_ns, n_genesing_s,
                             n_genesnvs, n_genesnvs_ns, n_genesnvs_s,
                             tot_snvs, tot_sing){
      
      res1 <- fisher.test(matrix(c(n_genesing, tot_sing,
                                   n_genesnvs, tot_snvs),
                                 ncol = 2))
      Ns <- n_genesnvs_s - n_genesing_s
      Nn <- n_genesnvs_ns - n_genesing_ns
      
      res2 <- fisher.test(matrix(c(n_genesing_ns, Nn,
                                   n_genesing_s, Ns),
                                 ncol = 2))
      
      tibble::tibble(gene_id = gene_id,
                     # n_genesing = n_genesing,
                     # n_genesnvs = n_genesnvs,
                     # tot_sing = tot_sing,
                     # tot_snvs = tot_snvs,
                     # Sn = n_genesing_ns,
                     # Ss = n_genesing_s,
                     Nn = Nn,
                     Ns = Ns,
                     snv_OR = res1$estimate,
                     snv_pval = res1$p.value,
                     sing_OR = res2$estimate,
                     sing_pval = res2$p.value)
      
    }, tot_snvs = sum(!is.na(info$gene_id)),
    tot_sing = nrow(Sing))
}
