# 0.1-1
1. Overall cleaning and NAMESPACE update
2. Adding focal group map pre-processing option for
midas_mktest function and mktest,r executable script.
3. Flipping mice imputation matrix. Added support for blocks.
NEED TO DOCUMENT block_size
4. Reverting to standard alphanumeric names for chromosomes
in BIMBAM format.
5. Eliminating annoying messages from read_midas_abun.
6. Adding toy example in MIDAS format (Issue
https://github.com/surh/HMVAR/issues/6).
7. Updating `read_midas_data` function, and adding examples.
8. Function `midas_mktest` now accepts either a file path to map
or a data frame witht the map. It also has examples.
9. Examples for functions that calculate MK contingency table.
10. Added unit testing support via testhat (Issue
https://github.com/surh/HMVAR/issues/25).
11. Fixed allele asignment by `determine_site_dist` (Issue
https://github.com/surh/HMVAR/issues/1).
12. Unit testing for functions that obtain MK contingency table.
13. Adding `determine_substitution_type`.
14. Adding `match_freq_and_depth`.
15. Adding `plotgg_sacked_columns`.
16. Added `variable_dist_per_site`.
17. Added `varsites_pipeline`. Updated `varsites.r` executable.
to use this function. (Issue https://github.com/surh/HMVAR/issues/4).
18. Added `plotgg_manhattan` function. (Issue https://github.com/surh/HMVAR/issues/10).
19. Fixed various issues with NAMESPACE. Including removing AMOR from
depends field (Issue https://github.com/surh/HMVAR/issues/13).
20. Generalizing column names for genomewide_ni.
21. Added `calculate_mktable`.
22. Added functions for DoS (Issue https://github.com/surh/HMVAR/issues/11).
23. Created `dos.r` executable (Issue
https://github.com/surh/HMVAR/issues/12).
24. Gene-set enrichment analysis (`gsea`) function created
(Issue https://github.com/surh/HMVAR/issues/18).
25. Function `test_go` added (Issue https://github.com/surh/HMVAR/issues/16
https://github.com/surh/HMVAR/issues/14).
26. Function `terms_enrichments`.
27. `metawas_enrichments.r` script converted into executable.
`annotation_enrichments.r`. Issue (https://github.com/surh/HMVAR/issues/15).
28. Adding `sign_test` function. Issue (https://github.com/surh/HMVAR/issues/19).

# 0.1-0
1. Added exectuable script to perform mktest.
2. determine_snp_effect can handle non-coding loci.
3. Added executable to calculate and plot the number of
variable and fixed positions within sample (Issue
https://github.com/surh/HMVAR/issues/4)
4. Added executable to subset SNPs from MIDAS after merging 
(Issue https://github.com/surh/HMVAR/issues/3).
6. read_midas_data now can keep non CDS positions.
(Issue https://github.com/surh/HMVAR/issues/7)
7. Added functions for metawas steps: midas_to_bimbam,
gemma_kinship, and gemma_lmm.
8. Created util set of functions for internal utility functions.
9. Added mice_impute() function
10. Added benchmark imputation function and executable script.
11. Added genome_metawas.r executable script.

# 0.0-4
1. Added ggtree import
2. Changed name to HMVAR
3. Inititiated moving to tidyverse
4. Updated genome_wide_ni
5. Fising documentation issues
6. calculate_neutrality_index.r has been updated.
7. Added pval_qqplot functions and documentation.
8. Recoment svglite package
9. Added functions for mktest from MIDAS output

# 0.0-3
1. Addded script to test strains for overall SNP
differentiation between aevery pair of comparisons.
2. Added script and functions for genome wide neutrality index
calculations.
3. Functions and script for checking and plotting p-value distribution
4. Numerous scripts and functions for evaluating mktest
results

# 0.0-2
1. Inputing HMQC data
2. Prevalence of HMQC data

# 0.0-1
1. First version.
2. Prevalence in HMMCP dataset
