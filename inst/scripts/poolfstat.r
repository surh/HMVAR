readcoverage <- (Dat$depth[1:10000,] %>%
                   select(-site_id) %>%
                   as.matrix())
refallele.readcount <- readcoverage - round((Dat$freq[1:10000,] %>%
                                               select(-site_id) %>%
                                               as.matrix()) * readcoverage)

readcoverage
refallele.readcount

groups <- setNames(map$Group, nm = map$sample)
groups <- groups[colnames(readcoverage)]
readcoverage <- AMOR::pool_samples(readcoverage, groups)
refallele.readcount <- AMOR::pool_samples(refallele.readcount, groups)

Dat.pooldata <- new("pooldata",
                    npools = 4,
                    nsnp = 10000,
                    refallele.readcount = refallele.readcount,
                    readcoverage = readcoverage,
                    snp.info = Dat$info[1:10000,] %>% select(ref_id, ref_pos, major_allele, minor_allele) %>% as.matrix(),
                    poolsizes = table(groups) %>% as.vector(),
                    poolnames = colnames(readcoverage))
FSTpool <- computeFST(pooldata = Dat.pooldata, method = "Anova")