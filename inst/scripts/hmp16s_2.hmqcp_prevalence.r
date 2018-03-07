library(HMAR)

# Mapping file is the same everywhere
Map <- read.table(file = "~/micropopgen/data/hmp_16S/HMQCP/v13_map_uniquebyPSN.txt.bz2",
                  sep = "\t", header = TRUE, comment.char = "", row.names = 1)
row.names(Map) <- paste("X", row.names(Map), sep = "")
head(Map)

Tab <- read.am(file = "~/micropopgen/data/hmp_16S/HMQCP/otu_table_psn_v13.txt.gz",
               format = 'qiime', taxonomy = "Consensus.Lineage")
# head(Tab)
Tab$Tab[1:5,1:5]
head(Tab$Map)
head(Tab$Tax)

to_remove <- setdiff(samples(Tab), row.names(Map))
Tab <- remove_samples(Tab, samples = to_remove, droplevels = TRUE)

Map[ setdiff(samples(Tab), row.names(Map)), ]

Dat <- create_dataset(Tab = Tab$Tab, Map = Map[ samples(Tab), ], Tax = Tab$Tax )
Dat <- subset(Dat, HMPbodysubsite %in% c("Buccal_mucosa", "Supragingival_plaque", "Tongue_dorsum", "Stool"),
              drop = TRUE, clean = TRUE)
Dat

prev <- calculate_prevalence(Dat = Dat, thres = 1, group = "HMPbodysubsite")
head(prev)

# Sort
prev <- prev[ order(prev$Group, prev$Proportion , decreasing = TRUE), ]
prev$Taxon <- factor(prev$Taxon, unique(prev$Taxon))