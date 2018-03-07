library(HMAR)

# List input files
files <- data.frame(Name = c("hmqcp.v13.otu",
                             "hmqcp.v35.otu"),
                    counts = c("~/micropopgen/data/hmp_16S/HMQCP/otu_table_psn_v13.txt.gz",
                               "~/micropopgen/data/hmp_16S/HMQCP/otu_table_psn_v35.txt.gz"),
                    map = c("~/micropopgen/data/hmp_16S/HMQCP/v13_map_uniquebyPSN.txt.bz2",
                            "~/micropopgen/data/hmp_16S/HMQCP/v35_map_uniquebyPSN.txt.bz2"),
                    stringsAsFactors = FALSE)
files



i <- 1
Dat <- format_input_hmqcp(name = files$Name[i],
                          counts_file = files$counts[i],
                          map_file = files$map[i],
                          collapse_level = 6)

# Calculate prevalence and change names to homgenize with HMMCP
prev <- calculate_prevalence(Dat = Dat, thres = 1, group = "HMPbodysubsite")
levels(prev$Group) <- c("Buccal mucosa", "Stool", "Supragingival plaque", "Tongue dorsum")
# head(prev)

# Sort
prev <- prev[ order(prev$Group, prev$Proportion , decreasing = TRUE), ]
prev$Taxon <- factor(prev$Taxon, unique(prev$Taxon))

# Plot
p1 <- ggplot(prev, aes(x = Taxon, y = Proportion, group = Group, col = Group )) +
  facet_wrap(~ Group, ncol = 1) +
  geom_line() +
  theme(axis.text.x = element_blank()) +
  theme_blackbox
p1


filename <- paste(files$Name[i], ".prevalence_by_site.svg", sep = "")
cat(filename, "\n")
ggsave(filename, p1, width = 4, height = 8)

filename <- paste(files$Name[i], ".topprev.txt", sep = "")
cat(filename, "\n")
write.table(subset(prev, Proportion >= 0.75), file = filename,
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

filename <- paste(files$Name[i], ".fullprev.txt", sep = "")
cat(filename, "\n")
write.table(prev, file = filename,
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

