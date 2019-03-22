# (C) Copyright 2019 Sur Herrera Paredes
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

library(HMAR)


# Mapping file is the same everywhere
Map <- read.table(file = "~/micropopgen/data/hmp_16S/HMMCP/pds.metadata.bz2",
                  sep = "\t", header = TRUE)
row.names(Map) <- paste(Map$nap_id, Map$dataset, sep = ".")
# head(Map)

files <- data.frame(Name = c("hmmcp.v13.hq.otu",
                             "hmmcp.v13.hq.phylotype",
                             "hmmcp.v35.hq.otu",
                             "hmmcp.v35.hq.phylotype",
                             "hmmcp.v69.hq.otu",
                             "hmmcp.v69.hq.phylotype"),
                    counts = c("~/micropopgen/data/hmp_16S/HMMCP/hmp1.v13.hq.otu.counts.bz2",
                               "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v13.hq.phylotype.counts.bz2",
                               "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v35.hq.otu.counts.bz2",
                               "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v35.hq.phylotype.counts.bz2",
                               "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v69.hq.otu.counts.bz2",
                               "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v69.hq.phylotype.counts.bz2"),
                    taxonomy = c("~/micropopgen/data/hmp_16S/HMMCP/hmp1.v13.hq.otu.lookup.bz2",
                                 "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v13.hq.phylotype.lookup.bz2",
                                 "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v35.hq.otu.lookup.bz2",
                                 "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v35.hq.phylotype.lookup.bz2",
                                 "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v69.hq.otu.lookup.bz2",
                                 "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v69.hq.phylotype.lookup.bz2"),
                    stringsAsFactors = FALSE)
files

for(i in 1:nrow(files)){
  Dat <- format_input(name = files$Name[i], counts_file = files$counts[i], 
                      taxonomy_file = files$taxonomy[i], Map = Map, collapse_level = 7)
  # Subset sites
  Dat <- subset(Dat, body_site %in% c("Buccal mucosa", "Supragingival plaque", "Tongue dorsum", "Stool"),
                clean = TRUE, drop = TRUE)
  
  prev <- calculate_prevalence(Dat = Dat, thres = 1, group = "body_site")
  head(prev)
  
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
  
  
}


# x <- as.character(Tax$Taxonomy[1:200])
# x
# system.time(phylotype2rdp(x = x))
# system.time(phylotype2rdp2(x = x))
# 
# dat <- NULL
# for(f in seq(from = 100, to = 600, by = 100)){
#   x <- as.character(Tax$Taxonomy[1:f])
#   
#   
#   t1 <- system.time(phylotype2rdp(x = x))
#   t2 <- system.time(phylotype2rdp2(x = x))
#   
#   t1 <- rbind(t1)
#   colnames(t1) <- paste("t1", colnames(t1), sep = ".")
#   row.names(t1) <- as.character(f)
#   
#   t2 <- rbind(t2)
#   colnames(t2) <- paste("t2", colnames(t2), sep = ".")
#   row.names(t2) <- as.character(f)
#   
#   dat <- rbind(dat,cbind(t1, t2))
# }
# dat





