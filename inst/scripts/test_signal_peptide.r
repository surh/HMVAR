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



#########################################


pos <- read.table(file = "gff_grampos.txt", sep = "\t", stringsAsFactors = FALSE)
head(pos)
neg <- read.table(file = "gff_gramneg.txt", sep = "\t", stringsAsFactors = FALSE)
head(neg)

signalp <- c(neg$V1, pos$V1)
signalp <- unique(signalp)

Tab <- read.table(file = "signifincant_mktest", header = TRUE)
head(Tab)
ftable(A ~ B, Tab)

Tab$signalP <- "No"
Tab$signalP[ Tab$gene %in% signalp ] <- "Yes"
table(Tab$signalP)

p1 <- ggplot(Tab, aes(x = signalP, y = log2(ratio))) +
  # geom_boxplot() +
  geom_violin()
p1

# summary(lm(log2(ratio) ~ signalP, data = Tab))
