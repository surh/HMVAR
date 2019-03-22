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

Tab <- read.table("mktest_fraserv/Haemophilus_parainfluenzae_57123/mk_results.Buccal.mucosa_Supragingival.plaque.txt",
                  sep = "\t", header = TRUE)
head(Tab)

summary(Tab$hg)
summary(Tab$ratio)
summary(Tab$ratio)


plot(Tab$hg)

plot(Tab$hg.pval, Tab$ratio.pval)
hist(Tab$hg.pval, breaks = 20)
hist(Tab$ratio.pval, breaks = 20)


tab <- subset(Tab, Ds > 0 & Pn > 0)
hist(tab$ratio.pval, breaks = 20)
summary(qvalue::qvalue(tab$ratio.pval))
tab[ order(tab$ratio.pval, decreasing = FALSE), ]

tests <- apply(tab,1,function(x, a){
  t <- as.numeric(x[5:8])
  t <- matrix(t, ncol = 2, byrow = TRUE)
  # print(t)
  f <- fisher.test(t, alternative = a)
  
  res <- data.frame(Dn = as.numeric(x[5]),
                    Ds = as.numeric(x[6]),
                    Pn = as.numeric(x[7]),
                    Ps = as.numeric(x[8]),
                    ratio = f$estimate,
                    pval = f$p.value)
  
  return(res)
}, a = 'greater')
tests <- do.call(rbind, tests)
tests

(tests$Dn * tests$Ps) / (tests$Ds * tests$Pn)


plot(tab$ratio.pval, tests$pval)
