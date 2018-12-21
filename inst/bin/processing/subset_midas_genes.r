#!/usr/bin/env Rscript

# (C) Copyright 2018 Sur Herrera Paredes
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

library(HMVAR)
library(tidyverse)
library(argparser)


args <- list(midas_dir = "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/Leptotrichia_shahii_58912/",
             genes = "/godot/users/sur/exp/fraserv/2018/today2/selected_genes.txt",
             outdir = "out/")


