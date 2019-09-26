#!/usr/bin/env nextflow
// Copyright (C) 2019 Sur Herrera Paredes

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Takes output of midas merge and metwas and calculates the proportion
// of significant results for snp type.

// Params
params.indir = ""
params.map = "map.txt"
params.method = "Fstpool"
params.outdir = "output"

println params.indir
INDIRS = Channel.fromPath("${params.indir}/*", type: 'dir')

process fst{
  label 'r'
  publishDir params.outdir, saveAs: 'rellink'

  input:
  file indir from INDIRS
  val method from params.method
  file map from params.map

  output:
  file "out/*.fst.txt"

  // exec:
  // printn indir
  """
  ${workflow.projectDir}/calculate_fst.r \
    $indir \
    --map_file $map \
    --outdir out/
    --type single
    --method $method
  """
}
