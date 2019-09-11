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

// Parameters
params.genomes = ''
params.midas_dir = ''
params.metwas_dir = ''
params.outdir = 'output/'
params.metawas_suffix = "_lmm.assoc.txt"

// Prepare files
genomes = file(params.genomes)
reader = genomes.newReader()
GENOMES = []
// INTERVALS = []
while(str = reader.readLine()){
  GENOMES = GENOMES + [tuple(str,
    file("${params.metawas_dir}/${str}${params.metwas_suffix}"),
    file("${params.midas_dir}/${str}/snps_info.txt"))]
}

process prop_sig{
  label 'r'
  publishDir "${params.outdir}/plots/scatter", mode: 'rellink',
    pattern: "freq_vs_pval.png", saveAs: {"${genome}.png"}
  publishDir "${params.outdir}/plots/facets", mode: 'rellink',
    pattern: "freq_vs_pval.facets.png", saveAs: {"${genome}.png"}
  publishDir "${params.outdir}/counts", mode: 'rellink',
    pattern: "sig_nums.txt", saveAs: {"${genome}.txt"}

  input:
  set genome, metawas, info from GENOMES

  output:
  file "freq_vs_pval.png"
  file "freq_vs_pval.facets.png"
  file "sig_nums.txt"

  """
  ${workflow.projectDir}/prop_significant.r $metawas $info
  """
}
