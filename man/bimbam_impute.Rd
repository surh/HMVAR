% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metawas.r
\name{bimbam_impute}
\alias{bimbam_impute}
\title{Impute genotypes with BIMBAM}
\usage{
bimbam_impute(
  geno_file,
  pheno_file,
  pos_file,
  bimbam = "bimbam",
  outdir = "imputed/",
  em_runs = 10,
  em_steps = 20,
  em_clusters = 15,
  prefix = "imputed"
)
}
\arguments{
\item{geno_file}{Path to mean genotype file. See BIMBAM docummentation
for details.}

\item{pheno_file}{Path to phenotype file. See BIMBAM docummentation
for details.}

\item{pos_file}{Path to SNP position file. See BIMBAM docummentation
for details.}

\item{bimbam}{BIMBAM executable.}

\item{outdir}{Output directory.}

\item{em_runs}{Number of EM algorithm runs.}

\item{em_steps}{Steps of each EM run.}

\item{em_clusters}{Number of clusters in EM algorithm.}

\item{prefix}{Prefix of all output files}
}
\description{
Impute genotypes with BIMBAM
}
\references{
http://www.haplotype.org/download/bimbam-manual.pdf
}
