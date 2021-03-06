% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequence.r
\name{get_sequence_fragment}
\alias{get_sequence_fragment}
\title{Get sequence fragment}
\usage{
get_sequence_fragment(seq, ref_id, start, end)
}
\arguments{
\item{seq}{A sequence fasta object from \link{seqinr}}

\item{ref_id}{The sequence ID from the fragment,}

\item{start}{1-indexed start position.}

\item{end}{1-indexed end position (included in fragment).}
}
\value{
A character vector
}
\description{
Internal. Simple wrapper for \link[seqinr]{getFrag}
}
