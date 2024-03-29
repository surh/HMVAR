% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.r
\name{format_input}
\alias{format_input}
\alias{format_input_hmmcp}
\title{Formats input form HMMCP 16S files}
\usage{
format_input(
  name,
  counts_file,
  taxonomy_file,
  Map = Map,
  collapse_level = NULL
)
}
\arguments{
\item{name}{name(not used anymore?)}

\item{counts_file}{A file with counts}

\item{taxonomy_file}{A file with taxonomy}

\item{Map}{A mapping data.frame}

\item{collapse_level}{Level to collapse taxonomy}
}
\value{
A Dataset object
}
\description{
Only to be used to process 16S files from MOTHUR \
pipeline of HMP (HMMCP). Count files must be pre-edited
to remove trailing tab
}
\author{
Sur Herrera Paredes
}
