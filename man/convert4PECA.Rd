% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert4PECA.R
\name{convert4PECA}
\alias{convert4PECA}
\title{Convert table into the format for ROPECA}
\usage{
convert4PECA(data)
}
\arguments{
\item{data}{A data frame containing SWATH data.}
}
\value{
Returns a data frame in the appropriate format for ROPECA.
}
\description{
This functions selects the columns necessary for ROPECA.
}
\note{
The table must not contain any technical replica, the intensity of
  technical replica is averaged. This function requires the package reshape2.
}
\examples{
{
 data("OpenSWATH_data", package="SWATH2stats")
 data("Study_design", package="SWATH2stats")
 data <- sample_annotation(OpenSWATH_data, Study_design)
 data.filtered.decoy <- filter_mscore(data, 0.01)
 data.PECA <- convert4PECA(data.filtered.decoy)
}
}
\references{
Suomi, T. and Elo L.L. (2017). "Enhanced differential expression
  statistics for data-independent acquisition proteomics" Scientific Reports
  7, Article number: 5869.doi:10.1038/s41598-017-05949-y
}
\author{
Peter Blattmann
}
