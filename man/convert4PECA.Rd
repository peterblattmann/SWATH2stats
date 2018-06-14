\name{convert4PECA}
\alias{convert4PECA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
convert4PECA: Convert table into the format for ROPECA
}
\description{
This functions selects the columns necessary for ROPECA.
}
\usage{
convert4PECA(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  A data frame containing SWATH data.
}
}
\value{
Returns a data frame in the appropriate format for ROPECA.
}
\note{
The table must not contain any technical replica, the intensity of technical replica is averaged. This function requires the package reshape2. 
}
\references{
Suomi, T. and Elo L.L. (2017). "Enhanced differential expression statistics for data-independent acquisition proteomics" Scientific Reports 7, Article number: 5869.doi:10.1038/s41598-017-05949-y
}
\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.filtered.decoy <- filter_mscore(data, 0.01)
data.PECA <- convert4PECA(data.filtered.decoy)
}
