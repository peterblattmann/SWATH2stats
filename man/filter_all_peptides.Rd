\name{filter_all_peptides}
\alias{filter_all_peptides}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select all proteins that are supported by peptides.
}
\description{
This functions counts all proteins that are supported by peptides (including non proteo-typic peptides). All peptides (incl. non proteotypic peptides are selected. For the proteins supproted by proteotypic peptide the "1/" in front of the identifier is removed to facilitate further data processing. 
}
\usage{
filter_all_peptides(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing SWATH data.
}
}
\value{Returns a data frame with the data from both proteotypic and non-proteotypic peptides.}

\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.filtered.decoy <- filter_mscore(data, 0.01) 
data.all <- filter_all_peptides(data.filtered.decoy)
}
