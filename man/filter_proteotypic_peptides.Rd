\name{filter_proteotypic_peptides}
\alias{filter_proteotypic_peptides}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter for proteins that are supported by proteotypic peptides.
}
\description{
Peptides can match to several proteins. With this function proteotypic peptides, peptides that are only contained in one protein are selected. Additionally the number of proteins are counted and printed.
}
\usage{
filter_proteotypic_peptides(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing SWATH data.
}
}
\value{Returns a data frame with only the data supported by proteotypic peptides.}

\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.filtered.decoy <- filter_mscore(data, 0.01) 
data.all <- filter_proteotypic_peptides(data.filtered.decoy)
}
