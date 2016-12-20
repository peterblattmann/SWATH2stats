\name{unifyProteinGroupLabels}
\alias{unifyProteinGroupLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unify the protein group labels
}
\description{
Unify the protein group labels (2/ProteinA/ProteinB and 2/ProteinB/ProteinA) to one common label (e.g. 2/ProteinA/ProteinB)}
\usage{
unifyProteinGroupLabels(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing SWATH data.
}
}
\value{Returns a data frame with the unififed protein labels.}

\author{
Moritz Heusel
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.filtered.decoy <- filter_mscore(data, 0.01) 
data.unified <- unifyProteinGroupLabels(data.filtered.decoy)
}
