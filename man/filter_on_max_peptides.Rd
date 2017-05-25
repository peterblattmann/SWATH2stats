\name{filter_on_max_peptides}
\alias{filter_on_max_peptides}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter only for the highest intense peptides
}
\description{
In order to reduce the data, the data is filtered only for the proteins with the highest intensity peptides.
}
\usage{
filter_on_max_peptides(data, n_peptides, rm.decoy = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing SWATH data with the column names: ProteinNames, PeptideSequence, PrecursorCharge, Intensity.
}
  \item{n_peptides}{
Maximum number of highest intense peptides to filter the data on.
}
  \item{rm.decoy}{
Option to remove the decoys during filtering.
}

}
\value{
Returns a data frame of the filtered data
}
\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.filtered <- filter_mscore_freqobs(data, 0.01,0.8)
data.max <- filter_on_max_peptides(data.filtered, 5)
}
