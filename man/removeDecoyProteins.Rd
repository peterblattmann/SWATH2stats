\name{removeDecoyProteins}
\alias{removeDecoyProteins}
\alias{rmDecoyProt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Removes decoy proteins from the protein group label
}
\description{
There exist peptides annotated as protein groups with 2/ProteinA/DECOY_ProteinB. However these are in principal proteotypic peptides and should be annoated 1/ProteinA. This function changes these labels accordingly. The subfunction rmDecoyProt removes the Decoy protein, calling removeDecoyProteins also changes the nubmer before the protein group accordingly. 
}

\usage{
removeDecoyProteins(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing SWATH data.
}
}
\value{Returns a data frame with changed protein labels.}

\author{
Moritz Heusel
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.filtered.decoy <- filter_mscore(data, 0.01) 
data.2 <- removeDecoyProteins(data.filtered.decoy)
}
