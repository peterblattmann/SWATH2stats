\name{filter_mscore}
\alias{filter_mscore}
\alias{filter_mscore_freqobs}
\alias{filter_mscore_condition}
\title{
filter_mscore: Filter openSWATH output table according to mscore
}
\description{
This function filters the SWATH data according to the m_score value, as well as to the number of occurence in the data (requant) and within a condition (condition)
}
\usage{
filter_mscore(data, mscore, rm.decoy=TRUE, mscore.col = "m_score")
filter_mscore_freqobs(data, mscore, percentage=NULL, 
rm.decoy = TRUE, mscore.col = "m_score")
filter_mscore_condition(data, mscore, n.replica, rm.decoy = TRUE, 
mscore.col = "m_score")
}
\arguments{
  \item{data}{
A data frame containing SWATH data.
}
  \item{mscore}{
Value that defines the mscore threshold according to which the data will be filtered.
}
  \item{n.replica}{
Number of measurements within at least one condition that have to pass the mscore threshold for this transition.
}
  \item{percentage}{
Percentage in which replicas the transition has to reach the mscore threshold
}
  \item{rm.decoy}{
Option to remove the decoys during filtering.
}
  \item{mscore.col}{
Defines the column from which to retrieve the m_score. If you use JPP (Rosenberger, Bludau et al. 2017) this can be used to select between Protein and transition_group m_score.
}

}
\value{Returns a data frame with the filtered data.}

\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.filtered <- filter_mscore(data, 0.01)
data.filtered <- filter_mscore_freqobs(data, 0.01, 0.8)
data.filtered <- filter_mscore_condition(data, 0.01, 3)

}
