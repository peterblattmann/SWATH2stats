\name{count_analytes}
\alias{count_analytes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
count_analytes: Counts analytes in different injections
}
\description{
This functions counts the number of different peakgroups, peptides and proteins in different injections
}
\usage{
count_analytes(data, column.levels = c("transition_group_id", "FullPeptideName",
"ProteinName"), column.by="run_id", rm.decoy=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing SWATH data.
}
  \item{column.levels}{
Columns in which different identifiers should be counted.
}
  \item{column.by}{
Column for which the  different identifiers should be counted for, e.g. for the different injections.
}
  \item{rm.decoy}{
Option to not remove decoy before counting.
}

}
\value{
Returns a data frame with the count of the different identifiers per e.g. injection.
}

\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
count_analytes(data)
}

