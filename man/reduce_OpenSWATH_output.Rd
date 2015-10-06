\name{reduce_OpenSWATH_output}
\alias{reduce_OpenSWATH_output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reduce columns of OpenSWATH data
}
\description{
This function selects the columns from the standard OpenSWATH output to column needed for MSstats, aLFQ and mapDIA.
}
\usage{
reduce_OpenSWATH_output(data, column.names=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing SWATH data.
}
  \item{column.names}{
A vector of column names that can be selected.
}
}
\value{
Returns a data frame with the selected columns.
}
\note{
A basic set of columns are defined in the function and are used if no column names are indicated
}

\author{
Peter Blattmann
}
\note{
The column.names can be omitted and then the following columns are selected that are needed for MSstats and mapDIA analysis: ProteinName, FullPeptideName, Sequence, Charge, aggr_Fragment_Annotation, aggr_Peak_Area, align_origfilename, m_score, decoy, Intensity, RT.
This function should be ommitted if the data is analyzed afterwards with the aLFQ or imsbInfer package that needs further columns. 
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.filtered <- reduce_OpenSWATH_output(data) 
}
