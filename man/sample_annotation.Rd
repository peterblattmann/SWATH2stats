\name{sample_annotation}
\alias{sample_annotation}
\title{
sample_annotation: Annotate the SWATH data with the sample information
}
\description{
For statistical analysis and filtering the measurements need to be annotated with Filename, Condition, BioReplicate, and Run. This functions takes this information from a txt file containing this meta-data.
}
\usage{
sample_annotation(data, sample.annotation, data.type="openSWATH", 
column.file = "align_origfilename", change.run.id = TRUE, verbose=FALSE)
}
\arguments{
  \item{data}{
A data frame containing SWATH data.
}
  \item{sample.annotation}{
A data frame containing the columns: Filename, Condition, BioReplicate, Run. The values contained in the column filename have to be present in the filename of the SWATH data.
}
\item{data.type}{
Option to specify the format of the table, if the column names from an OpenSWATH output or MSstats table are used.
}
\item{column.file}{
Option to specify the column name where the injection file is specified. Default is set to "align_origfilename".
}
\item{change.run.id}{
Option to choose if the run\_id column shall be reassigned to a unique value combining the values of Condition, BioReplicate and Run. (Option only possible if data is of format "OpenSWATH")
}

\item{verbose}{
Option to turn on reporting on which filename it is working on.
}
}
\value{
Returns a dataframe with each row annotated for the study design
}

\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
}
