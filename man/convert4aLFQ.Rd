% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert4aLFQ.R
\name{convert4aLFQ}
\alias{convert4aLFQ}
\title{Convert table into the format expected by aLFQ.}
\usage{
convert4aLFQ(data, annotation = TRUE, check_transitions = TRUE)
}
\arguments{
\item{data}{A data frame containing the SWATH data in transition-level
format}

\item{annotation}{Option to indicate if the data has been annotated, i.e. if
the columns Condition, Replicate, Run are present. If option
is set to true it will write a new run_id as a string of the
combination of these three columns.}

\item{check_transitions}{Option if number of transitions should be checked.
As input only transition-level data should be used and
therefore this is checked. However, this makes the function
slow and herewith be omitted.}
}
\value{
Returns a data frame in the appropriate format for aLFQ.
}
\description{
This function selects the columns necessary for the aLFQ R package.
}
\examples{
{
 data("OpenSWATH_data", package="SWATH2stats")
 data("Study_design", package="SWATH2stats")
 data <- SWATH2stats::sample_annotation(OpenSWATH_data, Study_design, verbose=TRUE)
 data.filtered.decoy <- filter_mscore(data, 0.01)
 raw <- disaggregate(data.filtered.decoy)
 data.aLFQ <- convert4aLFQ(raw)
 }
}
\author{
Peter Blattmann
}
