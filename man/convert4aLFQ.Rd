\name{convert4aLFQ}
\alias{convert4aLFQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
convert4aLFQ: Convert table into the format for aLFQ
}
\description{
This functions selects the columns necessary for the aLFQ R package.
}
\usage{
convert4aLFQ(data, annotation = TRUE, check_transitions = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing the SWATH data in transition-level format
}
  \item{annotation}{
Option to indicate if the data has been annotated, i.e. if the columns Condition, Replicate, Run are present. If option is set to true it will write a new run_id as a string of the combination of these three columns.
}
  \item{check_transitions}{
  Option if number of transitions should be checked. As input only transition-level data should be used and therefore the number of transitions can be checked to spot any irregulatirities. However, this makes the function slow. 
  }
}
\value{
Returns a data frame in the appropriate format for aLFQ.
}
\references{
Rosenberger G, Ludwig C, Rost HL, Aebersold R, Malmstrom L. aLFQ: an R-package for estimating absolute protein quantities from label-free LC-MS/MS proteomics data. Bioinformatics. 2014 Sep 1;30(17):2511-3. doi: 10.1093/bioinformatics/btu200.
}

\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.filtered.decoy <- filter_mscore(data, 0.01)
raw <- disaggregate(data.filtered.decoy)
data.aLFQ <- convert4aLFQ(raw)
}
