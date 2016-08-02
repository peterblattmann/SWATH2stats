\name{transform_MSstats_OpenSWATH}
\alias{transform_MSstats_OpenSWATH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
transform_MSstats_OpenSWATH: Transforms column names to OpenSWATH column names
}
\description{
This functions transforms the column names from a data frame in MSstats format to a data frame with column names used by the OpenSWATH output. The original table needs to contain at least the 10 columns defined by MSstats: ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity.)
}
\usage{
transform_MSstats_OpenSWATH(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing the SWATH data in the MSstats format
}

}
\value{
Returns the data frame in the appropriate format.
}
\references{
Choi M, Chang CY, Clough T, Broudy D, Killeen T, MacLean B, Vitek O. MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments.Bioinformatics. 2014 Sep 1;30(17):2524-6. doi: 10.1093/bioinformatics/btu305.
}

\author{
Peter Blattmann
}
\examples{
data("MSstats_data", package="SWATH2stats")
transform_MSstats_OpenSWATH(MSstats_data)
}
