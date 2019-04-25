\name{disaggregate}
\alias{disaggregate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
disaggregate: Transforms the SWATH data from a peptide- to a transition-level table.
}
\description{
If the SWATH data should be analyzed on transition-level the data needs to be tranformed from peptide-level table to a transition-level table (one row per transition instead of one row per peptide). The columns "aggr_Fragment_Annotation" and "aggr_Peak_Area" are disaggregated into the new columns "FragmentIon" and "Intensity".
The following columns are renamed if they exist: FullPeptideName -> PeptideSequence, Charge -> PrecursorCharge, Area -> Intensity, Fragment -> FragmentIon, Sequence -> NakedSequence. 
}
\usage{
disaggregate(data, all.columns = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing SWATH data.
}
  \item{all.columns}{
Option that all columns are processed. Otherwise only the typical columns needed for downstream analysis are processed. 
}
}
\value{Returns a data frame containing the SWATH data in a transition-level table.}

\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.filtered.decoy <- filter_mscore(data, 0.01) 
raw <- disaggregate(data.filtered.decoy)
}
