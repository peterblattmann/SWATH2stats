\name{convert4MSstats}
\alias{convert4MSstats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
convert4MSstats: Convert table into the format for MSstats
}
\description{
This functions selects the columns necessary for MSstats and renames them if necessary.
}
\usage{
convert4MSstats(data, replace.values = TRUE, replace.colnames = TRUE, 
replace.Unimod = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing SWATH data.  
}
  \item{replace.values}{
Option to indicate if negative and 0 values should be replaced with NA.
}
  \item{replace.colnames}{
Option to indicate if column names should be renamed and columns reduced to the necessary columns for MSstats
}
  \item{replace.Unimod}{
Option to indicate if Unimod Identifier should be replaced from ":" to "_".
}

}
\details{
The necessary columns are selected and three columns renamed: 
FullPeptideName -> PeptideSequence
Charge -> PrecursorCharge
filename -> File 

}
\value{
Returns a data frame in the appropriate format for MSstats.
}
\references{
Choi M, Chang CY, Clough T, Broudy D, Killeen T, MacLean B, Vitek O. MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments.Bioinformatics. 2014 Sep 1;30(17):2524-6. doi: 10.1093/bioinformatics/btu305.
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
data.mapDIA <- convert4MSstats(raw)
}

