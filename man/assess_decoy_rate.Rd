\name{assess_decoy_rate}
\alias{assess_decoy_rate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
assess_decoy_rate: Assess decoy rate
}
\description{
This function counts the number of decoy peptides.
}
\usage{
assess_decoy_rate(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame that contains at least a column named "FullPeptideName" and "decoy". 
}
}
\details{
A printout is generated to indicate the number of non-decoy, decoy peptides and the rate of decoy vs non-decoy peptides. Unique peptides are counted, so a precursor with different charge states is counted as one peptide. In the column "decoy" the values need to be 1,0 or TRUE and FALSE.  
}
\value{
Prints the decoy rate.
}

\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data <- OpenSWATH_data
assess_decoy_rate(data)
}
