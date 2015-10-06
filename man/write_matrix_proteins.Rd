\name{write_matrix_proteins}
\alias{write_matrix_proteins}
\title{
write_matrix_proteins: Writes out an overview matrix of summed signals per protein identifier (lines) over run_id(columns). 
}
\description{
Writes out an overview matrix on protein level of a supplied (unfiltered or filtered) OpenSWATH results data frame.
The protein quantification is achieved by summing the areas under all 6 transitions per precursor, summing all precursors per FullPeptideName and all FullPeptideName signals per ProteinName entry.

This function does not select consistently quantified or top peptides but sums all signals availabe that may or may not originate from the same set of peptides across different runs. A more detailed overview can be generated using the function write_matrix_peptides().

Peptide selection can be achieved upstream using e.g. the functions filter_mscore_requant(), filter_on_max_peptides() and filter_on_min_peptides().
}
\usage{
write_matrix_proteins(data, filename = "SWATH2stats_overview_matrix_proteinlevel.csv", 
rm.decoy = FALSE)
}
\arguments{
  \item{data}{
A data frame containing annotated OpenSWATH/pyProphet data.
}
  \item{filename}{
File base name of the .csv matrix written out to the working folder
}
  \item{rm.decoy}{
Logical whether decoys will be removed from the data matrix. Defaults to FALSE. It's sometimes useful to know how decoys behave across a dataset and how many you allow into your final table with the current filtering strategy.
}
}
\value{
No return value, output .csv matrix is written to the working folder.
}
\author{
Moritz Heusel
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
write_matrix_proteins(data)
}
