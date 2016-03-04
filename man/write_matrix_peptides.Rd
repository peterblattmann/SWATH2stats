\name{write_matrix_peptides}
\alias{write_matrix_peptides}
\title{
write_matrix_peptides: Writes out an overview matrix of peptides mapping to a FDR quality controlled protein master list at controlled global peptide FDR quality.
}
\description{
Writes out an overview matrix on peptide level of a supplied (unfiltered or prefiltered) OpenSWATH results data frame.
The peptide quantification is achieved by summing the areas under all 6 transitions per precursor and summing all precursors per FullPeptideName.
In order to keep the peptide-to-protein association, the FullPeptideName is joined with the ProteinName.
}
\usage{
write_matrix_peptides(data, write.csv=FALSE, 
filename = "SWATH2stats_overview_matrix_peptidelevel.csv",
rm.decoy = FALSE)
}
\arguments{
  \item{data}{
A data frame containing annotated OpenSWATH/pyProphet data.
}
  \item{write.csv}{
Option to determine if table should be written automatically into csv file.
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
write_matrix_peptides(data)
}
