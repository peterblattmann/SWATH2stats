\name{plot.fdr_cube}
\alias{plot.fdr_cube}
\title{
Plot functionality for FDR assessment result arrays as produced by e.g. the function assess_fdr_byrun()
}
\description{
This function creates standard plots from result arrays as produced by e.g. the function assess_fdr_byrun(), visualizig assay, peptide and protein level FDR for each run at m-score cutoffs 1e-2 and 1e-3. Furthermore, Target and Decoy ID numbers are visualized.
}
\usage{
\method{plot}{fdr_cube}(x, output = "Rconsole", filename = "FDR_report_byrun", \dots)
}
\arguments{
  \item{x}{
Array of by-run FDR assessment results as produced e.g. by the function assess_fdr_byrun() from this package.
}
  \item{output}{
Choose output type. "pdf_csv" creates the output as files in the working directory, "Rconsole" triggers delivery of the output to the console enabling further computation and/or custom plotting / output.
}
  \item{filename}{
Basename for output files to be created (if output = "pdf_csv" has been selected).
}
  \item{\dots}{
further arguments passed to method.
}
}
\value{
Plots in Rconsole or report files.
}
\author{
Moritz Heusel
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
x <- assess_fdr_byrun(data, FFT=0.7, output = "Rconsole", plot = FALSE)
plot.fdr_cube(x, output = "pdf_csv", filename = "Assess_fdr_byrun_testplot")
}
