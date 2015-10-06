\name{plot.fdr_table}
\alias{plot.fdr_table}
\title{
Plot functionality for results of class "fdr_table" as produced by e.g. the function assess_fdr_overall()
}
\description{
This function created standard plots from results of class "fdr_table" as produced by e.g. the function assess_fdr_overall() visualizig ID numbers in dependence of estimated FDR and also estimated FDR in dependence of m_score cutoff.
}
\usage{
\method{plot}{fdr_table}(x, output = "Rconsole", filename = "FDR_report_overall", \dots)
}
\arguments{
  \item{x}{
List of class "fdr_table" as produced e.g. by the function assess_fdr_overall() from this package.
}
  \item{output}{
Choose output type. "pdf_csv" creates the output as files in the working directory, "Rconsole" triggers delivery of the output to the console enabling further computation or custom plotting / output.
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
data("OpenSWATH_data_FDR", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data_FDR, Study_design) 
x <- assess_fdr_overall(data, FFT=0.7, output = "Rconsole", plot = FALSE)
plot.fdr_table(x, output = "pdf_csv", filename = "Assess_fdr_overall_testplot")
}
