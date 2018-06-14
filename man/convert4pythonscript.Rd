\name{convert4pythonscript}
\alias{convert4pythonscript}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
convert4bashscript: Convert data into the format for running a bash script
}
\description{
This functions selects the columns suggested to run a python script to change the data from peptide-level to transition-level.
}
\usage{
convert4pythonscript(data, replace.Unimod = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing SWATH data.
}
  \item{replace.Unimod}{
Option to indicate if Unimod Identifier should be replaced form ":"" to "_".
}
}
\details{
The necessary columns are selected and the run column is renamed to filename for the script. The intensities are taken from the column aggr_Peak_Area and therefore the Intensity column is not exported. 

}
\value{
Returns a data frame in the appropriate format to be used by a custom python script stored in the scripts folder.
}

\author{
Peter Blattmann
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.filtered.decoy <- filter_mscore(data,0.01) 
data.pythonscript <- convert4pythonscript(data.filtered.decoy)
}

