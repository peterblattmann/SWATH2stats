\name{convert4mapDIA}
\alias{convert4mapDIA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
convert4mapDIA: Convert table into the format for mapDIA
}
\description{
This functions selects the columns necessary for mapDIA.
}
\usage{
convert4mapDIA(data, RT=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  A data frame containing SWATH data. 
}
  \item{RT}{
  Option to export the retention times.  
}

}
\value{
Returns a data frame in the appropriate format for mapDIA.
}
\references{
Teo et al. submitted, mapDIA is available on Sourceforge \url{http://sourceforge.net/projects/mapdia/}
}

\note{
The table must not contain any technical replica, the intensity of technical replica is averaged. This function requires the package reshape2. 
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
data.mapDIA <- convert4mapDIA(raw, RT=TRUE)
}
