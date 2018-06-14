\name{load_mart}
\alias{load_mart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
load_mart: Load a biomart database
}
\description{
This function loads a biomart database
}
\usage{
load_mart(species, ensembl.path, mart, verbose = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species}{
  The species of the protein identifiers in the term used by biomaRt (e.g. "hsapiens_gene_ensembl", "mmusculus_gene_ensembl", "drerio_gene_ensembl", etc.)
}
  \item{ensembl.path}{
  Path of the biomaRt database. 
}
  \item{mart}{
  The species of the protein identifiers in the term used by biomaRt (e.g. "hsapiens_gene_ensembl", "mmusculus_gene_ensembl", "drerio_gene_ensembl", etc.)
}
\item{verbose}{
  Option to write a file containing the version of the database used. 
}
 
}
\value{
Returns the mart from biomaRt.
}
\note{
Protein identifiers from shared peptides should be separated by a forward slash. The host of archived ensembl databases can be introduced as well (e.g. "dec2017.archive.ensembl.org") 
}
\author{
Peter Blattmann
}
\examples{
organism.mart <- load_mart("hsapiens_gene_ensembl", "www.ensembl.org", "ENSEMBL_MART_ENSEMBL")
class(organism.mart)
}
