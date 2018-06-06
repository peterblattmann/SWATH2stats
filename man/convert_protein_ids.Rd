\name{convert_protein_ids}
\alias{convert_protein_ids}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
convert_protein_ids: Convert protein ids
}
\description{
This function renames protein ids in a data frame or file
}
\usage{
convert_protein_ids(data_table, 
                 column.name = "Protein", 
                 species = "hsapiens_gene_ensembl", 
                 host = "www.ensembl.org", 
                 mart = "ENSEMBL_MART_ENSEMBL",
                 ID1 = "uniprotswissprot", 
                 ID2 = "hgnc_symbol", 
                 verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_table}{
  A data frame or file name.
}
  \item{column.name}{
  The column name where the original protein identifiers are present. 
}
  \item{species}{
  The species of the protein identifiers in the term used by biomaRt (e.g. "hsapiens_gene_ensembl", "mmusculus_gene_ensembl", "drerio_gene_ensembl", etc.)
}
  \item{host}{
  Path of the biomaRt database. 
}
  \item{mart}{
  The species of the protein identifiers in the term used by biomaRt (e.g. "hsapiens_gene_ensembl", "mmusculus_gene_ensembl", "drerio_gene_ensembl", etc.)
}
  \item{ID1}{
  The type of the original protein identifiers.
}
  \item{ID2}{
  The type of the converted protein identifiers.
}
  \item{verbose}{
  Option to write a file containing the version of the database used. 
}
 
}
\value{
Returns the data frame with an added column of the converted protein identifiers.
}
\note{
Protein identifiers from shared peptides should be separated by a forward slash. The host of archived ensembl databases can be introduced as well (e.g. "dec2017.archive.ensembl.org") 
}
\author{
Peter Blattmann
}
\examples{
data_table <- data.frame(Protein = c("Q01581", "P49327", "2/P63261/P60709"), 
Abundance = c(100, 3390, 43423))
convert_protein_ids(data_table)
}
