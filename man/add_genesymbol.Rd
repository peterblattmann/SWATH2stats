\name{add_genesymbol}
\alias{add_genesymbol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
add_genesymbol: Adds gene symbols to a table
}
\description{
This function adds an additional gene or protein symbol to a table.
}
\usage{
add_genesymbol(data_table, gene.ID.table, column.name = "Protein", 
ID1 = "uniprotswissprot", ID2 = "hgnc_symbol", id.separator = "/", copy_nonconverted = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_table}{
  A data frame or file name.
}
  \item{gene.ID.table}{
  A dictionary table containing the original and new identifiers.
}
  \item{column.name}{
  The column name where the original protein identifiers are present. 
}
  \item{ID1}{
  The type of the original protein identifiers.
}
  \item{ID2}{
  The type of the converted protein identifiers.
}
  \item{id.separator}{
  Separator between protein identifiers of shared peptides.
}
  \item{copy_nonconverted}{
  Option defining if the identifiers that cannot be converted should be copied.
}


}
\value{
Returns the data frame with an added column of the converted protein identifiers.
}
\note{
Protein identifiers from shared peptides should be separated by a forward slash.
}
\author{
Peter Blattmann
}
\examples{
gene.ID.table <- data.frame(uniprotswissprot = c("Q01581", "P49327", "P60709"), 
hgnc_symbol = c("HMGCS1", "FASN", "ACTB"))
data_table <- data.frame(Protein = c("Q01581", "P49327", "2/P63261/P60709"), 
Abundance = c(100, 3390, 43423))
add_genesymbol(data_table, gene.ID.table)
}
