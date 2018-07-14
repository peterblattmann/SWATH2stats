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
                 id.separator = "/",
                 copy_nonconverted = TRUE,
                 verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_table}{
  A data frame or file name.
}
  \item{column.name}{
  The column name where the original protein identifiers are present. Default: "Protein".
}
  \item{species}{
  The species of the protein identifiers in the term used by biomaRt (e.g. "hsapiens_gene_ensembl", "mmusculus_gene_ensembl", "drerio_gene_ensembl", etc.). Default: "hsapiens_gene_ensembl".
}
  \item{host}{
  Path of the biomaRt database (e.g. "www.ensembl.org", "dec2017.archive.ensembl.org"). Default: "www.ensembl.org".
}
  \item{mart}{
  The type of mart (e.g. "ENSEMBL_MART_ENSEMBL", etc.). Default: "ENSEMBL_MART_ENSEMBL"
}
  \item{ID1}{
  The type of the original protein identifiers (e.g. "uniprotswissprot", "ensembl_peptide_id"). Default: "uniprotswissprot".
}
  \item{ID2}{
  The type of the converted protein identifiers (e.g. "hgnc_symbol", "mgi_symbol", "external_gene_name"). Default: "hgnc_symbol".
}
  \item{id.separator}{
  Separator between protein identifiers of shared peptides. Default: "/". 
}
  \item{copy_nonconverted}{
  Option defining if the identifiers that cannot be converted should be copied. Default: TRUE.
}

  \item{verbose}{
  Option to write a file containing the version of the database used. Default: FALSE.
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
