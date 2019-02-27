\name{filter_mscore_fdr}
\alias{filter_mscore_fdr}
\title{
Filter annotated OpenSWATH/pyProphet output table to achieve a high FDR quality data matrix with controlled overall protein FDR and quantitative values for all peptides mapping to these high-confidence proteins (up to a desired overall peptide level FDR quality).
}
\description{
This function controls the protein FDR over a multi-run OpenSWATH/pyProphet output table and filters all quantitative values to a desired overall/global peptide FDR level.

It first finds a suitable m-score cutoff to minimally achieve a desired global FDR quality on a protein master list based on the function mscore4protfdr.
It then finds a suitable m-score cutoff to minimally achieve a desired global FDR quality on peptide level based on the function mscore4pepfdr.
Finally, it reports all the peptide quantities derived based on the peptide level cutoff for only those peptides mapping to the protein master list.
It further summarizes the protein and peptide numbers remaining after the filtering.
It further evaluates the individual run FDR qualities of the peptides (and quantitation events) selected.
}
\usage{
filter_mscore_fdr(data, FFT = 1, overall_protein_fdr_target = 0.02,
upper_overall_peptide_fdr_limit = 0.05, rm.decoy = TRUE, mscore.col = "m_score")
}

\arguments{
  \item{data}{
Annotated OpenSWATH/pyProphet data table
}
  \item{FFT}{
Ratio of false positives to true negatives, q-values from [Injection_name]_full_stat.csv in pyProphet stats output. As an approximation, the q-values of multiple runs are averaged and supplied as argument FFT. Numeric from 0 to 1. Defaults to 1, the most conservative value (1 Decoy indicates 1 False target). For further details see the Vignette Section 1.3 and 4.1.
}
  \item{overall_protein_fdr_target}{
FDR target for the protein master list for which quantitative values down to the less strict peptide_fdr criterion will be kept/reported. Defaults to 0.02.
}
  \item{upper_overall_peptide_fdr_limit}{
FDR target for the quantitative values kept/reported for all peptides mapping to the high-confidence protein master list. Defaults to 0.05. If all values up to m_score 0.01 shall be kept, set = 1.
}
  \item{rm.decoy}{
Logical T/F, whether decoy entries should be removed after the analysis. Defaults to TRUE. Can be useful to disable to track the influence on decoy fraction by further filtering steps such as requiring 2 peptides per protein.
}
  \item{mscore.col}{
Defines the column from which to retrieve the m_score. If you use JPP (Rosenberger, Bludau et al. 2017) this can be used to select between Protein and transition_group m_score.
}

}
\value{
\item{data.filtered}{the filtered data frame}
}

\author{
Moritz Heusel
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.fdr.filtered<-filter_mscore_fdr(data, FFT=0.7, overall_protein_fdr_target=0.02,
upper_overall_peptide_fdr_limit=0.1)
}
