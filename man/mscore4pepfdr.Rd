\name{mscore4pepfdr}
\alias{mscore4pepfdr}
\title{
Find m_score cutoff to reach a desired FDR on peptide level (over the entire OpenSWATH/pyProphet output table)
}
\description{
This function estimates the m_score cutoff required in a dataset to reach a given overall peptide level FDR.
It counts target and decoy peptides (unique FullPeptideName) at high resolution across the m_score cutoffs and reports a useful m_score cutoff - peptide FDR pair close to the supplied fdr_target level over the entire dataset. The m_score cutoff is returned by the function and can be used in the context of the filtering functions, e.g.:

data.pepFDR2pc<-filter_mscore(data, mscore4pepfdr(data, fdr_target=0.02))

To arrive from decoy counts at an estimation of the false discovery rate (false positives among the targets remaining at a given mscore cutoff) the ratio of false positives to true negatives (decoys) (FFT) must be supplied. It is estimated for each run individually by pyProphet and contained in the pyProphet statistics [Injection_name]_full_stat.csv. As an approximation, the FFTs of multiple runs are averaged and supplied as argument FFT. For further details see the Vignette Section 1.3 and 4.1.

For FDR evaluations on assay and protein level, please refer to functions
mscore4assayfdr
mscore4protfdr
}
\usage{
mscore4pepfdr(data, FFT = 1, fdr_target = 0.01, mscore.col = "m_score")
}
\arguments{
  \item{data}{
Annotated OpenSWATH/pyProphet data table. See function sample_annotation from this package.
}
  \item{FFT}{
Ratio of false positives to true negatives, q-values from [Injection_name]_full_stat.csv in pyProphet stats output. As an approximation, the q-values of multiple runs are averaged and supplied as argument FFT. Numeric from 0 to 1. Defaults to 1, the most conservative value (1 Decoy indicates 1 False target).
}
  \item{fdr_target}{
FDR target, numeric, defaults to 0.01. An m_score cutoff achieving an FDR < fdr_target will be selected.
Calculated as FDR = (TN*FFT/T); TN=decoys, T=targets, FFT=see above.
}
  \item{mscore.col}{
Defines the column from which to retrieve the m_score. If you use JPP (Rosenberger, Bludau et al. 2017) this can be used to select between Protein and transition_group m_score.
}

}
\value{
Returns the m_score cutoff selected to arrive at the desired FDR
}


\author{
Moritz Heusel
}
\examples{
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
mscore4pepfdr(data, FFT=0.7, fdr_target=0.01)
}
