#' Assess overall FDR in annotated OpenSWATH/pyProphet output table in
#' dependence of m_score cutoff
#'
#' This function estimates the assay, peptide and protein FDR over a multi-run
#' OpenSWATH/pyProphet output table. It counts target and decoy assays (unique
#' transition_group_id), peptides (unique FullPeptideName) and proteins (unique
#' ProteinName) in dependence of the m-score cutoff (1e-2 to 1e-20).
#' To arrive from decoy counts at an estimation of the false discovery rate
#' (false positives among the targets remaining at a given mscore cutoff) the
#' ratio of false positives to true negatives (decoys) (FFT) must be
#' supplied. It is estimated for each run individually by pyProphet and
#' contained in the pyProphet statistics [Injection_name]_full_stat.csv. As an
#' approximation, the FFTs of multiple runs are averaged and supplied as
#' argument FFT. For further details see the Vignette Section 1.3 and
#' 4.1. Protein FDR control on peak group quality level is a very strict filter
#' and should be handled with caution.
#' FDR is calculated as FDR = (TN*FFT/T); TN=decoys, T=targets, FFT=see above
#'
#' @param data  Data table that is produced by the OpenSWATH/pyProphet workflow
#' @param FFT Ratio of false positives to true negatives, q-values from
#'   [Injection_name]_full_stat.csv in pyProphet stats output. As an
#'   approximation, the q-values of multiple runs are averaged and supplied as
#'   argument FFT. Numeric from 0 to 1. Defaults to 1, the most conservative
#'   value (1 Decoy indicates 1 False target).
#' @param  n_range  I am also not certain what this is, nor why 20 is the
#'   optimal default value, but I think the idea is to set up a series of mscore
#'   thresholds.
#' @param output Choose output type. "pdf_csv" creates the output as files in
#'   the working directory, "Rconsole" triggers delivery of the output to the
#'   console enabling further computation or custom plotting / output.
#' @param plot  Logical, whether or not to create plots from the results (using
#'   the associated method plot.fdr_table()
#' @param filename  Optional, modifying the basename of the result files if
#'   applicable.
#' @return Returns a list of class "fdr_table". If output "pdf_csv" and plot =
#'   TRUE were chosen, report files are written to the working folder.
#' @author Moritz Heusel
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  assess_fdr_overall(data, FFT=0.7, output="Rconsole", plot=TRUE,
#'                     filename="Testoutput_assess_fdr_overall")
#' }
#' @export
assess_fdr_overall <-function(data, FFT=1, n_range=20, output="pdf_csv",
                              plot=TRUE, filename="FDR_report_overall") {
  mscore_levels <- 10 ^ - seq(n_range)
  ## create vectors to store count results
  target_assays <- NULL
  decoy_assays <- NULL
  target_peptides <- NULL
  decoy_peptides <- NULL
  target_proteins <- NULL
  decoy_proteins <- NULL

  # loop over IDs with different m_score thresholds counting targets & decoys
  for (i in seq_len(length(mscore_levels))) {
    target_assays[i] <- length(unique(data[data[["decoy"]] == FALSE &
                                           data[["m_score"]] <= mscore_levels[i],
                                           c("transition_group_id")]))
    decoy_assays[i] <- length(unique(data[data[["decoy"]] == TRUE &
                                          data[["m_score"]] <= mscore_levels[i],
                                          c("transition_group_id")]))
    target_peptides[i] <- length(unique(data[data[["decoy"]] == FALSE &
                                             data[["m_score"]] <= mscore_levels[i],
                                             c("fullpeptidename")]))
    decoy_peptides[i] <- length(unique(data[data[["decoy"]] == TRUE &
                                            data[["m_score"]] <= mscore_levels[i],
                                            c("fullpeptidename")]))
    target_proteins[i] <- length(unique(data[data[["decoy"]] == FALSE &
                                             data[["m_score"]] <= mscore_levels[i],
                                             c("proteinname")]))
    decoy_proteins[i] <- length(unique(data[data[["decoy"]] == TRUE &
                                            data[["m_score"]] <= mscore_levels[i],
                                            c("proteinname")]))
  }

  # calculate false target fraction at cutoff (FDR) by decoy fraction * FFT
  assay_fdr <- (decoy_assays / target_assays) * FFT
  peptide_fdr <- (decoy_peptides / target_peptides) * FFT
  protein_fdr <- (decoy_proteins / target_proteins) * FFT
  # calculate estimated number of true identifications among the target hits
  true_target_assays <- target_assays - (decoy_assays * FFT)
  true_target_peptides <- target_peptides - (decoy_peptides * FFT)
  true_target_proteins <- target_proteins - (decoy_proteins * FFT)

  fdr_table <- structure(list(
    "mscore_cutoff" = mscore_levels,
    "target_assays" = target_assays,
    "decoy_assays" = decoy_assays,
    "assay_fdr" = assay_fdr,
    "true_target_assays" = true_target_assays,
    "target_peptides" = target_peptides,
    "decoy_peptides" = decoy_peptides,
    "peptide_fdr" = peptide_fdr,
    "true_target_peptides" = true_target_peptides,
    "target_proteins" = target_proteins,
    "decoy_proteins" = decoy_proteins,
    "protein_fdr" = protein_fdr,
    "true_target_proteins" = true_target_proteins),
    class="fdr_table")

  if (isTRUE(plot)) {
    plotret <- plot_fdr_table(fdr_table, filename=filename, output=output)
  }

  if (output == "Rconsole") {
    return(fdr_table)
  }

  if (output == "pdf_csv") {
    fdr_table_csv <- cbind(mscore_levels, target_assays, decoy_assays,
                           assay_fdr, true_target_assays, target_peptides,
                           decoy_peptides, peptide_fdr, true_target_peptides,
                           target_proteins, decoy_proteins, protein_fdr, true_target_proteins)
    write.csv(fdr_table_csv, file=paste(filename, "table.csv", sep="_"), row.names=TRUE)
    message(paste0(filename, "_table.csv written to working folder"))
  }
  return(fdr_table)
}
