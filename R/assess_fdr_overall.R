utils::globalVariables(c("write.csv"))

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
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  assess_fdr_overall(data, FFT=0.7, output="Rconsole", plot=TRUE,
#'                     filename="Testoutput_assess_fdr_overall")
#' @export
assess_fdr_overall <- function(data, FFT = 1, n.range = 20, output = "pdf_csv", plot = TRUE,
    filename = "FDR_report_overall", mscore.col = "m_score") {
    mscore.col <- JPP_update(data, mscore.col)

    mscore_levels <- 10^-seq(n.range)

    # create vectors to store count results
    target.assays <- NULL
    decoy.assays <- NULL
    target.peptides <- NULL
    decoy.peptides <- NULL
    target.proteins <- NULL
    decoy.proteins <- NULL

    # loop over IDs with different m_score thresholds counting targets & decoys
    for (i in seq_len(length(mscore_levels))) {
        target.assays[i] <- length(unique(data[data$decoy == FALSE & data[, mscore.col] <=
            mscore_levels[i], c("transition_group_id")]))
        decoy.assays[i] <- length(unique(data[data$decoy == TRUE & data[, mscore.col] <=
            mscore_levels[i], c("transition_group_id")]))
        target.peptides[i] <- length(unique(data[data$decoy == FALSE & data[, mscore.col] <=
            mscore_levels[i], c("FullPeptideName")]))
        decoy.peptides[i] <- length(unique(data[data$decoy == TRUE & data[, mscore.col] <=
            mscore_levels[i], c("FullPeptideName")]))
        target.proteins[i] <- length(unique(data[data$decoy == FALSE & data[, mscore.col] <=
            mscore_levels[i], c("ProteinName")]))
        decoy.proteins[i] <- length(unique(data[data$decoy == TRUE & data[, mscore.col] <=
            mscore_levels[i], c("ProteinName")]))
    }

    # calculate false target fraction at cutoff (FDR) by decoy fraction * FFT
    assay.fdr <- (decoy.assays/target.assays) * FFT
    peptide.fdr <- (decoy.peptides/target.peptides) * FFT
    protein.fdr <- (decoy.proteins/target.proteins) * FFT
    # calculate estimated number of true identifications among the target hits
    true.target.assays <- target.assays - (decoy.assays * FFT)
    true.target.peptides <- target.peptides - (decoy.peptides * FFT)
    true.target.proteins <- target.proteins - (decoy.proteins * FFT)


    fdr_table <- structure(list(mscore_cutoff = mscore_levels, target.assays = target.assays,
        decoy.assays = decoy.assays, assay.fdr = assay.fdr, true.target.assays = true.target.assays,
        target.peptides = target.peptides, decoy.peptides = decoy.peptides, peptide.fdr = peptide.fdr,
        true.target.peptides = true.target.peptides, target.proteins = target.proteins,
        decoy.proteins = decoy.proteins, protein.fdr = protein.fdr, true.target.proteins = true.target.proteins),
        class = "fdr_table")

    if (isTRUE(plot)) {
        plot.fdr_table(fdr_table, filename = filename, output = output)
    }

    if (output == "Rconsole") {
        return(fdr_table)
    }

    if (output == "pdf_csv") {
        fdr_table_csv <- cbind(mscore_levels, target.assays, decoy.assays, assay.fdr,
            true.target.assays, target.peptides, decoy.peptides, peptide.fdr, true.target.peptides,
            target.proteins, decoy.proteins, protein.fdr, true.target.proteins)
        write.csv(fdr_table_csv, file = paste(filename, "table.csv", sep = "_"),
            row.names = TRUE)
        message(filename, "table.csv written to working folder", "\n")
    }

}
