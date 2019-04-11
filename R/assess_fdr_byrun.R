#' Assess assay, peptide and protein level FDR by run (for each MS_injection
#' separately) in OpenSWATH output table
#'
#' This function estimates the assay, peptide and protein FDR by run in an
#' OpenSWATH result table in dependence of a range of m_score cutoffs. The
#' results can be visualized and summarized by the associated method
#' plot.fdr_table().
#' It counts target and decoy assays (unique transition_group_id), peptides
#' (unique FullPeptideName) and proteins (unique ProteinName) in the OpenSWATH
#' output table in dependence of m-score cutoff, the useful m_score cutoff range
#' is evaluated for each dataset individually on the fly.
#' To arrive from decoy counts at an estimation of the false discovery rate
#' (false positives among the targets remaining at a given mscore cutoff) the
#' ratio of false positives to true negatives (decoys) (FFT) must be
#' supplied. It is estimated for each run individually by pyProphet and
#' contained in the pyProphet statistics [Injection_name]_full_stat.csv. As an
#' approximation, the FFTs of multiple runs are averaged and supplied as
#' argument FFT. For further details see the Vignette Section 1.3 and 4.1.
#' To assess fdr over the entire dataset, please refer to function
#' assess_fdr_overall. FDR is calculated as FDR = (TN*FFT/T); TN=decoys,
#' T=targets, FFT=see above.
#'
#' @param data Annotated OpenSWATH/pyProphet output table. Refer to function
#'   sample_annotation from this package for further information.
#' @param FFT Ratio of false positives to true negatives, q-values from
#'   [Injection_name]_full_stat.csv in pyProphet stats output. As an
#'   approximation, the q-values of multiple runs are averaged and supplied as
#'   argument FFT. Numeric from 0 to 1. Defaults to 1, the most conservative
#'   value (1 Decoy indicates 1 False target).
#' @param n_range  Option to set the number of magnitude for which the m_score
#'   threshold is decreased (e.g. n.range = 10, m-score from 0.1 until
#'   10^-10)^.
#' @param output Choose output type. "pdf_csv" creates the output as files in
#'   the working directory, "Rconsole" triggers delivery of the output to the
#'   console enabling further computation or custom plotting / output.
#' @param plot Logical, whether or not to create plots from the results (using
#'   the associated method plot.fdr_cube()
#' @param filename  Modify the basename of the result files if set.
#' @param output_mscore_levels Define m-score levels to plot and write the
#'   estimated FDR results.
#' @return Returns an array of target/decoy identification numbers and
#'   calculated FDR values at different m-score cutoffs.
#' @author Moritz Heusel
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  assessed <- assess_fdr_byrun(data, FFT=0.7, output="pdf_csv", plot=TRUE,
#'                               filename="Testoutput_assess_fdr_byrun")
#'  summary(assessed)
#' @export
assess_fdr_byrun <- function(data, FFT=1, n_range=20, output="pdf_csv", plot=TRUE,
                             filename="FDR_report_byrun", output_mscore_levels=c(0.01, 0.001)) {
  ## create m_score intervals to be tested
  test_levels <- 10 ^ -seq(1:n_range)

  ## Identify the minimal m-score cutoff at which all runs still contain decoys
  decoy_count_lengths <- NULL
  for (i in seq(1:n_range)) {
    decoy_count_lengths[i] <- length(by(data[data[["decoy"]] == TRUE &
                                             data[["m_score"]] <= test_levels[i],
                                             c("transition_group_id")],
                                        data[data[["decoy"]] == TRUE &
                                             data[["m_score"]] <= test_levels[i],
                                             c("run_id")],
                                        length))
  }

  mscore_limit <- length(
    decoy_count_lengths[decoy_count_lengths == length(unique(data[["run_id"]]))])

  if (mscore_limit < 2) {
    mscore_limit <- 2
  }

  mscore_levels <- 10 ^ -c(seq(2, mscore_limit))
  output_mscore_levels <- output_mscore_levels[output_mscore_levels %in% mscore_levels]

  fdr_cube <- array(NA, dim=c(12, length(unique(data[["run_id"]])), length(mscore_levels)))

  dimnames(fdr_cube) <- list(c("target_assays", "decoy_assays", "false_target_assays",
        "assay_fdr", "target_peptides", "decoy_peptides", "false_target_peptides", "peptide_fdr",
        "target_proteins", "decoy_proteins", "false_target_proteins", "protein_fdr"),
        sort(unique(data[["run_id"]])), mscore_levels)

  length_unique <- function(X) {
    length(unique(X))
  }

  data.t <- data[data[["decoy"]] == 0, ]
  data.d <- data[data[["decoy"]] == 1, ]

  ## for each m_score cutoff in mscore_levels, count targets and decoys and calculate false targets and fdr
  for (i in seq_len(length(mscore_levels))) {
    ## for each run_id, count target (and decoy) assays identified ("id" column entries)
    ## and store in pane i /row 1 (targets) /row 2 (decoys) & calculate false targets & fdr based on FFT
    fdr_cube[1, , i] <- by(data.t[data.t[["m_score"]] <= mscore_levels[i], c("transition_group_id")],
                           data.t[data.t[["m_score"]] <= mscore_levels[i], c("run_id")], length)
    fdr_cube[2, , i] <- by(data.d[data.d[["m_score"]] <= mscore_levels[i], c("transition_group_id")],
                           data.d[data.d[["m_score"]] <= mscore_levels[i], c("run_id")], length)

    fdr_cube[3, , i] <- fdr_cube[2, , i] * FFT
    fdr_cube[4, , i] <- fdr_cube[3, , i] / fdr_cube[1, , i]

    ## for each run_id, count target (and decoy) peptides identified (unique "FullPeptideName" column entries)
    ## and store in pane i /row 1 (targets) /row 2 (decoys) & calculate false targets & fdr based on FFT
    fdr_cube[5, , i] <- by(data.t[data.t[["m_score"]] <= mscore_levels[i],
                                  c("fullpeptidename")],
                           data.t[data.t[["m_score"]] <= mscore_levels[i],
                                  c("run_id")], length_unique)
    if (nrow(data.d) > 0) {
      fdr_cube[6, , i] <- by(data.d[data.d[["m_score"]] <= mscore_levels[i],
                                    c("fullpeptidename")],
                             data.d[data.d[["m_score"]] <= mscore_levels[i],
                                    c("run_id")], length_unique)
    } else {
      tmp_cube <- fdr_cube[5, , i]
      for (c in 1:length(tmp_cube)) {
        tmp_cube[c] <- 0
      }
      fdr_cube[6, , i] <- tmp_cube
    }
    fdr_cube[7, , i] <- fdr_cube[6, , i] * FFT
    fdr_cube[8, , i] <- fdr_cube[7, , i] / fdr_cube[5, , i]

    ## for each run_id, count target (and decoy) proteins identified (unique "ProteinName" column entries)
    ## and store in pane i / row 1 (targets) / row 2 (decoys) & calculate false targets & fdr based on FFT
    fdr_cube[9, , i] <- by(data.t[data.t[["m_score"]] <= mscore_levels[i],
                                  c("proteinname")],
                           data.t[data.t[["m_score"]] <= mscore_levels[i],
                                  c("run_id")], length_unique)
    if (nrow(data.d) > 0) {
      fdr_cube[10, , i] <- by(data.d[data.d[["m_score"]] <= mscore_levels[i],
                                     c("proteinname")],
                              data.d[data.d[["m_score"]] <= mscore_levels[i],
                                     c("run_id")], length_unique)
    } else {
      tmp_cube <- fdr_cube[9, , i]
      for (c in 1:length(tmp_cube)) {
        tmp_cube[c] <- 0
      }
      fdr_cube[10, , i] <- tmp_cube
    }
    fdr_cube[11, , i] <- fdr_cube[10, , i] * FFT
    fdr_cube[12, , i] <- fdr_cube[11, , i] / fdr_cube[9, , i]
  }

  ## print fdr values to console (m_score cutoff 1e-2)
  message("The average FDR by run on assay level is ", round(mean(fdr_cube[4, , 1], na.rm=TRUE), digits=3))
  message("The average FDR by run on peptide level is ", round(mean(fdr_cube[8, , 1], na.rm=TRUE), digits=3))
  message("The average FDR by run on protein level is ", round(mean(fdr_cube[12, , 1], na.rm=TRUE), digits=3))

  if (output == "pdf_csv") {
    message("Individual run FDR qualities can be retrieved from ", paste0(filename, ".csv"))
    ## Write csv reports for mscore 1e-2 and 1e-3

    for (i in output_mscore_levels) {
      k.mscore <- which(dimnames(fdr_cube)[[3]] == i)
      k.mscore.label <- as.numeric(dimnames(fdr_cube)[[3]][k.mscore])
      write.csv(fdr_cube[, , k.mscore], file=paste0(filename, "table_mscore_",
                                                    format(k.mscore.label, scientific=TRUE), ".csv"))
    }
    message(filename, ".csv reports written to working folder.")
  }

  ##fdr_cube2 <- fdr_cube
  ##class(fdr_cube2) <- "fdr_cube"
  class(fdr_cube) <- "fdr_cube"

  if (isTRUE(plot)) {
    plot_fdr_cube(fdr_cube, output=output, filename=filename,
                  plot_mscore_levels=output_mscore_levels)
  }

  return(fdr_cube)
}
