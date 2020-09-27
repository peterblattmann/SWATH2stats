#' S3 plot function for FDR assessment result arrays as produced by e.g. the
#' function assess_fdr_byrun()
#'
#' This function creates standard plots from result arrays as produced by
#' e.g. the function assess_fdr_byrun(), visualizig assay, peptide and protein
#' level FDR for each run at m-score cutoffs 1e-2 and 1e-3. Furthermore, Target
#' and Decoy ID numbers are visualized.
#'
#' @param x Array of by-run FDR assessment results as produced e.g. by the
#'   function assess_fdr_byrun() from this package.
#' @param output Choose output type. "pdf_csv" creates the output as files in
#'   the working directory, "Rconsole" triggers delivery of the output to the
#'   console enabling further computation and/or custom plotting / output.
#' @param filename  Basename for output files to be created (if output =
#'   "pdf_csv" has been selected).
#' @param plot_mscore_levels  Define m-score levels to plot the estimated FDR
#'   results.
#' @param ...  Extra arguments passed on to functions inside this.
#' @return  Plots in Rconsole or report files.
#' @author Moritz Heusel
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  x <- assess_fdr_byrun(data, FFT=0.7, output="Rconsole", plot=FALSE)
#'  retlist <- plot.fdr_cube(x, output="pdf_csv", filename="Assess_fdr_byrun_testplot",
#'                           plot_mscore_levels=0.01)
#' @export
plot.fdr_cube <- function(x, output = "Rconsole", filename = "FDR_report_byrun",
    plot_mscore_levels = c(0.01, 0.001), ...) {
    for (i in plot_mscore_levels) {
        k.mscore <- which(dimnames(x)[[3]] == i)
        k.mscore.label <- as.numeric(dimnames(x)[[3]][k.mscore])

        # write pdf reports
        run_ids <- unlist(dimnames(x)[2])

        # Save previous par(mfrow) settings and restore upon exit (BioConductor
        # suggestion)
        par.mfrow.old <- par(mfrow = par()$mfrow)
        on.exit(par(par.mfrow.old), add = TRUE)

        if (output == "pdf_csv") {
            pdf(file = paste0(filename, "_", format(k.mscore.label, scientific = TRUE),".pdf"))
            par(mfrow = c(3, 2))
            title <- NULL
        }

        if (output == "Rconsole") {
            title <- paste0("mscore_", format(k.mscore.label, scientific = TRUE))
        }

        # Assay level plots
        barplot(x[seq_len(2), , k.mscore], 
                ylim = c(0, 1.2 * max(x[seq_len(2), , k.mscore], na.rm = TRUE)), 
                ylab = "# of assays", legend = TRUE, 
                args.legend = list(x = "topleft", cex = 0.5), names.arg = run_ids, 
                cex.names = 0.5, las = 2, main = title)
        barplot(x[4, , k.mscore], 
                ylim = c(0, 1.2 * max(x[4, , k.mscore], na.rm = TRUE)),
                ylab = row.names(x)[4], names.arg = run_ids, 
                cex.names = 0.5, las = 2, main = title)
        # Peptide level plots
        barplot(x[5:6, , k.mscore], 
                ylim = c(0, 1.2 * max(x[5:6, , k.mscore], na.rm = TRUE)),
                ylab = "# of peptides", legend = TRUE, 
                args.legend = list(x = "topleft",cex = 0.5), 
                names.arg = run_ids, cex.names = 0.5, las = 2, main = title)
        barplot(x[8, , k.mscore], 
                ylim = c(0, 1.2 * max(x[8, , k.mscore], na.rm = TRUE)),
                ylab = row.names(x)[8], names.arg = run_ids, cex.names = 0.5, las = 2,
                main = title)
        # Protein level plots
        barplot(x[9:10, , k.mscore], 
                ylim = c(0, 2 * max(x[9:10, , k.mscore], na.rm = TRUE)),
                ylab = "# of proteins", legend = TRUE, 
                args.legend = list(x = "topleft", cex = 0.5), 
                names.arg = run_ids, cex.names = 0.5, las = 2, main = title)
        barplot(x[12, , k.mscore], 
                ylim = c(0, 1.2 * max(x[12, , k.mscore], na.rm = TRUE)),
                ylab = row.names(x)[12], 
                names.arg = run_ids, cex.names = 0.5, las = 2,
                main = title)

        if (output == "pdf_csv") {
            par(mfrow = c(1, 1))
            mtext(paste0("FDR by run target/decoy ID report (m_score <= ", 
                        format(k.mscore.label, scientific = TRUE), ")"), 
                  line = 2, adj = 0.5)
            dev.off()
            message(filename, ".pdf report plots written to working folder", "\n")
        }
    }
}
