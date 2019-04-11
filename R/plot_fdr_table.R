#' Plot functionality for results of class "fdr_table" as produced by e.g. the function assess_fdr_overall()
#'
#' This function created standard plots from results of class "fdr_table" as
#' produced by e.g. the function assess_fdr_overall() visualizig ID numbers in
#' dependence of estimated FDR and also estimated FDR in dependence of m_score
#' cutoff.
#'
#' @param x  List of class "fdr_table" as produced e.g. by the function
#'   assess_fdr_overall() from this package.
#' @param output Choose output type. "pdf_csv" creates the output as files in
#'   the working directory, "Rconsole" triggers delivery of the output to the
#'   console enabling further computation or custom plotting / output.
#' @param filename  Basename for output files to be created (if output =
#'   "pdf_csv" has been selected).
#' @param ...  Extra arguments passed on to functions inside this.
#' @return  Originally this returned nothing, but now it makes a list of ggplot2
#'   plots which may be passed along and plotted as desired (with that in mind,
#'   I would like to remove the explicit plot() calls in this function).
#' @author Moritz Heusel
#' @examples
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  x <- assess_fdr_overall(data, FFT=0.7, output="Rconsole", plot=FALSE)
#'  plot.fdr_table(x, output="pdf_csv", filename="Assess_fdr_overall_testplot")
#' @export
plot_fdr_table <- function(x, output="Rconsole", filename="FDR_report_overall", ...) {
  ## Plot and create output from ID-FDR report (x) ##
  ## Plot 1: Target/true target curves
  ## as estimated by decoy counting + FFT-correction

  ## Save previous par(mfrow) settings and restore upon exit (BioConductor suggestion)
  par.mfrow.old <- par(mfrow=par()$mfrow)
  on.exit(par(par.mfrow.old), add=TRUE)
  par.mar.old <- par(mar=par()$mar)
  on.exit(par(par.mar.old), add=TRUE)

  ## Open pdf if output as pdf was desired by user
  if (output == "pdf_csv") {
    pdf(file=paste(filename, ".pdf", sep=""), height=6, width=10)
    par(mfrow=c(1,3), mar=c(13, 4, 15, 2) + 0.1)
  }

  ##plot 1.1 Assay-level sensitivity
  plot(x[["assay_fdr"]], x[["target_assays"]],
       xlab="assay FDR", ylab="# of assays",
       type = "l", lty=2,
       xlim=c(0, max(x[["assay_fdr"]], na.rm=TRUE)),
       ylim=(c(0, max(x[["target_assays"]]))))
  lines(x[["assay_fdr"]], x[["true_target_assays"]],
        xlim=c(0, max(x[["assay_fdr"]])))
  legend("topleft", legend=c("all targets", "true targets"), cex=0.5, lty=c(2, 1))

  ##plot 1.2 Peptide-level sensitivity
  plot(x[["peptide_fdr"]], x[["target_peptides"]],
       xlab="peptide FDR", ylab="# of peptides",
       type = "l", lty=2,
       xlim=c(0, max(x[["peptide_fdr"]], na.rm=TRUE)),
       ylim=(c(0, max(x[["target_peptides"]]))))
  lines(x[["peptide_fdr"]], x[["true_target_peptides"]],
        xlim=c(0, max(x[["peptide_fdr"]])))
  legend("topleft", legend=c("all targets", "true targets"), cex=0.5, lty=c(2, 1))

  if (output == "pdf_csv") {
    mtext("SWATH2stats global FDR & sensitivity report of OpenSWATH/pyProphet results", line=8, cex=1.5)
    mtext("Overall sensitivity:", line=3, adj=0.5, cex=1.1)
  }

  ## plot 1.3 Protein-level sensitivity
  plot(x[["protein_fdr"]], x[["target_proteins"]],
       xlab="protein FDR", ylab="# of proteins",
       type = "l", lty=2,
       xlim=c(0, max(x[["protein_fdr"]], na.rm=TRUE)),
       ylim=(c(0, max(x[["target_proteins"]]))))
  lines(x[["protein_fdr"]], x[["true_target_proteins"]],
        xlim=c(0, max(x[["protein_fdr"]])))
  legend("topleft", legend=c("all targets", "true targets"), cex=0.5, lty=c(2, 1))

  ## Plot 2: Global m_score adjustment and connectivity to global FDR quality levels
  ## as estimated by decoy counting + FFT-correction
  xlimit <- -1 * (length(x[["assay_fdr"]]) + 1)

  par(mar=c(5, 8, 4, 8) + 0.1, mfrow=c(1, 1))
  plot(log10(x[["mscore_cutoff"]]), x[["assay_fdr"]],
       axes=FALSE, ylim=c(0, 1.1 * max(x[["assay_fdr"]], na.rm=TRUE)),
       main="Global m-score cutoff connectivity to FDR quality",
       xlab="", ylab="", type="l",col="black", xlim=c(xlimit, 0))
  points(log10(x[["mscore_cutoff"]]), x[["assay_fdr"]], pch=20, col="black")
  axis(2, ylim=c(0, 1.1 * max(x[["assay_fdr"]], na.rm=TRUE)), col="black", lwd=2)
  mtext(2, text="Assay FDR", line=2)
  par(new=TRUE)
  plot(log10(x[["mscore_cutoff"]]), x[["peptide_fdr"]], axes=FALSE,
       ylim=c(0, 1.1 * max(x[["peptide_fdr"]], na.rm=TRUE)),
       xlab="", ylab="", type="l", lwd=2,
       col="red", main="", xlim=c(xlimit, 0), lty=2)
  points(log10(x[["mscore_cutoff"]]), x[["peptide_fdr"]], pch=20, col="red")
  axis(2, ylim=c(0, 1.1 * max(x[["peptide_fdr"]], na.rm=TRUE)),
       col="red", lwd=2, line=3.5)
  mtext(2, text="Peptide FDR", col="red", line=5.5)
  par(new=TRUE)
  plot(log10(x[["mscore_cutoff"]]), x[["protein_fdr"]], axes=FALSE,
       ylim=c(0, 1.1 * max(x[["protein_fdr"]], na.rm=TRUE)), col="blue",
       xlab="", ylab="", type="l", lty=3, main="", xlim=c(xlimit, 0), lwd=2)
  points(log10(x[["mscore_cutoff"]]), x[["protein_fdr"]], col="blue", pch=20)
  axis(4, ylim=c(0, 1.1 * max(x[["protein_fdr"]], na.rm=TRUE)), col="blue", lwd=2)
  mtext(4, text="Protein FDR", col="blue", line=2)
  axis(1, at=seq(min(log10(x[["mscore_cutoff"]])), max(log10(x[["mscore_cutoff"]]))))
  grid()
  abline(v=seq(min(log10(x[["mscore_cutoff"]])),
               max(log10(x[["mscore_cutoff"]]))), col="grey", lty=3)
  mtext(1, text="log10(m_score cutoff)", col="black", line=2)

  if (output == "pdf_csv") {
    dev.off()
    message(filename,".pdf written to working folder", "\n")
  }
}
