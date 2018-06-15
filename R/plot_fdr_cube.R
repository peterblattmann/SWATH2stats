#' Plot functionality for FDR assessment result arrays as produced by e.g. the
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
#' @return  Originally this returned nothing, but now it makes a list of ggplot2
#'   plots which may be passed along and plotted as desired (with that in mind,
#'   I would like to remove the explicit plot() calls in this function).
#' @author Moritz Heusel
#' @examples
#' \dontrun{
#'  data("OpenSWATH_data", package="SWATH2stats")
#'  data("Study_design", package="SWATH2stats")
#'  data <- sample_annotation(OpenSWATH_data, Study_design)
#'  x <- assess_fdr_byrun(data, FFT=0.7, output = "Rconsole", plot = FALSE)
#'  plot.fdr_cube(x, output = "pdf_csv", filename = "Assess_fdr_byrun_testplot",
#'                plot_mscore_levels=0.01)
#' }
#' @export
plot_fdr_cube <- function(x, output="Rconsole", filename="FDR_report_byrun",
                          plot_mscore_levels=c(0.01, 0.001), ...) {
  retlist <- list()
  count <- 0
  for (i in plot_mscore_levels) {
    count <- count + 1
    retlist[[count]] <- list()
    k.mscore <- which(dimnames(x)[[3]] == i)
    k.mscore.label <- as.numeric(dimnames(x)[[3]][k.mscore])
    ## write pdf reports
    run_ids <- unlist(dimnames(x)[2])
    ## Save previous par(mfrow) settings and restore upon exit (BioConductor suggestion)
    par.mfrow.old <- par(mfrow=par()$mfrow)
    on.exit(par(par.mfrow.old), add=TRUE)

    if (output == "pdf_csv") {
      pdf(file=paste(filename,  "_", format(k.mscore.label, scientific=TRUE), ".pdf", sep=""))
      par(mfrow=c(3,2))
      title <- NULL
    }

    if (output == "Rconsole") {
      title <- paste("mscore_", format(k.mscore.label, scientific=TRUE), sep="")
    }

    ## Assay level plots
    ##barplot(x[1:2, , k.mscore], ylim=c(0,1.2*max(x[1:2, , k.mscore], na.rm=TRUE)), ylab = "# of assays",
    ##        legend = TRUE, args.legend = list(x="topleft", cex=0.5),  names.arg=run_ids, cex.names=0.5, las=2,
    ##        main = title)
    assays <- x[1:2, , k.mscore]
    assays_df <- reshape2::melt(assays)
    assay_barplot <- ggplot2::ggplot(data=assays_df,
                                  ggplot2::aes_string(x="Var2", y="value")) +
      ggplot2::geom_col(position="identity", color="black") +
      ggplot2::xlab("Sample ID") +
      ggplot2::ylab("# of assays") +
      ggplot2::theme_bw(base_size=base_size) +
      ggplot2::geom_text(parse=FALSE, angle=90, size=4, color="white", hjust=1.2,
                         ggplot2::aes_string(label='prettyNum(as.character(assays_df$value), big.mark=",")')) +
      ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                     axis.text.x=ggplot2::element_text(angle=90, vjust=0.5)) ##, hjust=1.5, vjust=0.5))
    plot(assay_barplot)

    ##barplot(x[4, , k.mscore], ylim=c(0,1.2*max(x[4, , k.mscore], na.rm=TRUE)) , ylab = row.names(x)[4],
    ##        names.arg=run_ids, cex.names=0.5, las=2, main = title)
    assay_mscores <- x[4, , k.mscore]
    assay_mscores_df <- reshape2::melt(assay_mscores)
    assay_mscores_df[["names"]] <- rownames(assay_mscores_df)
    assay_mscores_df[["value"]] <- signif(assay_mscores_df[["value"]], digits=4)
    assay_mscoreplot <- ggplot2::ggplot(data=assay_mscores_df,
                                        ggplot2::aes_string(x="names", y="value")) +
      ggplot2::geom_col(position="identity", color="black") +
      ggplot2::xlab("Sample ID") +
      ggplot2::ylab("assay fdr") +
      ggplot2::theme_bw(base_size=base_size) +
      ggplot2::geom_text(parse=FALSE, angle=90, size=4, color="white", hjust=1.2,
                         ggplot2::aes_string(label='prettyNum(as.character(assay_mscores_df$value), big.mark=",")')) +
      ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                     axis.text.x=ggplot2::element_text(angle=90, vjust=0.5)) ##, hjust=1.5, vjust=0.5))
    plot(assay_mscoreplot)

    ## Peptide level plots
    ##barplot(x[5:6, , k.mscore], ylim=c(0,1.2*max(x[5:6, , k.mscore], na.rm=TRUE)), ylab = "# of peptides",
    ##        legend = TRUE, args.legend = list(x="topleft", cex=0.5),  names.arg=run_ids, cex.names=0.5, las=2
    ##        , main = title)
    peptides <- x[5:6, , k.mscore]
    peptides_df <- reshape2::melt(peptides)
    peptide_barplot <- ggplot2::ggplot(data=peptides_df,
                                  ggplot2::aes_string(x="Var2", y="value")) +
      ggplot2::geom_col(position="identity", color="black") +
      ggplot2::xlab("Sample ID") +
      ggplot2::ylab("# of peptides") +
      ggplot2::theme_bw(base_size=base_size) +
      ggplot2::geom_text(parse=FALSE, angle=90, size=4, color="white", hjust=1.2,
                         ggplot2::aes_string(label='prettyNum(as.character(peptides_df$value), big.mark=",")')) +
      ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                     axis.text.x=ggplot2::element_text(angle=90, vjust=0.5)) ##, hjust=1.5, vjust=0.5))
    plot(peptide_barplot)

    ##barplot(x[8, , k.mscore], ylim=c(0,1.2*max(x[8, , k.mscore], na.rm=TRUE)) , ylab = row.names(x)[8],
    ##        names.arg=run_ids, cex.names=0.5, las=2, main = title)
    peptide_mscores <- x[8, , k.mscore]
    peptide_mscores_df <- reshape2::melt(peptide_mscores)
    peptide_mscores_df[["names"]] <- rownames(peptide_mscores_df)
    peptide_mscores_df[["value"]] <- signif(peptide_mscores_df[["value"]], digits=4)
    peptide_mscoreplot <- ggplot2::ggplot(data=peptide_mscores_df,
                                          ggplot2::aes_string(x="names", y="value")) +
      ggplot2::geom_col(position="identity", color="black") +
      ggplot2::xlab("Sample ID") +
      ggplot2::ylab("peptide fdr") +
      ggplot2::theme_bw(base_size=base_size) +
      ggplot2::geom_text(parse=FALSE, angle=90, size=4, color="white", hjust=1.2,
                         ggplot2::aes_string(label='prettyNum(as.character(peptide_mscores_df$value), big.mark=",")')) +
      ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                     axis.text.x=ggplot2::element_text(angle=90, vjust=0.5)) ##, hjust=1.5, vjust=0.5))
    plot(peptide_mscoreplot)

    ## Protein level plots
    ##barplot(x[9:10, , k.mscore], ylim=c(0,2*max(x[9:10, , k.mscore], na.rm=TRUE)), ylab = "# of proteins",
    ##        legend = TRUE, args.legend = list(x="topleft", cex=0.5),  names.arg=run_ids, cex.names=0.5, las=2
    ##      , main = title)
    proteins <- x[9:10, , k.mscore]
    proteins_df <- reshape2::melt(proteins)
    protein_barplot <- ggplot2::ggplot(data=proteins_df,
                                       ggplot2::aes_string(x="Var2", y="value")) +
      ggplot2::geom_col(position="identity", color="black") +
      ggplot2::xlab("Sample ID") +
      ggplot2::ylab("# of proteins") +
      ggplot2::theme_bw(base_size=base_size) +
      ggplot2::geom_text(parse=FALSE, angle=90, size=4, color="white", hjust=1.2,
                         ggplot2::aes_string(label='prettyNum(as.character(proteins_df$value), big.mark=",")')) +
      ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                     axis.text.x=ggplot2::element_text(angle=90, vjust=0.5)) ##, hjust=1.5, vjust=0.5))
    plot(protein_barplot)

    ##barplot(x[12, , k.mscore], ylim=c(0,1.2*max(x[12, , k.mscore], na.rm=TRUE)) , ylab = row.names(x)[12],
    ##        names.arg=run_ids, cex.names=0.5, las=2, main = title)
    protein_mscores <- x[12, , k.mscore]
    protein_mscores_df <- reshape2::melt(protein_mscores)
    protein_mscores_df[["names"]] <- rownames(protein_mscores_df)
    protein_mscores_df[["value"]] <- signif(protein_mscores_df[["value"]], digits=4)
    protein_mscoreplot <- ggplot2::ggplot(data=protein_mscores_df,
                                         ggplot2::aes_string(x="names", y="value")) +
      ggplot2::geom_col(position="identity", color="black") +
      ggplot2::xlab("Sample ID") +
      ggplot2::ylab("protein fdr") +
      ggplot2::theme_bw(base_size=base_size) +
      ggplot2::geom_text(parse=FALSE, angle=90, size=4, color="white", hjust=1.2,
                         ggplot2::aes_string(label='prettyNum(as.character(protein_mscores_df$value), big.mark=",")')) +
      ggplot2::theme(axis.text=ggplot2::element_text(size=base_size, colour="black"),
                     axis.text.x=ggplot2::element_text(angle=90, vjust=0.5)) ##, hjust=1.5, vjust=0.5))
    plot(protein_mscoreplot)

    retlist[[count]][["assay_barplot"]] <- assay_barplot
    retlist[[count]][["assay_mscore"]] <- assay_mscoreplot
    retlist[[count]][["peptide_barplot"]] <- peptide_barplot
    retlist[[count]][["peptide_mscore"]] <- peptide_mscoreplot
    retlist[[count]][["protein_barplot"]] <- protein_barplot
    retlist[[count]][["protein_mscore"]] <- protein_mscoreplot

    if (output == "pdf_csv") {
      par(mfrow=c(1, 1))
      mtext(paste("FDR by run target/decoy ID report (m_score <= ",
                  format(k.mscore.label, scientific= TRUE), ")", sep=""), line=2, adj=0.5)
      dev.off()
      message(filename, ".pdf report plots written to working folder", "\n")
    }
  }
  return(retlist)
}
