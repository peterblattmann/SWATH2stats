plot.fdr_cube <- function(x, output = "Rconsole", filename = "FDR_report_byrun", ...){
  # write pdf reports
  run_ids<-unlist(dimnames(x)[2])
  
  # Save previous par(mfrow) settings and restore upon exit (BioConductor suggestion)
  par.mfrow.old <- par(mfrow= par()$mfrow)
  on.exit(par(par.mfrow.old), add = TRUE)
  
  # at m_score <= 1e-2
  if(output == "pdf_csv"){
    pdf(file = paste(filename, "mscore_1e-2.pdf", sep = "_"))
    par(mfrow=c(3,2))
    title <- NULL
  }
  
  if(output == "Rconsole"){
    title <- "mscore_1e-2"
  }
  
  
  # Assay level plots
  barplot(x[1:2, , 1], ylim=c(0,1.2*max(x[1:2, , 1])), ylab = "# of assays",
          legend = TRUE, args.legend = list(x="topleft", cex=0.5),  names.arg=run_ids, cex.names=0.5, las=2,
          main = title)
  barplot(x[4, , 1], ylim=c(0,1.2*max(x[4, , 1])) , ylab = row.names(x)[4],
          names.arg=run_ids, cex.names=0.5, las=2, main = title)
  # Peptide level plots
  barplot(x[5:6, , 1], ylim=c(0,1.2*max(x[5:6, , 1])), ylab = "# of peptides",
          legend = TRUE, args.legend = list(x="topleft", cex=0.5),  names.arg=run_ids, cex.names=0.5, las=2
          , main = title)
  barplot(x[8, , 1], ylim=c(0,1.2*max(x[8, , 1])) , ylab = row.names(x)[8],
          names.arg=run_ids, cex.names=0.5, las=2, main = title)
  # Protein level plots
  barplot(x[9:10, , 1], ylim=c(0,2*max(x[9:10, , 1])), ylab = "# of proteins",
          legend = TRUE, args.legend = list(x="topleft", cex=0.5),  names.arg=run_ids, cex.names=0.5, las=2
          , main = title)
  barplot(x[12, , 1], ylim=c(0,1.2*max(x[12, , 1])) , ylab = row.names(x)[12],
          names.arg=run_ids, cex.names=0.5, las=2, main = title)
  
  if(output == "pdf_csv"){
    par(mfrow=c(1,1))
    mtext("FDR by run target/decoy ID report (m_score <= 1e-2)", line=2, adj=0.5)
    dev.off()
    
    # at m_score <= 1e-3
    pdf(file = paste(filename, "mscore_1e-3.pdf", sep = "_"))
    par(mfrow=c(3,2))
    title <- NULL
  }
  
  if(output == "Rconsole"){
    title <- "mscore_1e-3"
  }
  
  
  # Assay level plots
  barplot(x[1:2, , 2], ylim=c(0,1.2*max(x[1:2, , 2])), ylab = "# of assays",
          legend = TRUE, args.legend = list(x="topleft", cex=0.5),  names.arg=run_ids, cex.names=0.5, las=2
          , main = title)
  barplot(x[4, , 2], ylim=c(0,1.2*max(x[4, , 2])) , ylab = row.names(x)[4],
          names.arg=run_ids, cex.names=0.5, las=2, main = title)
  # Peptide level plots
  barplot(x[5:6, , 2], ylim=c(0,1.2*max(x[5:6, , 2])), ylab = "# of peptides",
          legend = TRUE, args.legend = list(x="topleft", cex=0.5),  names.arg=run_ids, cex.names=0.5, las=2
          , main = title)
  barplot(x[8, , 2], ylim=c(0,1.2*max(x[8, , 2])) , ylab = row.names(x)[8],
          names.arg=run_ids, cex.names=0.5, las=2, main = title)
  # Protein level plots
  barplot(x[9:10, , 2], ylim=c(0,2*max(x[9:10, , 2])), ylab = "# of proteins",
          legend = TRUE, args.legend = list(x="topleft", cex=0.5),  names.arg=run_ids, cex.names=0.5, las=2
          , main = title)
  barplot(x[12, , 2], ylim=c(0,1.2*max(x[12, , 2])) , ylab = row.names(x)[12],
          names.arg=run_ids, cex.names=0.5, las=2, main = title)
  
  if(output == "pdf_csv"){
    par(mfrow=c(1,1))
    mtext("FDR by run target/decoy ID report (m_score <= 1e-3)", line=2, adj=0.5)
    dev.off()
  
    message(filename, ".pdf report plots written to working folder", "\n")
  }
}