utils::globalVariables(c("par", "pdf", "barplot", "mtext", "dev.off"))

plot.fdr_cube <- function(x, output = "Rconsole", filename = "FDR_report_byrun", 
                          plot_mscore_levels = c(0.01, 0.001), ...){
  for(i in plot_mscore_levels){
    k.mscore <- which(dimnames(x)[[3]] == i)  
    k.mscore.label <- as.numeric(dimnames(x)[[3]][k.mscore])
    
    # write pdf reports
    run_ids<-unlist(dimnames(x)[2])
    
    # Save previous par(mfrow) settings and restore upon exit (BioConductor suggestion)
    par.mfrow.old <- par(mfrow= par()$mfrow)
    on.exit(par(par.mfrow.old), add = TRUE)
    
    if(output == "pdf_csv"){
      pdf(file = paste(filename,  "_", format(k.mscore.label, scientific= TRUE), ".pdf", sep = ""))
      par(mfrow=c(3,2))
      title <- NULL
    }
    
    if(output == "Rconsole"){
      title <- paste("mscore_", format(k.mscore.label, scientific= TRUE), sep="")
    }
    
    # Assay level plots
    barplot(x[seq_len(2), , k.mscore], ylim=c(0,1.2*max(x[seq_len(2), , k.mscore], na.rm=TRUE)), ylab = "# of assays",
            legend = TRUE, args.legend = list(x="topleft", cex=0.5),  names.arg=run_ids, cex.names=0.5, las=2,
            main = title)
    barplot(x[4, , k.mscore], ylim=c(0,1.2*max(x[4, , k.mscore], na.rm=TRUE)) , ylab = row.names(x)[4],
            names.arg=run_ids, cex.names=0.5, las=2, main = title)
    # Peptide level plots
    barplot(x[5:6, , k.mscore], ylim=c(0,1.2*max(x[5:6, , k.mscore], na.rm=TRUE)), ylab = "# of peptides",
            legend = TRUE, args.legend = list(x="topleft", cex=0.5),  names.arg=run_ids, cex.names=0.5, las=2
            , main = title)
    barplot(x[8, , k.mscore], ylim=c(0,1.2*max(x[8, , k.mscore], na.rm=TRUE)) , ylab = row.names(x)[8],
            names.arg=run_ids, cex.names=0.5, las=2, main = title)
    # Protein level plots
    barplot(x[9:10, , k.mscore], ylim=c(0,2*max(x[9:10, , k.mscore], na.rm=TRUE)), ylab = "# of proteins",
            legend = TRUE, args.legend = list(x="topleft", cex=0.5),  names.arg=run_ids, cex.names=0.5, las=2
            , main = title)
    barplot(x[12, , k.mscore], ylim=c(0,1.2*max(x[12, , k.mscore], na.rm=TRUE)) , ylab = row.names(x)[12],
            names.arg=run_ids, cex.names=0.5, las=2, main = title)
    
    if(output == "pdf_csv"){
      par(mfrow=c(1,1))
      mtext(paste("FDR by run target/decoy ID report (m_score <= ", format(k.mscore.label, scientific= TRUE), ")", sep=""), line=2, adj=0.5)
      dev.off()
      message(filename, ".pdf report plots written to working folder", "\n")
    }
  }
}