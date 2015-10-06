plot.fdr_table <- function(x, output="Rconsole", filename = "FDR_report_overall", ...){
## Plot and create output from ID-FDR report (x) ##
# Plot 1: Target/true target curves
# as estimated by decoy counting + FFT-correction

# Save previous par(mfrow) settings and restore upon exit (BioConductor suggestion)
par.mfrow.old <- par(mfrow= par()$mfrow)
on.exit(par(par.mfrow.old), add = TRUE)
par.mar.old <- par(mar= par()$mar)
on.exit(par(par.mar.old), add = TRUE)

# Open pdf if output as pdf was desired by user
if(output == "pdf_csv"){
  pdf(file=paste(filename, ".pdf", sep=""), height=6, width=10)
  par(mfrow=c(1,3), mar=c(13, 4, 15, 2) + 0.1)
}

#plot 1.1 Assay-level sensitivity
plot(x$assay.fdr, x$target.assays, xlab="assay FDR", ylab="# of assays", type = "l", 
     lty=2, xlim=c(0,max(x$assay.fdr, na.rm=TRUE)), ylim=(c(0,max(x$target.assays))))
lines(x$assay.fdr, x$true.target.assays, xlim=c(0,max(x$assay.fdr)))
legend("topleft", legend=c("all targets", "true targets"), cex=0.5, lty=c(2,1))

#plot 1.2 Peptide-level sensitivity
plot(x$peptide.fdr, x$target.peptides, xlab="peptide FDR", ylab="# of peptides", type = "l", 
     lty=2, xlim=c(0,max(x$peptide.fdr, na.rm=TRUE)), ylim=(c(0,max(x$target.peptides))))
lines(x$peptide.fdr, x$true.target.peptides, xlim=c(0,max(x$peptide.fdr)))
legend("topleft", legend=c("all targets", "true targets"), cex=0.5, lty=c(2,1))

if(output == "pdf_csv"){
  mtext("SWATH2stats global FDR & sensitivity report of OpenSWATH/pyProphet results", line=8, cex=1.5)
  mtext("Overall sensitivity:", line=3, adj=0.5, cex=1.1)
}

# plot 1.3 Protein-level sensitivity
plot(x$protein.fdr, x$target.proteins, xlab="protein FDR", ylab="# of proteins", type = "l",
     lty=2, xlim=c(0,max(x$protein.fdr, na.rm=TRUE)), ylim=(c(0,max(x$target.proteins))))
lines(x$protein.fdr, x$true.target.proteins, xlim=c(0,max(x$protein.fdr)))
legend("topleft", legend=c("all targets", "true targets"), cex=0.5, lty=c(2,1))

# Plot 2: Global m_score adjustment and connectivity to global FDR quality levels
# as estimated by decoy counting + FFT-correction
par(mar=c(5, 8, 4, 8) + 0.1, mfrow=c(1,1))

plot(log10(x$mscore_cutoff), x$assay.fdr, axes=FALSE, 
     ylim=c(0,1.1*max(x$assay.fdr, na.rm=TRUE)), 
     main="Global m-score cutoff connectivity to FDR quality", 
     xlab="", ylab="",type="l",col="black", xlim=c(-21,0))
points(log10(x$mscore_cutoff), x$assay.fdr, pch=20,col="black")
axis(2, ylim=c(0,1.1*max(x$assay.fdr, na.rm=TRUE)),col="black",lwd=2)
mtext(2,text="Assay FDR",line=2)
par(new=TRUE)
plot(log10(x$mscore_cutoff), x$peptide.fdr, axes=FALSE, 
     ylim=c(0,1.1*max(x$peptide.fdr, na.rm=TRUE)), 
     xlab="", ylab="",type="l", lwd=2, col="red", main="",xlim=c(-21,0), lty=2)
points(log10(x$mscore_cutoff), x$peptide.fdr, pch=20,col="red")
axis(2, ylim=c(0,1.1*max(x$peptide.fdr, na.rm=TRUE)),col="red",lwd=2, line =3.5)
mtext(2,text="Peptide FDR", col="red",line=5.5)
par(new=TRUE)
plot(log10(x$mscore_cutoff), x$protein.fdr, axes=FALSE, 
     ylim=c(0,1.1*max(x$protein.fdr, na.rm=TRUE)), col="blue",
     xlab="", ylab="", type="l",lty=3, main="",xlim=c(-21,0),lwd=2)
points(log10(x$mscore_cutoff), x$protein.fdr, col="blue",pch=20)
axis(4, ylim=c(0,1.1*max(x$protein.fdr, na.rm=TRUE)), col="blue",lwd=2)
mtext(4,text="Protein FDR", col="blue",line=2)
axis(1, at=seq(min(log10(x$mscore_cutoff)), max(log10(x$mscore_cutoff))))
grid()
abline(v=seq(min(log10(x$mscore_cutoff)), max(log10(x$mscore_cutoff))), col="grey", lty=3)
mtext(1,text="log10(m_score cutoff)", col="black", line =2)

if(output == "pdf_csv"){
  dev.off()
  message(filename,".pdf written to working folder", "\n")
}

}
