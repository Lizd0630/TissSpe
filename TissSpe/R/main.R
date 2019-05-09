# Function to detect tissue specific gene or AS events
# df: data.frame without
# names(meths) <- c("Tau", "Gini", "Tsi", "Counts", "Ee", "Pem", "Hg", "Z", "Spm")
# binary: "seq", "quant", "none"

ts_sp <- function(df,
                  method = "Tau",
                  binary = "seq",
                  n = 10,
                  min = 0.2,
                  max = 2,
                  norm_meth = "",
                  cutoff = 1,
                  name_vect,
                  identifier) {
  df <- na.omit(df)
  if(norm_meth == "log_QN") {
    dt <- df[, c(-1)]
    dt[dt < cutoff] <- 1
    dt <- log2(dt)
    cutoff <- log2(cutoff)
    df[, c(-1)] <- quant_norm(dt)
  } else if (norm_meth == "QN") {
    dt <- df[, c(-1)]
    dt[dt < cutoff] <- 0
    df[, c(-1)] <- quant_norm(dt)
  } else if (norm_meth == "log") {
    dt <- df[, c(-1)]
    dt[dt < cutoff] <- 1
    df[, c(-1)] <- log2(dt)
    cutoff <- log2(cutoff)
  } else {
    dt <- df[, c(-1)]
    dt[dt < cutoff] <- 0
    df[, c(-1)] <- dt
  }
  df <- rep_mean(df, name_vect, identifier)
  df$Max <- apply(df[,c(-1)], c(1), expr_max)
  df <- df[df$Max > cutoff,]
  df <- df[, c(-length(colnames(df)))]



  df_bin <- switch (binary,
              seq = seq_cut(df = df, n = n),
              quant = quant_cut(df = df, n = n, min = min, max = max),
              none = df
            )

  switch (meth,
    Tau = ts_Tau(df),
    Gini = ts_Gini(df),
    Tsi = ts_Tsi(df),
    Counts = ts_Counts(df),
    Ee = ts_Ee(df),
    Pem = ts_Pem(df),
    Hg = ts_Hg(df),
    Z = ts_Z(df),
    Spm = ts_Spm(df)
    )
}





###+++###
#Calculate and save tissue specificity parameters
#orgPSI = data set, PSI = cutt off, add = number of tissues, tNames = tissues to use, tNamesNew = tissues to name, RNAseq = how to normalise (log_QN, QN, log, NA)
#Only genes with Ensembl IDs are used, or for Drosophila
#Normalization is done on all tissues, not dependent which tissues are selected later
#1. Data are normalized
#2. All expression under PSI is set to 0
#3. Replicates mean is calculated (fReplicateMean)
#4. Genes that not expressed in any tissue are removed
#5. Tissue specificity parameters are calculated
fTS <- function(orgPSI, PSI, add, tNames, tNamesNew, RNAseq)
{
  orgPSI <- na.omit(orgPSI)
  print(summary(orgPSI))
  if(RNAseq == "log_QN"){
    x <- orgPSI[,c(-1)]
    x[x < PSI] <- 1
    x <- log2(x)
    PSI <- log2(PSI)
    orgPSI[,c(-1)] <- fQN(x)
  } else if (RNAseq == "QN")  {
    x <- orgPSI[,c(-1)]
    x[x < PSI] <- 0
    orgPSI[,c(-1)] <- fQN(x)
  } else if (RNAseq == "log")  {
    x <- orgPSI[,c(-1)]
    x[x < PSI] <- 1
    orgPSI[,c(-1)] <- log2(x)
    PSI <- log2(PSI)
  } else {
    x <- orgPSI[,c(-1)]
    x[x < PSI] <- 0
    orgPSI[,c(-1)] <- x
  }
  orgPSI <- fReplicateMean(orgPSI, organism, paste("Averaged.PSI.",tissuesNames, sep=""))
  orgPSI$Max <- apply(orgPSI[,c(-1)], c(1), fmax)
  orgPSI <- orgPSI[orgPSI$Max > PSI,]
  orgPSI <- orgPSI[,c(-length(colnames(orgPSI)))]
  print(summary(orgPSI))
  fPlotExpression(orgPSI[,-1], paste("Normalized expression (cutoff", 2^PSI, "PSI)", sep=" "), paste("NormalizedQN_", 2^PSI,"PSI", sep=""), tissuesPrintNames)

  orgPSI <- orgPSI[,c("SpliceJunction", paste("Averaged.PSI.", tNames,sep="")) ]
  colnames(orgPSI) <- c("SpliceJunction", paste("Averaged.PSI.", tNamesNew, sep=""))
  nTissues <- length(tNamesNew)
  tissuesNames <- tNamesNew
  print(paste("Analysis done on", nTissues, "tissue:", sep=" "))
  print(tissuesNames)

  orgPSI$Tau <- apply(orgPSI[,c(paste("Averaged.PSI.", tissuesNames[1:nTissues], sep=""))], 1, fTau)
  orgPSI$Gini <- apply(orgPSI[,c(paste("Averaged.PSI.", tissuesNames[1:nTissues], sep=""))], 1, fGini)
  orgPSI$Tsi <- apply(orgPSI[,c(paste("Averaged.PSI.", tissuesNames[1:nTissues], sep=""))], 1, fTsi)
  orgPSI$Counts <- apply(orgPSI[,c(paste("Averaged.PSI.", tissuesNames[1:nTissues], sep=""))], 1, function(x){x <- fCounts(x, PSI)})
  orgPSI$Hg <- apply(orgPSI[,c(paste("Averaged.PSI.", tissuesNames[1:nTissues], sep=""))], 1, fHg)
  orgPSI$Zscore <- fZ(orgPSI[,c(paste("Averaged.PSI.", tissuesNames[1:nTissues], sep=""))])
  orgPSI$Spm <- apply(orgPSI[,c(paste("Averaged.PSI.", tissuesNames[1:nTissues], sep=""))], 1, fSpm)
  orgPSI$Ee <- fEe(orgPSI[,c(paste("Averaged.PSI.", tissuesNames[1:nTissues], sep=""))])
  orgPSI$Pem <- fPem(orgPSI[,c(paste("Averaged.PSI.", tissuesNames[1:nTissues], sep=""))])

  orgPSI$Mean <- apply(orgPSI[,c(paste("Averaged.PSI.", tissuesNames[1:nTissues], sep=""))], 1, fmean)
  orgPSI$Max <- apply(orgPSI[,c(paste("Averaged.PSI.", tissuesNames[1:nTissues], sep=""))], 1, fmax)

  print(summary(orgPSI))

  p <- c("Tau", "Gini", "Tsi", "Counts", "Ee", "Hg", "Zscore", "Spm", "Pem")
  x <- as.matrix(orgPSI[,p])
  xs <- cor(x, method="spearman")
  xp <- cor(x, method="pearson")
  capture.output(c("Spearman correlation"),file=paste(folder, organism,"CorrelationTS_", add, ".txt", sep=""))
  capture.output(xs, append=TRUE, file=paste(folder, organism, "CorrelationTS_", add, ".txt", sep=""))
  capture.output(c("Pearson correlation"), append=TRUE, file=paste(folder, organism,"CorrelationTS_", add, ".txt", sep=""))
  capture.output(xp, append=TRUE, file=paste(folder, organism,"CorrelationTS_", add, ".txt", sep=""))

  # dev.new(height=9, width=12)
  #pdf(file=paste(folder, organism, expDataSource, "TScomparison_9_",  add,".pdf", sep=""), height=9, width=12)
  par(cex.main=0.95, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
  palette(rev(rich.colors(10)))
  #palette(rev(blues9))

  plot(density(orgPSI[,"Tau"],n=1000), main = " ", xlab="Tissue specificity",col=(1), lwd=4, lty=1
       ,ylim=c(0,8), xlim=c(-0.1,1.1)
  )
  lines(density(orgPSI[,"Gini"],n = 1000), col=(2), lwd=4, lty=2)
  lines(density(orgPSI[,"Tsi"],n = 1000), col=(3), lwd=4, lty=1)
  lines(density(orgPSI[,"Counts"],n = 1000), col=(4), lwd=4, lty=2)
  lines(density(orgPSI[,"Ee"],n = 1000), col=(5), lwd=4, lty=1)
  lines(density(orgPSI[,"Hg"],n = 1000), col=(6), lwd=4, lty=2)
  lines(density(orgPSI[,"Zscore"],n = 1000), col=(7), lwd=4, lty=1)
  lines(density(orgPSI[,"Spm"],n = 1000), col=(8), lwd=4, lty=2)
  lines(density(orgPSI[,"Pem"],n = 1000), col=(9), lwd=4, lty=1)

  legend("topright",c("Tau", "Gini", "TSI", "Counts", "EE", "Hg", "Zscore", "SPM", "PEM"),col=(1:11), lwd=4, lty=c(1,2), bty="n", seg.len=4)

  # dev.copy2pdf(device=quartz, file=paste(folder, organism, expDataSource, "TScomparison_PSI_",  add,".pdf", sep=""),onefile=TRUE)#,paper="A4r"
  #dev.off()

  write.table(orgPSI, file=paste(folder, organism,"TScomparisonTable_PSI_",  add,".txt",sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)

  fScatPlot(orgPSI, "Tau", c("Gini", "Tsi", "Counts", "Hg", "Zscore","Spm", "Ee", "Pem"), add, c(0, 0.94, 4, 0.5))
  fScatPlot2(orgPSI, "Mean", c("Tau", "Gini", "Tsi", "Counts", "Hg", "Zscore", "Spm", "Ee", "Pem"), add, c(ceiling(max(orgPSI$Mean)), 0.95, 2, 0.5), ceiling(max(orgPSI$Mean)))
  fScatPlot2(orgPSI, "Max", c("Tau", "Gini", "Tsi", "Counts", "Hg", "Zscore", "Spm", "Ee", "Pem"), add, c(ceiling(max(orgPSI$Max)), 0.95, 2, 0.5), ceiling(max(orgPSI$Max)))

  return()
}
###***###***###

