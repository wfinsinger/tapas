# Identify charcoal peaks based on possible thresholds ####

# based on 
# function [Charcoal, CharThresh] = CharPeakID (Charcoal,Pretreatment,PeakAnalysis,CharThresh);
#   Identifies charcoal samples that exceeds threshold value(s) determined 
#   in CharThrshLocal or CharThreshGlobal, and screens these value
#   according to the minimum-count criterion selected.

SeriesPeaks <- function(x, minCountP = NULL) {
  
  if (is.null(minCountP)) minCountP <- 0.05
  
  # Get data needed:
  Charcoal.peak <- as.vector(x$Charcoal.peak)
  thresh.values <- x$thresh.values
  nThresholds <- length(thresh.values)
  thresh.pos <- x$thresh.pos
  yr.interp <- x$yrInterp
  ybpI <- x$ybpI
  countI <- x$countI
  volI <- x$volI
  
  ## PEAK IDENTIFICATION ALGORITHM
  # Create space for peaks, Charcoal.charPeaks
  # For local thresholds only!
  Charcoal.charPeaks <- as.data.frame(matrix(data=0,
                                             nrow=length(Charcoal.peak),
                                             ncol=length(thresh.pos)))
  
  
  # Flag values exceeding thresholds
  for (i in 1:length(Charcoal.peak)) {  # For each value in Charcoal.peak
    for (j in 1:nThresholds) {          # For each threshold value
      if (Charcoal.peak[i] > thresh.pos[i, j]) {  # If Charcoal.peak exceeds threshold...
        Charcoal.charPeaks[i, j] <- 2                  # Charcoal.charPeaks = 2
      } else {
        Charcoal.charPeaks[i, j] <- 0                 # else Charcoal.charPeaks = 0
      }
    }
  }
  
  # Remove consecutive Charcoal.charPeaks
  for (i in 1:(length(Charcoal.peak)-1)) { # For each value in Charcoal.peak
    for (j in 1:nThresholds) {             # For each threshold value
      if (Charcoal.charPeaks[i, j] > 0
          && Charcoal.charPeaks[i+1, j] > 0) {  # if two consecutive values > 0 
        Charcoal.charPeaks[i, j] <- 1           # keep first as 2, mark subsequent as 1
      }
    }
  }
  
  for (i in 1:length(Charcoal.peak)) {
    for (j in 1:nThresholds) {
      if (Charcoal.charPeaks[i, j] < 2) {    # if value < 2
        Charcoal.charPeaks[i, j] <- 0        # mark sample as 0 (unflag Charcoal.charPeak)
      } else {
        Charcoal.charPeaks[i, j] <- 1        # else (if value=2) mark sample as 1 (flag as Charcoal.charPeak)
      }
    }
  }
  
  
  # Make a variable to hold the threshold value for each peak identified.
  # i.e. instead of 1 / 0 for peak / no peak, 1 is replaced with the 
  # threshold value used to identify the peak. This is purely for graphing 
  # purposes.
  Charcoal.charPeaksThresh <- Charcoal.charPeaks * 0
  for (i in 1:length(Charcoal.peak)) {
    for (j in 1:nThresholds) {
      Charcoal.charPeaksThresh[i, j] <- Charcoal.charPeaks[i,j] * thresh.pos[i,j]
    }
  }
  
  
  ## Minimum-count Analysis
  # Set [yr] Years before and after a peak to look for the min. and max. value
  mcWindow <- round(150/yr.interp)*yr.interp
  
  # Create space
  d <- as.data.frame(matrix(data=0, nrow=length(x$accI), ncol=nThresholds))
  CharThresh.minCountP <- as.data.frame(matrix(data=NA, nrow=length(x$acc), ncol=nThresholds))
  
  for (j in 1:nThresholds) {
    peakIndex <- which(Charcoal.charPeaks[ ,j] == 1) # Index to find peak samples
    
    if (length(peakIndex) > 1) {                     # Only proceed if there is > 1 peak
      for (i in 1:length(peakIndex)) {               # For each peak identified...
        peakYr <- ybpI[peakIndex[i]]      # Find age of peak and window around peak
        windowTime <- c(max(ybpI[which(ybpI <= peakYr+mcWindow)]),
                        min(ybpI[which(ybpI >= peakYr-mcWindow)]))
        windowTime_in <- c(which(ybpI == windowTime[1]), # Index to find range of window ages
                           which(ybpI == windowTime[2]))
        if (i == 1) {  # find the year of two adjacent Charcoal.charPeaks, unless first peak,
          #then use windowTime[2] as youngest
          windowPeak_in <- c(which(ybpI == ybpI[peakIndex[i+1]]),
                             which(ybpI == windowTime[2]))
        }
        if (i == length(peakIndex)) {  # if last peak, then use windowTime[1] as oldest age
          windowPeak_in <- c(which(ybpI == windowTime[1]),
                             which(ybpI == ybpI[peakIndex[i-1]]))
        }
        if (i > 1 && i < length(peakIndex)) {
          windowPeak_in <- c(which(ybpI == ybpI[peakIndex[i+1]]),
                             which(ybpI == ybpI[peakIndex[i-1]]))
        }
        if (windowTime_in[1] > windowPeak_in[1]) { # thus, if a peak falls within the time window defined by mcWindow
          windowTime_in[1] <- windowPeak_in[1] # replace the windowTime_in with the windowPeak_in
        }
        if (windowTime_in[2] < windowPeak_in[2]) { # thus, if a peak falls within the time window defined by mcWindow
          windowTime_in[2] <- windowPeak_in[2] # replace the windowTime_in with the windowPeak_in
        }
        
        # Final index value for search window: window (1) defines oldest sample,
        # window (2) defines youngest sample
        windowSearch <- c(windowTime_in[1], windowTime_in[2])
        
        # search for max and min Charcoal.charPeaks within this window.
        # [# cm^-3] Max charcoal concentration after peak.
        countMax <- round(max(countI[ windowSearch[2]:peakIndex[i] ]))
        # Index for location of max count.
        countMaxIn <- windowSearch[2]-1 + max(which(round(countI[windowSearch[2]:peakIndex[i]]) == countMax))
        # [# cm^-3] Min charcoal concentration before peak.
        countMin <- round(min(countI[peakIndex[i]:windowSearch[1] ]))
        # Index for location of Min count
        countMinIn <- peakIndex[i]-1 + min(which(round(countI[peakIndex[i]:windowSearch[1]]) == countMin))
        
        volMax <- volI[countMaxIn]
        volMin <- volI[countMinIn]
        d[peakIndex[i], j] <- (abs(countMin-(countMin+countMax)*
                                     (volMin/(volMin+volMax)))-0.5)/(sqrt((countMin+countMax)*
                                                                            (volMin/(volMin+volMax))*(volMax/(volMin+volMax))))
        
        # Test statistic
        CharThresh.minCountP[peakIndex[i], j] <- 1-pt(d[peakIndex[i], j], df=Inf)
        # Inverse of the Student's T cdf at 
        # CharThresh.minCountP, with Inf degrees of freedom.
        # From Charster (Gavin 2005):
        # This is the expansion by Shuie and Bain (1982) of the equation by 
        # Detre and White (1970) for unequal 'frames' (here, sediment 
        # volumes). The significance of d is based on the t distribution 
        # with an infinite degrees of freedom, which is the same as the 
        # cumulative normal distribution.
      }
    }  
  }
  
  # Clean Environment
  rm(mcWindow, d, countMax, countMaxIn, countMin, countMinIn, peakIndex, peakYr, volMax, volMin,
     windowPeak_in, windowSearch, windowTime, windowTime_in)
  
  
  # Take note of and remove peaks that do not pass the minimum-count screening-peak test
  Charcoal.charPeaks.insig <- as.data.frame(matrix(data=0, nrow=length(Charcoal.peak), ncol=length(thresh.pos)))
  for (j in 1:nThresholds) {
    insig.peaks <- intersect(which(Charcoal.charPeaks[ ,j] > 0),
                             which(CharThresh.minCountP[ ,j] > minCountP)) # Index for
    # Charcoal.charPeaks values that also have p-value > minCountP...thus insignificant
    Charcoal.charPeaks.insig[insig.peaks, j] <- 1
    Charcoal.charPeaks[insig.peaks, j] <- 0 # set insignificant peaks to 0
    Charcoal.charPeaksThresh[insig.peaks, j] <- 0
  }
  
  # Calculate sensitivity indices:
  # - FRI with different thresholds
  Charcoal.peaksTotal <- matrix(0, nrow=1, ncol=4)
  Charcoal.threshFRI <- matrix(NA, nrow=max(colSums(Charcoal.charPeaks))-1, ncol=4)
  for (j in 1:nThresholds) {
    Charcoal.peaksTotal[j] <- sum(Charcoal.charPeaks[ ,j])
    inFRI <- diff(ybpI[which(Charcoal.charPeaks[ ,j] == 1)])
    if (length(inFRI) > 0) {
      Charcoal.threshFRI[1:length(inFRI), j] <- inFRI
    }
  }
  
  # Write output to Environment
  x <- unclass(x)
  output <- structure(list(CHAR.sm=x, Charcoal.charPeaks=Charcoal.charPeaks,
                           Charcoal.insigPeaks=Charcoal.charPeaks.insig,
                           threshFRI=Charcoal.threshFRI, minCountP=minCountP))
  class(output) <- "CharPeaks"
  return(output)
}




plot.CharPeaks <- function(x, xlim=NULL, ylim=NULL,xlab=NULL,...) {
  
  # x <- CO.peak
  # xlim=NULL
  # ylim=NULL
  
  # Get data
  ageI <- x$CHAR.sm$ybpI
  accI <- x$CHAR.sm$accI
  cBack <- x$CHAR.sm$charBack
  SNI <- x$CHAR.sm$SNI
  SNIsm <- x$CHAR.sm$SNIsm
  thresh.pos <- x$CHAR.sm$thresh.pos$X4
  thresh.neg <- x$CHAR.sm$thresh.neg$X4
  Charcoal.peak <- x$CHAR.sm$Charcoal.peak
  sig.peaks <- x$Charcoal.charPeaks$V4
  sig.peaks <- which(sig.peaks > 0)
  insig.peaks <- x$Charcoal.insigPeaks$V4
  insig.peaks <- which(insig.peaks > 0)
  
  # Set axis limits
  if(is.null(xlim)){
    x.min <- floor(min(ageI)/100) * 100
    x.max <- ceiling(max(ageI)/100) * 100
    if(x.max - x.min > 1000) {
      x.by <- 1000
    } else {
      x.by <- round((x.max - x.min)/100) * 100
    }
    x.lim <- c(x.max, x.min)
  }
  if(is.null(ylim)) {
    y.lim <- c(min(accI), 1.2*max(accI))
  }
  peaks.ylim <- c(min(Charcoal.peak), max(Charcoal.peak))
  
  
  
  par(mfrow=c(2,1), mar=c(0.5, 4, 0.5, 1), oma=c(5,1,1,1), cex=0.7)
  plot(ageI, accI, type="h", col=grey(0.7), xlim=x.lim, ylim=y.lim,
       ylab="CHAR (# cm^-2 yr^-1)", axes=F)
  lines(ageI, cBack, lwd=1.5, col="red")
  lines(ageI, cBack+thresh.pos, lwd=0.5, col="red")
  axis(1, at=seq(0, x.max, by=x.by), labels=F)
  axis(2)
  legend("topleft", inset=c(0,0.2),
         legend=c("CHAR interpolated", "CHAR background", "Thresholds"),
         lwd=c(1.5, 1.5, 0.5), col=c("grey", "red", "red"),
         border="NA", bty="n", horiz=F, cex=0.7)
  
  plot(ageI, Charcoal.peak, type="h", col=grey(0.7), xlim=x.lim,
       xlab="", ylab="CHAR residuals (# cm^-2 yr^-1)", lwd=1, axes=F)
  lines(ageI, thresh.pos, type="l", col="red", lwd=1)
  lines(ageI, thresh.neg, type="l", col="red", lwd=1)
  points(ageI[insig.peaks], rep(0.85*peaks.ylim[2], length(insig.peaks)), type="p",
         pch=19, col="grey")
  points(ageI[sig.peaks], rep(0.9*peaks.ylim[2], length(sig.peaks)), type="p", pch=3)
  legend("topleft", inset=c(0,0.2),
         legend=c("CHARi residuals", "Thresholds"),
         lwd=c(1.5, 1), col=c("grey", "red"),
         border="NA", bty="n", horiz=F, cex=0.7)
  axis(1, at=seq(0, x.max, by=x.by))
  axis(2)
  mtext(paste("b) CHARinterpolated and different options for a yr CHARbackground",
              sep=''), cex=0.8, line=-2)
}


plot.SNI <- function(x, xlim=NULL, ylim=NULL,xlab=NULL,...) {
  
  # Get data
  ageI <- x$CHAR.sm$ybpI
  accI <- x$CHAR.sm$accI
  cBack <- x$CHAR.sm$charBack
  SNI <- x$CHAR.sm$SNI
  SNIsm <- x$CHAR.sm$SNIsm
  thresh.pos <- x$CHAR.sm$thresh.pos$X4
  thresh.neg <- x$CHAR.sm$thresh.neg$X4
  Charcoal.peak <- x$CHAR.sm$Charcoal.peak
  sig.peaks <- x$Charcoal.charPeaks$V4
  sig.peaks <- which(sig.peaks > 0)
  insig.peaks <- x$Charcoal.insigPeaks$V4
  insig.peaks <- which(insig.peaks > 0)
  
  # Set axis limits
  if(is.null(xlim)){
    x.min <- floor(min(ageI)/100) * 100
    x.max <- ceiling(max(ageI)/100) * 100
    if(x.max - x.min > 1000) {
      x.by <- 1000
    } else {
      x.by <- round((x.max - x.min)/100) * 100
    }
    x.lim <- c(x.max, x.min)
  }
  if(is.null(ylim)) {
    y.lim <- c(min(accI), 1.2*max(accI))
  }
  peaks.ylim <- c(min(Charcoal.peak), max(Charcoal.peak))
  
  
  par(mfrow=c(2,1), mar=c(0.5, 4, 0.5, 1), oma=c(5,1,1,1), cex=0.7)
  plot(ageI, accI, type="h", col=grey(0.7), xlim=x.lim, ylim=y.lim,
       xlab="time (cal yr. BP)", ylab="CHAR (# cm^-2 yr^-1)", axes=F)
  lines(ageI, cBack, lwd=1.5, col="red")
  lines(ageI, cBack+thresh.pos, lwd=0.5, col="red")
  points(ageI[insig.peaks], rep(0.9*y.lim[2], length(insig.peaks)), type="p",
         pch=19, col="grey")
  points(ageI[sig.peaks], rep(0.99*y.lim[2], length(sig.peaks)), type="p", pch=3)
  axis(1, at=seq(0, x.max, by=x.by), labels=F)
  axis(2)
  
  plot(ageI, SNI, type="l", col=grey(0.7), xlim=x.lim, ylim=c(0, max(SNI)),
       xlab="time (cal yr. BP)", ylab="SNI", axes=F)
  lines(ageI, SNIsm, lwd=1.5, col="black")
  abline(h=3, lty="dashed")
  axis(1, at=seq(0, x.max, by=x.by), labels=T)
  axis(2)
  
  boxplot(x$threshFRI[ ,1:3], notch=T, boxwex=0.2, names=x$CHAR.sm$thresh.values[1:3],
          ylab="FRI (years)")
}

