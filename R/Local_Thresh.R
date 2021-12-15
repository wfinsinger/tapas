#  ***************************************************************************
#   Local_Thresh.R
#  ---------------------
#   Date                 : December 2021
#   Copyright            : (C) 2021 by Walter Finsinger
#   Email                : walter.finsinger@umontpellier.fr
#  ---------------------
#
#  ***************************************************************************
#
#   The script determines threshold values to decompose a detrended timeseries into
#     a 'noise' and a 'signal' component using a 2-component Gaussian Mixture Model (GMM).
#     Based on Phil Higuera's CharThreshLocal.m Matlab code.
#   Determines a positive and a negative threshold value for each interpolated sample, based on the
#     distribution of values within the selected window (thresh.yr).
#   The procedure uses a Gaussian mixture model on the assumption that the noise component
#     is normally distributed around 0 (the values were detrended!).
#
#   Requires an output from the SeriesDetrend() function.
#
# The user-defined parameters are as follows:
#   series  ->  the output from the SeriesDetrend() function
#   
#   proxy   ->  Set proxy="VariableName" select the variable for the peak-detection analysis.
#               If the dataset includes only one variable, proxy does not need to be specified. 
#   
#   t.lim   ->  allows defining a portion of the time series.
#     With t.lim=NULL (as per default), the analysis will be performed using
#     the entire timeseries.
#   
#   thresh.value  ->  Determines the threshold as the nth-percentile of the
#     Gaussian Model of the noise component. Default thresh.value = 0.95
#   
#   thresh.yr     =>  determines the length of the window width (in years) from which
#                       values are selected to determine the local threshold. By default,
#                       this value is inherited from the smoothing.yr value set in
#                       the SeriesDetrend() function.
#   
#   smoothing.yr  =>  Width of moving window for computing SNI
#                     
#   keep_consecutive => logical. If FALSE (default), consecutive peak samples exceeding the
#                       threshold will be removed and only the first (older) sample is retained.
#                       
#   MinCountP       => Defines the probability (default = 0.05) that two resampled counts
#                     could arise from the same Poisson distribution. This is used to
#                     screen peak samples and remove any that fail to pass the 
#                     minimum-count test. If MinCountP = NULL, the test will not be performed.
#
#   MinCountP_window = Defines the width (in years) of the search window used for the 
#                     minimum-count test. Default: MinCountP_window=150.
#   
#   out.dir = path to the output directory where figures will be written as *.pdf
#
#  ***************************************************************************


local_thresh <- function(series = NA, proxy = NULL, t.lim = NULL, thresh.yr = NULL,
                         thresh.value = 0.95, smoothing.yr = NULL,
                         keep_consecutive = F,
                         minCountP = 0.05, MinCountP_window = 150,
                         out.dir = "Figures") {
  
  require(mclust)
  
  # Initial check-up of input parameters ####
  if (keep_consecutive == T & is.null(minCountP) == F) {
    print('Fatal error: inconsistent choice of arguments. If keep_consecutive=T, set minCountP=NULL.')
    return()
  }
  
  # Determine path to output folder for Figures
  out.path <- paste0("./", out.dir, "/")
  
  if (dir.exists(out.path) == F) {
    dir.create(out.path)
  }
  
  
  # Extract parameters from input list (i.e. detrended series) ####
  
  # Determines for which of the variables the threshold analysis should be made
  if (is.null(proxy) == T) { # if proxy = NULL, use the data in series$detr$detr
    if (dim(series$detr$detr)[2] > 2) {
      print('Fatal error: please specify which proxy you want to use')
      return()
    } else {
      proxy <- colnames(series$detr$detr) [2]
      a <- series$detr$detr
    }
  }
  
  if (is.null(proxy) == F) { # if we want to select a specific variable
    a <- data.frame(series$detr$detr$age)
    colnames(a) <- "age"
    a <- cbind(a, series$detr$detr[[proxy]])
    colnames(a)[2] <- proxy
  }
  
  # Checks if smoothing.yr and thresh.yr were specified,
  # else use smoothing.yr used for detrending
  if (is.null(smoothing.yr) == T) {
    smoothing.yr <- series$detr$smoothing.yr
  }
  
  if (is.null(thresh.yr) == T) {
    thresh.yr <- series$detr$smoothing.yr
  }
  
  # Further extract parameters from input dataset
  a.names <- colnames(a)
  s.name <- series$detr$series.name
  yr.interp <- series$int$yr.interp
  
  # Determine whether to limit the analysis to a portion of the record
  if (is.null(t.lim) == T) {
    ageI <- a$age
    t.lim <- c(min(ageI), max(ageI))
  }
  if (is.null(t.lim) == F) {
    a <- a[which(a$age <= max(t.lim) & a$age >= min(t.lim)), ]
    ageI <- a$age[which(a$age <= max(t.lim) & a$age >= min(t.lim))]
  }
  x.lim <- c(max(ageI), min(ageI))
  
  
  # Determine the proportion of datapoints used to smooth local thresholds with loess()
    n.smooth <- round(smoothing.yr/yr.interp)
    span.sm <- n.smooth/dim(a)[1]

  
  # Create empty list where output data will be stored
  a.out <- list(span.sm = span.sm, thresh.value = thresh.value, yr.interp = yr.interp)
  
  
  # Create space for local variables
  v <- a[ ,2] # variable
  v.name <- colnames(a)[2]
  v.gmm <- v[which(complete.cases(v))]
  y.lim <- c(min(v.gmm), max(v.gmm))
  
  
  thresh.pos <- matrix(data = NA, nrow = length(v)) # threshold values
  thresh.neg <- matrix(data = NA, nrow = length(v)) # threshold values
  muHat <- data.frame(matrix(data = NA, nrow = length(v), ncol = 2)) # mean of noise distribution
  sigmaHat <- data.frame(matrix(data = NA, nrow = length(v), ncol = 2)) # standard deviation of noise distribution
  propN <- data.frame(matrix(data = NA, nrow = length(v), ncol = 2)) # proportion of each Cluster-identified distribution
  Peaks.pos <- matrix(data = 0, nrow = length(v), ncol = 1) # positive peaks
  Peaks.neg <- matrix(data = 0, nrow = length(v), ncol = 1) # negative peaks
  
  # Define number of plots to print for evaluation of local threshold
  j <- 1
  
  num.plots <- seq(from = round(n.smooth), to = length(v), by = round(n.smooth/2))
  my.plots <- vector(length(num.plots), mode = 'list')
  
  
  
  # Determine local thresholds with Gaussian mixture models (2 components) ####
  
  # SELECT peak VALUES TO EVALUATE, BASED ON Smoothing.yr
  for (i in 1:length(v)) {  #For each value in the detrended dataseries, find the threshold.
    #i=3
    age.i <- a[i, 1]
    #cat(paste0("Calculating ", i, "th local threshold of ", length(v)))
    #i=28
    if (i < round(0.5*(thresh.yr/yr.interp)) + 1) { # First 'thresh.yr' samples.
      #         X = v(1:round(0.5*(thresh.yr/r))); % Pre June 2009.
      #         X = v(1:round(thresh.yr/r)); % Modified, June 2009, PEH.
      X <- v[1:round(0.5*(thresh.yr/yr.interp)) + i] # Modified, % June 2009, PEH.
    }
    if (i > (length(v) - round(0.5*(thresh.yr/yr.interp)))) {  # Last 'thresh.yr' samples.
      #             X = v(length(v)-...
      #                 round((thresh.yr/r)):end);   % Pre June 2009.
      X <- v[(i - round(0.5*(thresh.yr/yr.interp))):length(v)]   # Modified, June 2009, PEH.
      # As recommended by RK, this uses samples from 
      # a half-window before i, all the way to end of record.
    }
    if (i >= round(0.5*(thresh.yr/yr.interp)) + 1 && i <= (length(v) - round(0.5*(thresh.yr/yr.interp)))) {
      # All samples between first and last 'thrshYr' samples.
      X <- v[(i - round(0.5*(thresh.yr/yr.interp))):(i + round(0.5*(thresh.yr/yr.interp)))]
    }
    
    

    ## Estimate local noise distribution with Guassian mixture model ####
    
    X.gmm <- X[which(complete.cases(X))]
    
    if (length(X) < 3) {
      cat("NOTE: Less than 3 peak values in moving window; cannot fit noise distribution.")
      cat("\n      Mean and standard deviation forced to equal 0.") 
      cat("\n      Consider longer smoothing window (thresh.yr).")
      muHat[i, ] <- 0
      sigmaHat[i, ] <- 10^-100
      propN[i, ] <- 0
      thresh[i] <- 0
      Thresh.SNI[i] <- 0
    } else {
      m <- mclust::densityMclust(data = X.gmm, G = 2, verbose = F)
      # plot(x=m, what="density", data=X)
      # summary(m, parameters=T, classification=T)
      # plot.Mclust(m, what="classification")
      
      # Get parameters from Gaussian mixture models
      muHat[i, ] <- m$parameters$mean
      sigmaHat[i, 1] <- sqrt(m$parameters$variance$sigmasq[1])
      sigmaHat[i, 2] <- sqrt(m$parameters$variance$sigmasq[2])
      propN[i, ] <- m$parameters$pro
      
      # Check
      if (muHat[i, 1] == muHat[i, 2]) {
        warning('Poor fit of Gaussian mixture model')
      }
      
      ## Define local threshold
      if (!is.na(m$parameters$variance$sigmasq[1]) == T) {
        thresh.pos[i] <- m$parameters$mean[1] + qnorm(p = thresh.value) *
          sqrt(m$parameters$variance$sigmasq)[1]
        thresh.neg[i] <- m$parameters$mean[1] - qnorm(p = thresh.value) *
          sqrt(m$parameters$variance$sigmasq)[1]
        
        if (m$parameters$mean[1] < 0) {
          thresh.pos[i] <- 0 + qnorm(p = thresh.value) *
            sqrt(m$parameters$variance$sigmasq)[1]
        }
        
        #sig_i.pos <- X[which(X > thresh.pos[i])]
        #noise_i <- X[which(X <= thresh.pos[i] & X >= thresh.neg[i])]
        #sig_i.neg <- X[which(X < thresh.neg[i])]
        
        
      }
      
      if (is.na(m$parameters$variance$sigmasq[1]) == T) {
        thresh.pos[i] <- 0
        thresh.neg[i] <- 0
      }
      
    }
    
    
    # Plot some of the noise and signal distributions
    if (any(num.plots == i)) {
      # Print plots for positive peaks
      par(mfrow = c(2,1), mar = c(5,4,1,1))
      
      plot.Mclust(m, what = "classification")
      
      h <- hist(x = X, breaks = 50, plot = F)
      dens <- hist(X, breaks = 50, plot = F)$density
      plot(h, freq = F, col = "grey", border = "grey",
           xlim = c(min(X, na.rm = T), max(X, na.rm = T)),
           ylab = "Density", xlab = '', main = "Local threshold")
      # plot(h, freq = F, col = "grey", border = "grey", xlim = c(min(h$breaks), max(h$breaks)),
      #      ylab = "Density", xlab = '', main = "Local threshold")
      par(new = T)
      pdf2 <- curve(dnorm(x = x, mean = muHat[i, 2], sd = sigmaHat[i, 2]),
                    from = min(h$breaks), to = max(h$breaks),
                    ylim = c(0, max(dens)), type = "l", col = "blue", lwd = 1.5,
                    axes = F, ylab = '', xlab = '')
      par(new = T)
      pdf1 <- curve(dnorm(x = x, mean = muHat[i, 1], sd = sigmaHat[i, 1]),
                    from = min(h$breaks), to = max(h$breaks),
                    ylim = c(0, max(dens)), type = "l", col = "orange", lwd = 1.5,
                    axes = F, ylab = '', xlab = '')
      par(new = T)
      lines(x = c(thresh.pos[i], thresh.pos[i]), y = c(0, max(dens)),
            type = "l", col = "red", lwd = 1.5)
      lines(x = c(thresh.neg[i], thresh.neg[i]), y = c(0, max(dens)),
            type = "l", col = "red", lwd = 1.5)
      mtext(paste0(age.i, " years", "; thresh.value = ", thresh.value),
            side = 3, las = 0, line = -1)
      # mtext(text=paste0("SNIi pos.= ", round(Thresh.SNI.pos[i], digits=2),
      #                   "SNIi neg.= ", round(Thresh.SNI.neg[i], digits=2)), 
      #       side=3, las=0, line=-2)
      
      my.plots[[j]] <- recordPlot()
      j <- j + 1
    }
  }
  
  # Print pdf with selected plots that were saved at the end of the loop above
  pdf(paste0(out.path, s.name, '_', v.name, '_Local_GMM_Evaluation.pdf'),
      onefile = TRUE, paper = "a4")
  par(mfrow = (c(5,5)), mar = c(0.5,1,0.5,1), oma = c(1,1,0.5,1), cex = 0.7)
  for (k in 1:length(num.plots)) {
    replayPlot(my.plots[[k]])
  }
  dev.off()
  
  #rm(my.plots, num.plots)
  
  
  ## Calculate SNI ####
  ## 
  ## Get resampled values for selected variable (proxy) for selected time interval (t.lim)
  SNI_in_index <- which(series$int$series.int$age <= max(t.lim) & 
                          series$int$series.int$age >= min(t.lim))
  
  SNI_in <- series$int$series.int[[proxy]] [SNI_in_index]
  
  
  SNI_pos <- SNI(ProxyData = cbind(ageI, SNI_in, thresh.pos),
                 BandWidth = smoothing.yr)
  
  ## If SNI = NA (in case all values are > threshold, or all values < threshold)
  if (all(is.na(SNI_pos$SNI_raw))) {
    SNI_pos$SNI_raw <- 0
    SNI_pos$SNI_sm <- 0
  }
  
  
  SNI_neg <- SNI(ProxyData = cbind(ageI, -1 * SNI_in, -1 * thresh.neg),
                 BandWidth = smoothing.yr)
  
  ## If SNI = NA (in case all values are > threshold, or all values < threshold)
  if (all(is.na(SNI_neg$SNI_raw))) {
    SNI_neg$SNI_raw <- 0
    SNI_neg$SNI_sm <- 0
  }
  
  rm(SNI_in, SNI_in_index)
  
  
  
  # Smooth local threshold values ####
  ## Smooth thresholds Lowess smoother
  # thresh.pos.sm <- stats::loess(thresh.pos ~ ageI, data = data.frame(ageI, thresh.pos),
  #                        span = span.sm, family = "gaussian")$fitted
  # thresh.neg.sm <- stats::loess(thresh.neg ~ ageI, data = data.frame(ageI, thresh.neg),
  #                        span = span.sm, family = "gaussian")$fitted
  
  thresh.pos.sm <- stats::lowess(ageI, thresh.pos, f = span.sm, iter = 0)$y
  thresh.neg.sm <- stats::lowess(ageI, thresh.neg, f = span.sm, iter = 0)$y
  
  
  
  ## Get peaks ####
  
  ## If keep_consecutive == F ####
  ## 
  if (keep_consecutive == F) {  # if consecutive peak samples should be removed
    
    # Flag values exceeding thresholds
    Peaks.pos[which(v > thresh.pos.sm)] <- 2
    Peaks.neg[which(v < thresh.neg.sm)] <- 2
    
    # Then remove consecutive peaks
    # For positive peaks
    for (i in 1:(length(Peaks.pos) - 1)) { # For each value in Charcoal.peak
      if (Peaks.pos[i] > 0
          && Peaks.pos[i + 1] > 0) {  # if two consecutive values > 0 
        Peaks.pos[i] <- 1           # keep first as 2, mark subsequent (earlier) as 1
      }
    }
    
    for (i in 1:length(Peaks.pos)) {
      if (Peaks.pos[i] < 2) {    # if value < 2
        Peaks.pos[i] <- 0        # mark sample as 0 (unflag Peak)
      } else {
        Peaks.pos[i] <- 1        # else (if value=2) mark sample as 1 (flag as Peak)
      }
    }
    
    
    # For negative peaks
    for (i in 1:(length(Peaks.neg) - 1)) { # For each value in Charcoal.peak
      if (Peaks.neg[i] > 0
          && Peaks.neg[i + 1] > 0) {  # if two consecutive values > 0 
        Peaks.neg[i] <- 1           # keep first as 2, mark subsequent (earlier) as 1
      }
    }
    
    for (i in 1:length(Peaks.neg)) {
      if (Peaks.neg[i] < 2) {    # if value < 2
        Peaks.neg[i] <- 0        # mark sample as 0 (unflag Peak)
      } else {
        Peaks.neg[i] <- 1        # else (if value=2) mark sample as 1 (flag as Peak)
      }
    }
  } else {
    Peaks.pos[which(v > thresh.pos.sm)] <- 1
    Peaks.neg[which(v < thresh.neg.sm)] <- 1
  }
  
  
  
  ## Minimum-count Analysis ####
  ## If minCountP is not NULL, screen positive peaks with minimum count test ####
  if (is.null(minCountP) == F) {
    ## Minimum-count Analysis
    
    # Set [yr] Years before and after a peak to look for the min. and max. value
    if (is.null(MinCountP_window) == T) {
      MinCountP_window <- round(150/yr.interp) * yr.interp
    }
    
    countI_index <- which(colnames(series$int$series.int) == proxy)
    countI <- series$int$series.conI[ ,countI_index]
    volI <- series$int$volI
    
    # Create space
    d <- rep_len(NA, length.out = dim(a) [1])
    Thresh_minCountP <- rep_len(NA, length.out = dim(a) [1])
    
    peakIndex <- which(Peaks.pos == 1) # Index to find peak samples
    
    if (length(peakIndex) > 1) {                     # Only proceed if there is > 1 peak
      for (i in 1:length(peakIndex)) {               # For each peak identified...
        peakYr <- ageI[peakIndex[i]]      # Find age of peak and window around peak
        windowTime <- c(max(ageI[which(ageI <= peakYr + MinCountP_window)]),
                        min(ageI[which(ageI >= peakYr - MinCountP_window)]))
        windowTime_in <- c(which(ageI == windowTime[1]), # Index to find range of window ages
                           which(ageI == windowTime[2]))
        if (i == 1) {  # find the year of two adjacent Peaks.pos, unless first peak,
          #then use windowTime[2] as youngest
          windowPeak_in <- c(which(ageI == ageI[peakIndex[i + 1]]),
                             which(ageI == windowTime[2]))
        }
        if (i == length(peakIndex)) {  # if last peak, then use windowTime[1] as oldest age
          windowPeak_in <- c(which(ageI == windowTime[1]),
                             which(ageI == ageI[peakIndex[i - 1]]))
        }
        if (i > 1 && i < length(peakIndex)) {
          windowPeak_in <- c(which(ageI == ageI[peakIndex[i + 1]]),
                             which(ageI == ageI[peakIndex[i - 1]]))
        }
        if (windowTime_in[1] > windowPeak_in[1]) { # thus, if a peak falls within the time window defined by MinCountP_window
          windowTime_in[1] <- windowPeak_in[1] # replace the windowTime_in with the windowPeak_in
        }
        if (windowTime_in[2] < windowPeak_in[2]) { # thus, if a peak falls within the time window defined by MinCountP_window
          windowTime_in[2] <- windowPeak_in[2] # replace the windowTime_in with the windowPeak_in
        }
        
        # Final index value for search window: window (1) defines oldest sample,
        # window (2) defines youngest sample
        windowSearch <- c(windowTime_in[1], windowTime_in[2])
        
        # search for max and min Peaks.pos within this window.
        # [# cm^-3] Max charcoal concentration after peak.
        countMax <- round(max(countI[ windowSearch[2]:peakIndex[i] ]))
        
        # Index for location of max count.
        countMaxIn <- windowSearch[2] - 1 + max(which(round(countI[windowSearch[2]:peakIndex[i]]) == countMax))
        
        # [# cm^-3] Min charcoal concentration before peak.
        countMin <- round(min(countI[peakIndex[i]:windowSearch[1] ]))
        
        # Index for location of Min count
        countMinIn <- peakIndex[i] - 1 + min(which(round(countI[peakIndex[i]:windowSearch[1]]) == countMin))
        
        volMax <- volI[countMaxIn]
        volMin <- volI[countMinIn]
        d[ peakIndex[i] ] <- (abs(countMin - (countMin + countMax) *
                                    (volMin/(volMin + volMax))) - 0.5)/(sqrt((countMin + countMax) *
                                                                               (volMin/(volMin + volMax)) * 
                                                                               (volMax/(volMin + volMax))))
        
        # Test statistic
        Thresh_minCountP[peakIndex[i]] <- 1 - pt(q = d[peakIndex[i]], df = Inf)
        # Inverse of the Student's T cdf at 
        # Thresh_minCountP, with Inf degrees of freedom.
        # From Charster (Gavin 2005):
        # This is the expansion by Shuie and Bain (1982) of the equation by 
        # Detre and White (1970) for unequal 'frames' (here, sediment 
        # volumes). The significance of d is based on the t distribution 
        # with an infinite degrees of freedom, which is the same as the 
        # cumulative normal distribution.
      }
    }
    
    # Clean Environment
    rm(MinCountP_window, d, countMax, countMaxIn, countMin, countMinIn, peakIndex,
       peakYr, volMax, volMin, windowPeak_in, windowSearch, windowTime, windowTime_in)
    
    
    # Take note of and remove peaks that do not pass the minimum-count screening-peak test
    Peaks.pos.insig <- rep_len(0, length.out = dim(a) [1])
    
    insig.peaks <- intersect(which(Peaks.pos > 0),
                             which(Thresh_minCountP > minCountP)) # Index for
    # Peaks.pos values that also have p-value > minCountP...thus insignificant
    Peaks.pos.insig[insig.peaks] <- 1
    Peaks.pos[insig.peaks] <- 0 # set insignificant peaks to 0
    #Peaks.posThresh[insig.peaks] <- 0
    
  } else {
    insig.peaks <- NULL
    Peaks.pos.insig <- rep_len(NA, length.out = dim(a) [1])
  }
  
  
  
  
  
  ## Gather data to calculate return intervals ####    
    
  ## Plot series with trend + threshold + peaks
  Peaks.pos.final <- which(Peaks.pos == 1)
  Peaks.neg.final <- which(Peaks.neg == 1)
  
  peaks.pos.ages <- ageI[which(Peaks.pos == 1)]
  peaks.neg.ages <- ageI[which(Peaks.neg == 1)]
  
  
  ## Calculate Event Return Intervals
  RI_neg <- c(diff(peaks.neg.ages), NA)
  RI_pos <- c(diff(peaks.pos.ages), NA)
  
  
  
  ## Get peak magnitude ####
  # Peak magnitude (pieces cm-2 peak-1) is the sum of all samples exceeding
  # threshFinalPos for a given peak.
  # The units are derived as follows:
  # [pieces cm-2 yr-1] * [yr peak-1] = [pieces cm-2 peak-1].  
  
  ## Get peak-magnitude values for positive peaks
  if (length(Peaks.pos.final) > 0) {
    PeakMag_pos_val <- v - thresh.pos
    PeakMag_pos_val[PeakMag_pos_val < 0] <- 0
    PeakMag_pos_index <- which(PeakMag_pos_val > 0)
    
    Peaks_pos_groups <- split(PeakMag_pos_index,
                              cumsum(c(1, diff(PeakMag_pos_index) != 1)))
    
    PeakMag_pos <- vector(mode = "numeric", length = length(Peaks_pos_groups))
    PeakMag_pos_ages <- vector(mode = "numeric", length = length(Peaks_pos_groups))
    Peak_duration <- vector(mode = "numeric", length = length(Peaks_pos_groups))
    
    for (i in 1:length(Peaks_pos_groups)) {
      #i = 2
      n_groups <- length(Peaks_pos_groups[[i]])
      Peak_duration[i] <- ageI[ Peaks_pos_groups[[i]] [n_groups] + 1] -
        ageI[ Peaks_pos_groups[[i]] [1]]
      PeakMag_pos[i] <- sum(c(PeakMag_pos_val[ Peaks_pos_groups[[i]] ])) *
        Peak_duration[i]
      PeakMag_pos_ages[i] <- ageI[Peaks_pos_groups[[i]] [n_groups]]
    }
    
    PeakMag_pos <- as.data.frame(cbind(PeakMag_pos_ages, PeakMag_pos))
    colnames(PeakMag_pos) <- c("ageI", "peak_mag")
    
    rm(PeakMag_pos_val, PeakMag_pos_index, Peak_duration, Peaks_pos_groups)
  } else {
    PeakMag_pos <- as.data.frame(matrix(NA, nrow = 1, ncol = 2))
    colnames(PeakMag_pos) <- c("ageI", "peak_mag")
  }
  
  ## Get peak-magnitude values for negative peaks
  if (length(Peaks.neg.final) > 0) {
    PeakMag_neg_val <- thresh.neg - v
    PeakMag_neg_val[PeakMag_neg_val < 0] <- 0
    PeakMag_neg_index <- which(PeakMag_neg_val > 0)
    
    Peaks_neg_groups <- split(PeakMag_neg_index,
                              cumsum(c(1, diff(PeakMag_neg_index) != 1)))
    
    PeakMag_neg <- vector(mode = "numeric", length = length(Peaks_neg_groups))
    PeakMag_neg_ages <- vector(mode = "numeric", length = length(Peaks_neg_groups))
    Peak_duration <- vector(mode = "numeric", length = length(Peaks_neg_groups))
    
    for (i in 1:length(Peaks_neg_groups)) {
      n_groups <- length(Peaks_neg_groups[[i]])
      Peak_duration[i] <- ageI[ Peaks_neg_groups[[i]] [n_groups] + 1] -
        ageI[ Peaks_neg_groups[[i]] [1]]
      PeakMag_neg[i] <- sum(c(PeakMag_neg_val[ Peaks_neg_groups[[i]] ])) *
        Peak_duration[i]
      PeakMag_neg_ages[i] <- ageI[Peaks_neg_groups[[i]] [n_groups]]
    }
    
    PeakMag_neg <- as.data.frame(cbind(PeakMag_neg_ages, PeakMag_neg))
    colnames(PeakMag_neg) <- c("ageI", "peak_mag")
    
    rm(PeakMag_neg_val, PeakMag_neg_index, Peak_duration, Peaks_neg_groups)
  } else {
    PeakMag_neg <- as.data.frame(matrix(NA, nrow = 1, ncol = 2))
    colnames(PeakMag_neg) <- c("ageI", "peak_mag")
  }
  
  
  
  ## Make summary plots ####
  
  # Set y-axis limits
  y.lim <- c(min(v, na.rm = T), max(v, na.rm = T))
  
  pdf(paste0(out.path, s.name, '_', v.name, '_Local_ThresholdSeries.pdf'))
  par(mfrow = c(2,1), oma  =  c(3,1,0,0), mar = c(1,4,0.5,0.5))
  
  plot(ageI, a[ ,2], type = "s", axes = F, ylab = a.names[2], xlab = a.names[1],
       xlim = x.lim)
  abline(h = 0)
  lines(ageI, thresh.pos.sm, col = "red")
  lines(ageI, thresh.neg.sm, col = "blue")
  points(ageI[Peaks.pos.final], rep(0.9*y.lim[2], length(Peaks.pos.final)),
         pch = 3, col = "red", lwd = 1.5)
  points(x = ageI[insig.peaks], y = rep(0.9*y.lim[2], length(insig.peaks)),
         pch = 16, col = "darkgrey", lwd = 1.5)
  points(ageI[Peaks.neg.final], rep(0.8*y.lim[2], length(Peaks.neg.final)),
         pch = 3, col = "blue", lwd = 1.5)
  axis(side = 1, labels = T, tick = T)
  axis(2)
  mtext(paste0("thresh.value  =  ", thresh.value), side = 3, las = 0, line = -0.5)
  
  plot(ageI, SNI_pos$SNI_raw, type = "p", xlim = x.lim, col = "grey",
       ylim = c(0, 1.2*max(SNI_pos$SNI_raw, na.rm = T)), axes  =  F,
       xlab  =  "Age (cal yrs BP)", ylab  =  "SNI")
  lines(ageI, SNI_pos$SNI_sm, col = "black", lwd = 2)
  abline(h  =  3, lty  =  "dashed")
  axis(side = 1, labels = T, tick = T)
  axis(2)
  
  dev.off()
  
  
  
  ## Prepare output
  out1 <- structure(list(proxy = proxy, ages.thresh = ageI,
                         thresh.value = thresh.value,
                         SNI_pos = SNI_pos, SNI_neg = SNI_neg,
                         thresh.pos = thresh.pos, thresh.neg = thresh.neg,
                         thresh.pos.sm = thresh.pos.sm, thresh.neg.sm = thresh.neg.sm,
                         peaks.pos = Peaks.pos, peaks.neg = Peaks.neg,
                         peaks.pos.insig = Peaks.pos.insig,
                         peaks.pos.ages = peaks.pos.ages, peaks.neg.ages = peaks.neg.ages,
                         PeakMag_pos = PeakMag_pos, PeakMag_neg = PeakMag_neg,
                         RI_neg = RI_neg, RI_pos = RI_pos,
                         span.sm = span.sm,
                         x.lim = x.lim))
  
  a.out <- append(series, list(out1))
  names(a.out) [5] <- "thresh"
  
  return(a.out)
  
}