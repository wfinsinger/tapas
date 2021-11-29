#  ***************************************************************************
#   Local_Thresh.R
#  ---------------------
#   Date                 : July 2020
#   Copyright            : (C) 2020 by Walter Finsinger
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
#   proxy   ->  Set proxy="VariableName" to determine threshold for one single variable.
#     In this case, a Figure will be generated and saved to the "output" directory on the
#     hard drive and a list will be returned with the threshold data for the analysed proxy.
#     With the beta-version default option (proxy=NULL), the analysis uses all variables
#     that are in the output from the SeriesDetrend() function;
#     Figures for each variable will be generated and saved to the "output" directory, and
#     a list will be returned with the threshold data for the last proxy only.
#   
#   t.lim   ->  allows defining a portion of the time series.
#     Should be t.lim=c(older age limit, younger age limit).
#     With t.lim=NULL (as per default) the analysis will be performed using the entire timeseries.
#   
#   thresh.value  ->  Determines the location of the threshold as the percentile of the
#     Gaussian Model of the noise component
#   
#   noise.gmm     =>  Determines which of the two GMM components should be considered as
#     the noise component. By default noise.gmm=1.
#   
#   thresh.yr     =>  determines the length of the window width from which
#                       values are selected to determine the local threshold (if gm.local=T)
#   
#   smoothing.yr  =>  determines the smoothing-window width to smooth the threshold values
#   
#   span.sm       =>  determines the smoothing of the threshold values.
#                     If span.sm=NULL, the span is calculated based on smoothing.yr
#
#  ***************************************************************************


local_thresh <- function(series=NA, proxy=NULL, t.lim=NULL, thresh.yr=1000,
                         thresh.value=0.95, noise.gmm=1, smoothing.yr=500, span.sm=NULL,
                         gm.local=T, keep_consecutive=F, out.dir="Figures") {
  
  
  series = co.detr1
  proxy = "charAR"
  t.lim = NULL
  thresh.yr = 500
  thresh.value = 0.99
  noise.gmm = 1
  smoothing.yr = 500
  span.sm = NULL
  #span.sm = 0.2
  gm.local = T
  keep_consecutive = F
  out.dir = "Figures"
  
  
  require(mclust)
  
  # Determine path to output folder for Figures
  out.path <- paste0("./", out.dir, "/")
  
  if (dir.exists(out.path) == F) {
    dir.create(out.path)
  }
  
  
  # Extract parameters from input list (detrended series) ####
  
  # Determines if analysis of all proxies or of one single proxy
  if (is.null(proxy) == T) { # if proxy = NULL, use the data in series$detr$detr
    if (dim(series$detr$detr)[2] > 2) {
      print('Fatal error: please specify which proxy you want to use')
      return()
    } else {
      a <- series$detr$detr
    }
  }
  
  if (is.null(proxy) == F) { # if we want to select a specific variable
    a <- data.frame(series$detr$detr$age)
    colnames(a) <- "age"
    a <- cbind(a, series$detr$detr[[proxy]])
    colnames(a)[2] <- proxy
  }
  
  # Further extract parameters from input dataset
  a.names <- colnames(a)
  s.name <- series$detr$series.name
  yr.interp <- series$int$yr.interp
  
  # Determine whether to limit the analysis to a portion of the record
  if (is.null(t.lim) == T) {
    ageI <- a$age
  }
  if (is.null(t.lim) == F) {
    a <- a[which(a$age <= max(t.lim) & a$age >= min(t.lim)), ]
    ageI <- a$age[which(a$age <= max(t.lim) & a$age >= min(t.lim))]
  }
  x.lim <- c(max(ageI), min(ageI))
  
  
  if (noise.gmm == 1) {
    signal.gmm <- 2
  }
  if (noise.gmm == 2) {
    signal.gmm <- 1
  }
  
  # Determine the proportion of datapoints used to smooth local thresholds with loess()
  if (gm.local == T) {
    n.smooth <- round(smoothing.yr/yr.interp)
    if (is.null(span.sm)) {
      span.sm <- n.smooth/dim(a)[1]
    }
  }
  
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
    
    
    ## ESTIMATE LOCAL NOISE DISTRIBUTION
    # Estimate noise distribution with Guassian mixture model
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
      sigmaHat[i, noise.gmm] <- sqrt(m$parameters$variance$sigmasq[noise.gmm])
      sigmaHat[i, signal.gmm] <- sqrt(m$parameters$variance$sigmasq[signal.gmm])
      propN[i, ] <- m$parameters$pro
      
      # Check
      if (muHat[i, noise.gmm] == muHat[i, signal.gmm]) {
        warning('Poor fit of Gaussian mixture model')
      }
      
      ## Define local threshold
      if (!is.na(m$parameters$variance$sigmasq[noise.gmm]) == T) {
        thresh.pos[i] <- m$parameters$mean[noise.gmm] + qnorm(p = thresh.value) *
          sqrt(m$parameters$variance$sigmasq)[noise.gmm]
        thresh.neg[i] <- m$parameters$mean[noise.gmm] - qnorm(p = thresh.value) *
          sqrt(m$parameters$variance$sigmasq)[noise.gmm]
        
        if (m$parameters$mean[noise.gmm] < 0) {
          thresh.pos[i] <- 0 + qnorm(p = thresh.value) *
            sqrt(m$parameters$variance$sigmasq)[noise.gmm]
        }
        
        #sig_i.pos <- X[which(X > thresh.pos[i])]
        #noise_i <- X[which(X <= thresh.pos[i] & X >= thresh.neg[i])]
        #sig_i.neg <- X[which(X < thresh.neg[i])]
        
        
      }
      
      if (is.na(m$parameters$variance$sigmasq[noise.gmm]) == T) {
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
      pdf2 <- curve(dnorm(x = x, mean = muHat[i, signal.gmm], sd = sigmaHat[i, signal.gmm]),
                    from = min(h$breaks), to = max(h$breaks),
                    ylim = c(0, max(dens)), type = "l", col = "blue", lwd = 1.5,
                    axes = F, ylab = '', xlab = '')
      par(new = T)
      pdf1 <- curve(dnorm(x = x, mean = muHat[i, noise.gmm], sd = sigmaHat[i, noise.gmm]),
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
  
  
  ## Calculate SNI:
  if (is.null(proxy)) {
    SNI_pos <- SNI(ProxyData = cbind(series$int$series.int$age,
                                     series$int$series.int[ ,2],
                                     thresh.pos),
                   BandWidth = smoothing.yr)
    SNI_neg <- SNI(ProxyData = cbind(series$int$series.int$age,
                                     -1 * series$int$series.int[ ,2],
                                     -1 * thresh.neg),
                   BandWidth = smoothing.yr)
  } else {
    SNI_pos <- SNI(ProxyData = cbind(series$int$series.int$age,
                                     series$int$series.int[[proxy]],
                                     thresh.pos),
                   BandWidth = smoothing.yr)
    SNI_neg <- SNI(ProxyData = cbind(series$int$series.int$age,
                                     -1 * series$int$series.int[[proxy]],
                                     -1 * thresh.neg),
                   BandWidth = smoothing.yr)
  }
  
  
  # Smooth local thresholds
  ## Smooth thresholds and SNI with Lowess smoother
  # thresh.pos.sm <- stats::loess(thresh.pos ~ ageI, data = data.frame(ageI, thresh.pos),
  #                        span = span.sm, family = "gaussian")$fitted
  # thresh.neg.sm <- stats::loess(thresh.neg ~ ageI, data = data.frame(ageI, thresh.neg),
  #                        span = span.sm, family = "gaussian")$fitted
  
  thresh.pos.sm <- stats::lowess(ageI, thresh.pos, f = span.sm, iter = 0)$y
  thresh.neg.sm <- stats::lowess(ageI, thresh.neg, f = span.sm, iter = 0)$y
  
  
  
  ## Get peaks
  if (keep_consecutive == F) { # if consecutive peaks should be removed
    
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
  
  
  ## Plot series with trend + threshold + peaks
  Peaks.pos.plot <- which(Peaks.pos > 0)
  Peaks.neg.plot <- which(Peaks.neg > 0)
  
  peaks.pos.ages <- ageI[which(Peaks.pos == 1)]
  peaks.neg.ages <- ageI[which(Peaks.neg == 1)]
  
  
  ## Calculate Event Return Intervals
  RI_neg <- c(diff(peaks.neg.ages), NA)
  RI_pos <- c(diff(peaks.pos.ages), NA)
  
  
  # Set y-axis limits
  y.lim <- c(min(v, na.rm = T), max(v, na.rm = T))
  
  pdf(paste0(out.path, s.name, '_', v.name, '_Local_ThresholdSeries.pdf'))
  par(mfrow = c(2,1), oma  =  c(3,1,0,0), mar = c(1,4,0.5,0.5))
  
  plot(ageI, a[ ,2], type = "s", axes = F, ylab = a.names[2], xlab = a.names[1],
       xlim = x.lim)
  abline(h = 0)
  lines(ageI, thresh.pos.sm, col = "red")
  lines(ageI, thresh.neg.sm, col = "blue")
  points(ageI[Peaks.pos.plot], rep(0.9*y.lim[2], length(Peaks.pos.plot)),
         pch = 3, col = "red", lwd = 1.5)
  points(ageI[Peaks.neg.plot], rep(0.8*y.lim[2], length(Peaks.neg.plot)),
         pch = 3, col = "blue", lwd = 1.5)
  axis(side = 1, labels = T, tick = T)
  axis(2)
  mtext(paste0("thresh.value  =  ", thresh.value), side = 3, las = 0, line = -1)
  
  plot(ageI, SNI_pos$SNI_raw, type = "p", xlim = x.lim, col = "grey",
       ylim = c(0, 1.2*max(SNI_pos$SNI_raw, na.rm = T)), axes  =  F,
       xlab  =  "Age", ylab  =  "SNI")
  lines(ageI, SNI_pos$SNI_sm, col = "black", lwd = 2)
  abline(h  =  3, lty  =  "dashed")
  axis(side = 1, labels = T, tick = T)
  axis(2)
  
  dev.off()
  
  
  
  ## Prepare output
  out1 <- structure(list(proxy = proxy, ages.thresh = ageI,
                         gm.local = gm.local, thresh.value = thresh.value,
                         SNI_pos = SNI_pos, SNI_neg = SNI_neg,
                         thresh.pos = thresh.pos, thresh.neg = thresh.neg,
                         thresh.pos.sm = thresh.pos.sm, thresh.neg.sm = thresh.neg.sm,
                         peaks.pos = Peaks.pos, peaks.neg = Peaks.neg,
                         peaks.pos.ages = peaks.pos.ages, peaks.neg.ages = peaks.neg.ages,
                         RI_neg = RI_neg, RI_pos = RI_pos,
                         x.lim = x.lim))
  
  a.out <- append(series, list(out1))
  names(a.out) [4] <- "thresh"
  
  return(a.out)
  
}