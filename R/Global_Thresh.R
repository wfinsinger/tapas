#' Determine threshold values to extract signal from a detrended timeseries.
#'
#' The script determines threshold values
#' to decompose a detrended timeseries
#' into a *noise* and a *signal* component
#' using a 2-component Gaussian Mixture Model (GMM).
#' This is based on Phil Higuera's CharThreshLocal.m Matlab code.
#' It determines a positive and a negative threshold value
#' for each interpolated sample,
#' based on the distribution of values within the entire record.
#' The procedure uses a Gaussian mixture model
#' with the assumption that the noise component
#' is normally distributed around 0 (because input values were detrended!).
#' A figure is generated and saved to the \code{output} directory
#' and a list is returned with the threshold data for the analyzed proxy.
#'
#' Requires output from the \code{\link{SeriesDetrend}()} function.
#'
#' @param series The output of the \code{SeriesDetrend()} function.
#' @param proxy Set \code{proxy = "VariableName"} to select the variable
#'              for the peak-detection analysis.
#'              If the dataset includes only one variable,
#'              \code{proxy} does not need to be specified.
#' @param t.lim Restricted portion of the timeseries.
#'              With \code{t.lim = NULL} (by default),
#'              the analysis is performed with the entire timeseries.
#' @param thresh.value Threshold: the nth-percentile of the Gaussian Model
#'                     of the noise component.
#'                     Defaults to \code{thresh.value = 0.95}.
#' @param noise.gmm Specifies which of the two GMM components
#'                 should be considered as the noise component.
#'                 By default \code{noise.gmm = 1}.
#' @param smoothing.yr Width of the moving window
#'                     for computing \code{\link{SNI}}.
#'                     By default, this value is inherited
#'                     from the \code{smoothing.yr} value
#'                     set in the \code{SeriesDetrend()} function.
#' @param keep_consecutive Logical. When \code{FALSE} (by default),
#'                         consecutive peak samples exceeding the threshold
#'                         will be removed,
#'                         and only the first sample (the oldest) is retained.
#' @param minCountP Probability that two resampled counts
#'                  could arise from the same Poisson distribution
#'                  (defaults to \code{0.05}).
#'                  This is used to screen peak samples and remove
#'                  any that fails to pass the minimum-count test.
#'                  If \code{MinCountP = NULL}, the test is not performed.
#' @param MinCountP_window Width of the search window (in years)
#'                         used for the minimum-count test.
#'                         Defaults to \code{MinCountP_window = 150}.
#' @param out.dir Path to the folder where figures will be written to.
#'                Use \code{out.dir = NULL} to plot to current device instead.
#' @param plot.global_thresh Logical. If \code{TRUE}, then \code{*.pdf} files
#'                           are produced and written to \code{out.dir} folder.
#'
#' @return A list similar to \code{series} with additional data appended.
#'
#' @importFrom stats complete.cases dnorm pt qnorm
#' @importFrom graphics curve hist points
#'
#' @export
global_thresh <- function(series = NA, proxy = NULL, t.lim = NULL,
                          thresh.value = 0.95, noise.gmm = 1, smoothing.yr = NULL,
                          keep_consecutive = F,
                          minCountP = 0.05, MinCountP_window = 150,
                          out.dir = NULL, plot.global_thresh = T) {

  ## Put the user’s layout settings back in place when the function is done
  opar <- par("mfrow", "mar", "oma", "cex")
  on.exit(par(opar))

  # initial check-up of input parameters ####
  if (keep_consecutive == T & is.null(minCountP) == F) {
    stop("Fatal error: inconsistent choice of arguments. If keep_consecutive=T, set minCountP=NULL.")
  }


  # Determine path to output folder for Figures
  if (plot.global_thresh == T && !is.null(out.dir)) {
    out.path <- paste0("./", out.dir, "/")

    if (dir.exists(out.path) == F) {
      dir.create(out.path)
    }
  }

  # Extract parameters from input list (i.e. the detrended series) ####

  # Determines for which of the variables (proxy) in the dataset the analysis should be made
  #
  if (is.null(smoothing.yr) == T) {
    smoothing.yr <- series$detr$smoothing.yr
  }

  if (is.null(proxy) == T) { # if proxy = NULL, use the data in series$detr$detr
    if (dim(series$detr$detr)[2] > 2) {
      stop("Fatal error: please specify which proxy you want to use")
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


  # Determine which of the two GMM components should be considered as 'noise'
  if (noise.gmm == 1) {
    signal.gmm <- 2
  }
  if (noise.gmm == 2) {
    signal.gmm <- 1
  }


  # Create empty list where output data will be stored
  a.out <- list(thresh.value = thresh.value, yr.interp = yr.interp)


  # Create space for local variables
  v <- a[ ,2] # variable
  v.name <- colnames(a)[2]
  v.gmm <- v[which(complete.cases(v))]
  y.lim <- c(min(v.gmm), max(v.gmm))



  # Determine global threshold with 2 components ####
  #
  # Make space in empty matrices
  Peaks.pos <- matrix(data = 0, nrow = length(v), ncol = 1) # Space for positive component
  Peaks.neg <- matrix(data = 0, nrow = length(v), ncol = 1) # Space for negative component

  # Get 2 components with Gaussian mixture models for positive values
  m2 <- mclust::densityMclust(data = v.gmm, G = 2,
                              verbose = FALSE, plot = FALSE)

  if (m2$parameters$mean[1] == m2$parameters$mean[2]) {
    warning('Poor fit of Gaussian mixture model')
  }

  if (length(m2$parameters$variance$sigmasq) == 1) {
    m2.sigmasq <- m2$parameters$variance$sigmasq
  } else {
    m2.sigmasq <- m2$parameters$variance$sigmasq[noise.gmm]
  }

  ## Define global thresholds ####
  thresh.pos <- m2$parameters$mean[noise.gmm] + qnorm(p = thresh.value) * sqrt(m2.sigmasq)
  thresh.neg <- m2$parameters$mean[noise.gmm] - qnorm(p = thresh.value) * sqrt(m2.sigmasq)

  ## Get data for the 2 components
  #sig.pos <- v[which(v > thresh.pos)]
  #noise_i <- v[which(v <= thresh.pos & v >= thresh.neg)]
  #sig.neg <- v[which(v < thresh.neg)]



  ## Calculate SNI ####

  ## Get resampled values for selected variable (proxy) for selected time interval (t.lim)
  SNI_in_index <- which(series$int$series.int$age <= max(t.lim) &
                          series$int$series.int$age >= min(t.lim))

  SNI_in <- series$int$series.int[[proxy]] [SNI_in_index]

  thresh_pos_sni <- SNI_in - series$detr$detr[[proxy]] + thresh.pos
  SNI_pos <- SNI(ProxyData = cbind(ageI, SNI_in, thresh_pos_sni),
                 BandWidth = smoothing.yr)

  thresh_neg_sni <- SNI_in - series$detr$detr[[proxy]] + thresh.neg
  SNI_neg <- SNI(ProxyData = cbind(ageI, -1 * SNI_in, -1 * thresh_neg_sni),
                 BandWidth = smoothing.yr)
  rm(SNI_in, SNI_in_index)


  ## Get peaks ####

  ## If keep_consecutive == F ####
  if (keep_consecutive == F) { # if consecutive peak samples should be removed

    Peaks.pos[which(v > thresh.pos)] <- 2
    Peaks.neg[which(v < thresh.neg)] <- 2

    # Then remove consecutive peaks
    # For positive peaks
    for (i in 1:(length(Peaks.pos) - 1)) { # For each value in Peaks.pos
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
    for (i in 1:(length(Peaks.neg) - 1)) { # For each value in Peaks.neg
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
    # Get all peak samples
    Peaks.pos[which(v > thresh.pos)] <- 1
    Peaks.neg[which(v < thresh.neg)] <- 1
  }




  ## Minimum-count Analysis ####
  ## If minCountP is not NULL, screen positive peaks with minimum count test
  if (is.null(minCountP) == F) {

    # Set [yr] Years before and after a peak to look for the min. and max. value
    if (is.null(MinCountP_window) == T) {
      MinCountP_window <- round(150/yr.interp) * yr.interp
    }

    countI_index <- which(colnames(series$int$series.int) == proxy)
    countI <- series$int$series.countI[ ,countI_index]
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



  ## Gather data for plots ####

  # Use these for plotting significant peaks
  Peaks.pos.final <- which(Peaks.pos == 1)
  Peaks.neg.final <- which(Peaks.neg == 1)

  peaks.pos.ages <- ageI[which(Peaks.pos == 1)]
  peaks.neg.ages <- ageI[which(Peaks.neg == 1)]


  ## Calculate Event Return Intervals ####
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

  # Plot Gaussian Mixture Model with thresholds and classification (from densityMclust)

  ### Get parameters for the PDFs obtained from Gaussian mixture models
  muHat <- m2$parameters$mean
  sigmaHat <- sqrt(m2$parameters$variance$sigmasq)
  propN <- m2$parameters$pro

  ## Make pdf
  if (plot.global_thresh == T) {
    if (!is.null(out.dir)) {
    pdf(paste0(out.path, s.name, '_', v.name, '_Global_01_GMM_Evaluation.pdf'),
        onefile = TRUE, paper = "a4")
    }
    par(mfrow = c(2,1), mar = c(6,4,1,1), oma = c(1,1,1,1))

    # Plot classification
    mclust::plot.Mclust(m2, what = "classification",
                        xlab = paste(proxy, "(detrended)"))

    # Gather data to plot GMMs
    h <- hist(v, breaks = 50, plot = F)
    gmm_sig <- dnorm(x = h$breaks,
                     mean = muHat[signal.gmm],
                     sd = sigmaHat[signal.gmm])
    gmm_noise <- dnorm(x = h$breaks,
                       mean = muHat[noise.gmm],
                       sd = sigmaHat[noise.gmm])

    # Plot GMMs and thresholds
    plot(h, freq = F, col = "grey", border = "grey",
         xlim = c(min(h$breaks), max(h$breaks)),
         ylab = "Density", xlab = paste(proxy, "(detrended)"),
         main = "Global threshold")
    par(new = T)
    plot(h$breaks, gmm_sig,
         xlim = c(min(h$breaks), max(h$breaks)),
         ylim = c(0, max(h$density)), type = "l", col = "blue", lwd = 1.5,
         axes = F, ylab = '', xlab = '')
    par(new = T)
    plot(h$breaks, gmm_noise,
         xlim = c(min(h$breaks), max(h$breaks)),
         ylim = c(0, max(h$density)), type = "l", col = "orange", lwd = 1.5,
         axes = F, ylab = '', xlab = '')
    par(new = T)
    lines(c(thresh.pos, thresh.pos), y = c(0, max(h$density)), type = "l",
          col = "red", lwd = 1.5)
    lines(c(thresh.neg, thresh.neg), y = c(0, max(h$density)), type = "l",
          col = "red", lwd = 1.5)
    mtext(paste0("thresh.value  =  ", thresh.value), side = 3, las = 0, line = -1)
    mtext(paste0("propN = ", round(propN[1], 2), ", ", round(propN[2], 2)), side = 3,
          las = 0, line = -2)
    if (!is.null(out.dir)) {
    dev.off()
    }


    ## Plot series with trend + threshold + significant peaks
    if (!is.null(out.dir)) {
    pdf(paste0(out.path, s.name, '_', v.name, '_Global_02_ThresholdSeries.pdf'))
    }
    par(mfrow = c(2,1), oma = c(3,1,0,0), mar = c(1,4,0.5,0.5))
    plot(ageI, a[ ,2], type = "s", axes = F, ylab = a.names[2], xlab = "",
         xlim = x.lim)
    abline(h = 0)
    abline(h = thresh.pos, col = "red")
    abline(h = thresh.neg, col = "blue")
    points(ageI[Peaks.pos.final], rep(0.9*y.lim[2], length(Peaks.pos.final)),
           pch = 3, col = "red", lwd = 1.5)
    points(ageI[insig.peaks], rep(0.9*y.lim[2], length(insig.peaks)),
           pch = 16, col = "darkgrey", lwd = 1.5)
    points(ageI[Peaks.neg.final], rep(0.8*y.lim[2], length(Peaks.neg.final)),
           pch = 3, col = "blue", lwd = 1.5)
    axis(side = 1, labels = F, tick = T)
    axis(2)
    mtext(paste0("thresh.value = ", thresh.value), side = 3, las = 0, line = -0.5)

    plot(ageI, SNI_pos$SNI_raw, type = "p", xlim = x.lim, col = "grey",
         ylim = c(0, 1.2*max(SNI_pos$SNI_raw, na.rm = T)), axes  =  F,
         xlab  =  "Age (cal yrs BP)", ylab  =  "SNI")
    lines(ageI, SNI_pos$SNI_sm, col = "black", lwd = 2)
    abline(h  =  3, lty  =  "dashed")
    axis(side = 1, labels = T, tick = T)
    axis(2)
    if (!is.null(out.dir)) {
    dev.off()
    }
  }

  ## Gather data for output
  out1 <- structure(list(proxy = proxy, ages.thresh = ageI,
                         thresh.value = thresh.value,
                         SNI_pos = SNI_pos, SNI_neg = SNI_neg,
                         thresh.pos = thresh.pos, thresh.neg = thresh.neg,
                         thresh.pos.sm = thresh.pos, thresh.neg.sm = thresh.neg,
                         peaks.pos = Peaks.pos, peaks.neg = Peaks.neg,
                         peaks.pos.insig = Peaks.pos.insig,
                         peaks.pos.ages = peaks.pos.ages, peaks.neg.ages = peaks.neg.ages,
                         PeakMag_pos = PeakMag_pos, PeakMag_neg = PeakMag_neg,
                         RI_neg = RI_neg, RI_pos = RI_pos,
                         x.lim = x.lim))

  ## Merge output into input list
  a.out <- append(series, list(out1))
  names(a.out) [5] <- "thresh"

  ## Return final output
  return(a.out)
}
