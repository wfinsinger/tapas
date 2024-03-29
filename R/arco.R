#' Screen peaks when input data is charcoal areas
#'
#' @description
#'  Screens peaks when input data is charcoal-area records. Modified from
#'  the ARCO v1 function, available at url{https://github.com/wfinsinger/ARCO}.
#'
#'
#' @param Seedle.file A data frame with charcoal-particle areas. Should have as
#'  many rows as the number of charcoal particles and two columns in the
#'  following order:
#'    - Column 1: Top depth of samples;
#'    - Column 2: Charcoal-particle areas.
#' @param Smpl.file A data frame with as many rows as the number of samples and
#'   seven columns in the following order: \code{CmTop, CmBot, AgeTop, AgeBot,
#'   volume, charcoal counts, charcoal areas}.
#' @param FireA.file Output from peak-detection analysis based on charcoal areas.
#'    Generate it using first the \code{tapas::peak_detection} function (set
#'    argument \code{min_CountP = NULL}), and thereafter the
#'    \code{tapas::tapas_export()} function.
#' @param FireC.file Output from peak-detection analysis based on charcoal counts.
#'    Generate it using first the \code{tapas::peak_detection} function (set
#'    argument \code{0 < min_CountP < 1}), and thereafter the
#'    \code{tapas::tapas_export()} function.
#' @param n.boot Number of bootstrap samples generated by the function to
#'  obtain a distribution of simulated charcoal-areas. By default,
#'  \code{n.boot = 10000}.
#' @param thresh.prob The pth percentile threshold used to separate significant
#'   charcoal-area peaks. By default \code{thresh.prob = 0.95}.
#' @param win.width The temporal span of the window from which bootstrap samples
#'   are generated. For each peak to be screened, particles are randomly drawn
#'   (with replacement) from all samples within a focal window, which is
#'   centered on the peak and has a full span of win.width.
#'   By default \code{win.width = 1000}.
#' @param plotit Logical. If \code{plotit = FALSE} (default), plots are not
#'   sent to the device.
#' @param breakage Logical. If \code{breakage = FALSE}, plots also
#'   C#/CA-ratios in one of the diagnostic plots. By default
#'   \code{breakage = FALSE}.
#' @param storedat Logical. If \code{storedat = TRUE}, the function returns
#'   also other output data related to charcoal-particle areas.
#'
#' @details
#' This screening procedure is specific for data sets comprising charcoal
#' numbers (counts) and areas. It screens the charcoal-area estimates with
#' respect to the number and size of charcoal particles.
#'
#' The method begins with a charcoal-area data set analysed by existing
#' methods to identify peaks representing fire episodes, e.g. the
#' \code{peak_detection()} function with the argument
#' \code{min_CountP = NULL}).
#' To screen these peaks, the method uses bootstrap resampling of
#' charcoal-particle areas observed in a user-defined subsection of the data
#' set around each peak (by default \code{win.width = 1000} to obtain the range
#' of likely charcoal areas for different counts. Peaks with total area within
#' the likely range of bootstrapped samples (e.g. p > 0.05) are flagged as
#' potentially unreliable, whereas peak samples with total area
#' significantly greater than expected by chance are deemed robust indicators
#' of past fire episodes.
#'
#' @returns
#' By default (as with \code{storedat = FALSE}), the function returns a
#' list with the same structure as the input list (FireA.file).
#'
#' @references
#' Finsinger, W., R. Kelly, J. Fevre, and E.K. Magyari. 2014. A guide to
#'  screening charcoal peaks in macrocharcoal-area records for fire episode
#'  reconstructions. The Holocene, 24:1002–1008. doi: 10.1177/0959683614534737
#'
#' Higuera, P.E., L.B. Brubaker, P.M. Anderson, F.S. Hu, and T.A. Brown. 2009.
#'  Vegetation mediated the impacts of postglacial climate change on fire
#'  regimes in the south-central Brooks Range, Alaska. Ecological Monographs,
#'  79:201–219.
#'
#' @author Ryan Kelly
#' @author Walter Finsinger
#'
#' @importFrom stats quantile lm
#' @importFrom dplyr select all_of count relocate
#' @importFrom graphics boxplot legend
#' @importFrom tidyselect last_col
#' @importFrom grDevices grey
#' @importFrom rlang .data
#'
#' @export

arco <- function(Seedle.file, Smpl.file, FireA.file, FireC.file,
                     n.boot = 10000, thresh.prob = 0.95, win.width = 1000,
                     plotit = FALSE, breakage = FALSE, storedat = FALSE) {


  # -------------------- SETUP -------------------- #

  ## Put the user’s layout settings back in place when the function is done
  opar <- par("mfrow", "mar", "mai", "oma", "cex", "cex.lab")
  on.exit(par(opar))

  # ----- LOAD files with Charcoal per sample and of Seedles
  # file with charcoal-particle areas
  Seedle <- Seedle.file
  colnames(Seedle) <- c("Depth", "SdlArea")

  # file with charcoal counts and charcoal areas
  positions <- c(1,3,6,7)
  Smpl <- Smpl.file %>% dplyr::select(all_of(positions))
  colnames(Smpl) <- c("Depth", "Age_calBP", "SmplCount", "SmplArea")

  # LOAD files with Fire history derived from tapas
  ca_temp <- tapas::tapas_export(FireA.file)
  CA.dat <- ca_temp[complete.cases(ca_temp$cm_i), ]

  # tapas output, with Charcoal counts WITH pMinCount
  cc_temp <- tapas::tapas_export(FireC.file)
  CC.dat <- cc_temp[complete.cases(cc_temp$cm_i), ]

  # Define minimum-area bootstrapping function
  b.stat <- function(seedle.areas, sample.size, n.boot, thresh.prob) {
    rj.resample <- matrix(sample(seedle.areas, sample.size*n.boot,
                                 replace = T), nrow = n.boot)
    rj.sum <- apply(rj.resample, 1, sum)
    return(stats::quantile(rj.sum, probs = thresh.prob))
  }

  # Add Seedle counts to the Seedle data.frame
  if (min(Seedle$SdlArea) == 0) {      # thus if Seedle file contains zeros
    S <- Seedle[which(Seedle$SdlArea > 0), ]
  } else {
    S <- Seedle
  }

  sdl_counts <- S %>% dplyr::count(.data$Depth, name = "SdlCounts")
  sdl <- merge(Seedle.file, sdl_counts, by = "Depth", all.x = TRUE)
  sdl$SdlCounts[is.na(sdl$SdlCounts)] <- 0
  Seedle <- sdl %>% dplyr::relocate(.data$SdlArea, .after = last_col())

  rm(S, sdl_counts, sdl)

  # Add ages to each Seedle (used later for local bootstrapping)
  Age <- Smpl[ ,-c(3,4)]
  Seedle <- merge(Seedle, Age, by.x = "Depth", all.y = FALSE)


  # Extract data to obtain identified peaks and samples above threshold
  CA.age    <- CA.dat$age_top_i # [ ,2]
  CA.cpeak  <- CA.dat$bkg_trend # [ ,8]
  CA.thresh <- CA.dat$thresh_final_pos_sm # [,12]
  CA.peaks  <- CA.dat$peaks_pos_sig # [,19]
  CA.res    <- mean(diff(CA.age))
  CC.age    <- CC.dat$age_top_i # [,2]
  CC.cpeak  <- CC.dat$bkg_trend # [,8]
  CC.thresh <- CC.dat$thresh_final_pos_sm # [,12]


  ## SCREENING --------------------

  # Find all samples contributing to identified peaks
  peak.ind <- which(CA.peaks == 1)
  n.peak <- length(peak.ind)

  # Create some variables
  overthresh.ind = peaksamples.ind = peaksamples.count =
    peaksamples.area = peaksamples.age = list()
  overthresh.interval = matrix(NA, nrow = n.peak, ncol = 2)

  # Loop over all identified peaks to find all samples contributing to each peak.
  # (Because when multiple samples in a row exceed the threshold, only the
  # lowermost one is marked as a peak)
  for (i in 1:n.peak) {
    # Start with the sample identified as the peak
    overthresh.ind[[i]] <- peak.ind[i]

    # If the first peak is in the first sample, don't look for other
    # samples belonging to the peak
    if (peak.ind[i] == 1) {
      overthresh.ind[[1]] <- peak.ind[1]
    }
    if (peak.ind[i] == 2) {
      ind.cur <- peak.ind[i] - 1    # start with the next sample above...
      if (CA.cpeak[ind.cur] > CA.thresh[ind.cur]) { # if it's also above the threshold...
        overthresh.ind[[i]] <- c(overthresh.ind[[i]], ind.cur)  # add it to 'overthresh.ind'...
      } else {
        overthresh.ind[[i]] <- peak.ind[i]  # add it to
      }
    }
    if (peak.ind[i] > 2) {
      # Else Look for other samples belonging to this peak
      ind.cur <- peak.ind[i] - 1   # start with the next sample above...
      while (CA.cpeak[ind.cur] > CA.thresh[ind.cur]) {  # if it's also above the threshold...
        overthresh.ind[[i]] <- c(overthresh.ind[[i]], ind.cur)  # add it to 'overthresh.ind'...
        ind.cur <- ind.cur - 1  # and get ready to check the next sample above.
      }
    }


    # For each peak, find the full age interval spanned by all samples
    # exceeding the threshold.
    overthresh.interval[i,] <- c(min(CA.age[overthresh.ind[[i]]]),
                                 max(CA.age[overthresh.ind[[i]]]) + CA.res)

    # Now go to the original data (i.e. not the CharAnalysis-interpolated data)
    # and find all samples that contribute to each peak identified in CharAnalysis
    peaksamples.ind[[i]] <- which(Smpl$Age_calBP >= overthresh.interval[i,1] &
                                    Smpl$Age_calBP <= overthresh.interval[i,2])

    # Add the sample that has its top depth just above the top of the peak
    # interval, since it partially contributes to the peak as well.
    peaksamples.ind[[i]] <- union(
      peaksamples.ind[[i]], which.min(
        sapply(overthresh.interval[i,1] - Smpl$Age_calBP,
               function(x) ifelse(x < 0, NA, x))))

    # Finally, extract the original count, area, and age data for all of these samples.
    peaksamples.count[[i]] <- Smpl$SmplCount[ peaksamples.ind[[i]] ]
    peaksamples.area[[i]] <- Smpl$SmplArea[ peaksamples.ind[[i]] ]
    peaksamples.age[[i]] <- Smpl$Age_calBP[ peaksamples.ind[[i]] ]
  }


  # ----- Separate the seedle areas for peak vs. non-peak samples
  # First, find the age intervals not spanned by peaks,
  # i.e. anything not covered by overthresh.interval.
  underthresh.interval <- matrix(NA, nrow = n.peak + 1, ncol = 2)
  underthresh.interval[1, 1] <- min(CA.age)
  underthresh.interval[2:(n.peak + 1), 1] <- overthresh.interval[ ,2]
  underthresh.interval[1:n.peak, 2] <- overthresh.interval[ ,1]
  underthresh.interval[(n.peak + 1), 2] <- max(CA.age) + CA.res

  # Then loop over each non-peak interval to find all samples in the original
  # data that did not contribute to a CharAnalysis-identified peak.
  nonpeaksamples.ind <- list()
  for (i in 1:(n.peak + 1)) {
    nonpeaksamples.ind[[i]] <- which(
      Smpl$Age_calBP >= underthresh.interval[i,1] &
        Smpl$Age_calBP <= underthresh.interval[i,2])

    # Add the sample that has its top depth just above the top of the
    # non.peak interval, since it partially contributes to the non.peak as well.
    nonpeaksamples.ind[[i]] <- union(
      nonpeaksamples.ind[[i]], which.min(sapply(
        underthresh.interval[i,1] - Smpl$Age_calBP,
        function(x) ifelse(x < 0, NA, x))))
  }

  # Finally, get the seedle data associated with peak and non-peak samples.
  # Data in 'Seedle' table is linked to samples in the 'Smpl' table by depth.
  nonpeak.Seedle <- Seedle[
    sapply(Seedle$Depth, function(x) {x %in% Smpl$Depth[
      unlist(nonpeaksamples.ind)]}), ]
  # peak.Seedle <- Seedle[
  #   sapply(Seedle$Depth, function(x) {x %in% Smpl$Depth[
  #     unlist(peaksamples.ind)]}), ]


  # ----- Perform screening
  # Create space for output
  #
  # # begin with the assumption that no peaks pass the screening
  screen.pass <- rep(0, n.peak)
  area.thresh <- rep(NA, nrow(Smpl))

  # Loop over each CharAnalysis-identified peak...
  for (i in 1:n.peak) {
    # ...which contains one or more samples from the original data
    # (loop over each of those too)...
    for (j in 1:length(peaksamples.ind[[i]])) {
      # Find all the non-peak seedle data that is within the specified
      # temporal window surrounding the peak sample being analyzed on this iteration.
      seedles.j.ind <- which(
        nonpeak.Seedle$Age_calBP >= (peaksamples.age[[i]][j] - win.width/2) &
          nonpeak.Seedle$Age_calBP <= (peaksamples.age[[i]][j] + win.width/2)
      )

      # Obtain the seedle areas for the identified non-peak samples in the window
      seedles.j <- Seedle$SdlArea[ seedles.j.ind ]

      # Bootstrap an area threshold from those areas
      area.thresh[peaksamples.ind[[i]][j]] <- b.stat(seedles.j,
                                                     peaksamples.count[[i]][j],
                                                     n.boot,
                                                     thresh.prob)

      # If the peak sample being analyzed exceeds the bootstrap threshold value,
      # then increment the value of 'screen.pass'
      if (peaksamples.area[[i]][j] > area.thresh[peaksamples.ind[[i]][j]]) {
        screen.pass[i] <- screen.pass[i] + 1
      }
    }
  }

  # If at least one of the original samples contributing to a CharAnalysis peak
  # is beyond its bootstrapped threshold, then that CharAnalysis peak passed the screening.
  peak.ind.screened <- peak.ind[which(screen.pass > 0)]
  peak.ind.notpass <- peak.ind[which(screen.pass < 1)]


  # ----- Prepare output files
  # Results related to area screening
  comp.perc <- cbind(Smpl, area.thresh)
  comp.perc$differ <- comp.perc$SmplArea - comp.perc$area.thresh
  comp.perc$val <- as.numeric( comp.perc$differ > 0 ) # as.numeric converts T/F to 1/0
  comp.perc$frag <- comp.perc$SmplCount/comp.perc$SmplArea


  ### New output file that reflects peak-area screening ###

  # First change "peaks Final", "peaks Insig.", and "peak Mag"
  CA.dat.out <- CA.dat       # Start with the original CA data
  CA.dat.out[peak.ind.notpass, 13] <- 0  # Remove screened peaks from "peaks Final"
  CA.dat.out[peak.ind.notpass, 12] <- 1  # Add them to "peaks Insig."


  ## Extract Seedles of Screened peak samples ####
  ## (called 'val' samples)
  valSmpl.perc <- comp.perc[which(comp.perc$val == 1), ]
  valSdl.perc <- merge(Seedle, valSmpl.perc[ ,-2],
                       by.x = "Depth", all.x = FALSE)

  # To compare between screened and unscreened samples #####
  valAllSdl.perc <- merge(Seedle, comp.perc[ ,-2],
                          by.x = "Depth", all.y = TRUE)

  # Add fragmentation in Smpl data frame
  Smpl$frag <- Smpl$SmplCount/Smpl$SmplArea
  is.na(Smpl) <- sapply(Smpl, is.infinite)


  ## Calculates FRIs -----------
  ## for CA, CAscreen, CC, and CCscreen fire-episode histories
  FRI.CA <- CA.dat[CA.dat[ ,13] == 1,1:2]
  FRI.CA$FRI <- c(diff(FRI.CA[ ,2]), NA)

  FRI.CAs <- CA.dat[CA.dat.out[ ,13] == 1,1:2]
  FRI.CAs$FRI <- c(diff(FRI.CAs[ ,2]), NA)

  FRI.CC <- CC.dat[ (CC.dat[,13] == 1 | CC.dat[ ,12] == 1), 1:2]
  FRI.CC$FRI <- c(diff(FRI.CC[ ,2]), NA)

  FRI.CCs <- CC.dat[CC.dat[,13] == 1, 1:2]
  FRI.CCs$FRI <- c(diff(FRI.CCs[ ,2]), NA)


  ## PLOTS --------------------

  if (plotit == TRUE) {
    ## Plots Seedle_area, C#/CA-ratio, and Sample Area vs Char_counts ####
    par(mfrow = c(3,1), cex = 1, mar = c(0.5, 6, 0.5, 1), oma = c(5,1,1,1))
    plot(Seedle$SdlCounts, Seedle$SdlArea,
         ylab = expression(Particle~area~mm^{2}),
         pch = 20, xlim = (range(Smpl$SmplCount)), xaxt = "n")
    plot(Smpl$SmplCount, Smpl$SmplArea,
         ylab = expression(Charcoal~area~(C[A]~mm^{2})), pch = 20,
         xlim = (range(Smpl$SmplCount)), xaxt = "n")
    abline(stats::lm(Smpl$SmplArea ~ Smpl$SmplCount))
    plot(Smpl$SmplCount, Smpl$SmplCount/Smpl$SmplArea,
         ylab = expression(C["#"]/C[A] - ratio), pch = 20,
         xlim = (range(Smpl$SmplCount)))
    mtext(expression(Pieces~sample^{-1}),
          side = 1, outer = TRUE, line = 3, cex = 1.5)


    ## Boxplot of seedle areas for validated peak samples ####
    par(mfrow = c(1,1), cex.lab = 1.2, mar = c(0.5,6,0.5,5), oma = c(5,1,1,1))
    graphics::boxplot(SdlArea ~ Age_calBP, data = valSdl.perc, las = 2,
                      varwidth = TRUE, notch = FALSE,
                      ylab = expression(paste("Particle area (mm"^"2",")")))
    mtext("Age (cal yrs BP)", side = 1, outer = TRUE, line = 3, cex = 1.2)


    # Boxplot comparison screened and unscreened samples #####
    all_sdl_val <- valAllSdl.perc %>% dplyr::filter(.data$SdlCounts > 0)
    all_sdl_val$val[is.na(all_sdl_val$val)] <- 0

    par(mai = c(1,1,0.5,0.5), cex.axis = 1.2, cex.lab = 1.3)
    graphics::boxplot(SdlArea ~ val, data = all_sdl_val,
                      names = c("non signif. peaks", "signif. peaks"),
                      xlab = "",
                      ylab = expression(paste("Particle area (mm"^"2",")")),
                      notch = TRUE, varwidth = TRUE)


    # Main Diagnostic plot ####
    # (Screened CA and CC peaks that failed screening test = grey dots)
    if (breakage == TRUE) {
      par(mfrow = c(3,1), mar = c(0.5, 5, 0.5, 5), oma = c(5,1,1,1))
      y.lim <- c(min(CA.cpeak), 1.2*max(CA.cpeak))
      plot(0,0, type = 'n', xlim = rev(range(CA.age)), ylim = y.lim,
           ylab = expression(CHAR[A]*~"residuals"), xaxt = "n", cex = 1)
      for (i in 1:n.peak) {
        polygon(c(overthresh.interval[i,], rev(overthresh.interval[i,])),
                rep(y.lim, each = 2), col = "mistyrose", border = NA)
      }
      for (i in 1:(n.peak + 1)) {
        polygon(c(underthresh.interval[i,], rev(underthresh.interval[i,])),
                rep(y.lim, each = 2), col = "lightcyan", border = NA)
      }

      lines(CA.age, CA.cpeak, type = 's', col = grDevices::grey(0.5))
      abline(h = 0, col = grDevices::grey(0.5))
      lines(CA.age, CA.thresh, col = 2)

      points(CA.age[peak.ind.notpass] + CA.res/2,
             rep(0.9*y.lim[2], length(peak.ind.notpass)),
             col = "gray", pch = 16, cex = 1)
      points(CA.age[peak.ind.screened] + CA.res/2,
             rep(0.9*y.lim[2], length(peak.ind.screened)),
             col = 2, pch = 3, lwd = 2)


      y.lim <- c(min(Smpl$SmplCount), 1.2*max(Smpl$SmplCount))
      plot(0, 0, type = 'n', xlim = rev(range(CA.age)),
           ylim = y.lim, ylab = expression(Char.~counts~(pieces*~sample^{"-1"})),
           xaxt = "n", cex = 1)
      for (i in 1:nrow(overthresh.interval)) {
        polygon(c(overthresh.interval[i,], rev(overthresh.interval[i,])),
                rep(y.lim, each = 2), col = "mistyrose", border = NA)
      }
      for (i in 1:(n.peak + 1)) {
        polygon(c(underthresh.interval[i,], rev(underthresh.interval[i,])),
                rep(y.lim, each = 2), col = "lightcyan", border = NA)
      }

      lines(Smpl$Age_calBP, Smpl$SmplCount, type = 's',
            col = grDevices::grey(0.5))
      abline(h = 0, col = grDevices::grey(0.5))

      points(CA.age[peak.ind.notpass] + CA.res/2,
             rep(0.9*y.lim[2], length(peak.ind.notpass)),
             col = "gray", pch = 16, cex = 1)
      points(CA.age[peak.ind.screened] + CA.res/2,
             rep(0.9*y.lim[2], length(peak.ind.screened)),
             col = 2, pch = 3, lwd = 2)

      y.lim <- c(min(CC.cpeak) + 0.01, 1.2*max(CC.cpeak))
      plot(0, 0, type = 'n', xlim = rev(range(CA.age)), ylim = y.lim,
           ylab = expression(CHAR["#"]*~"residuals"), cex = 1)

      lines(CC.age, CC.cpeak, type = 's', col = grDevices::grey(0.5))
      abline(h = 0, col = grDevices::grey(0.5))
      lines(CC.age,CC.thresh, col = 2)

      ind <- which(CC.dat[ ,12] == 1)
      points(CC.age[ind], rep(0.7*y.lim[2], length(ind)),
             col = "gray", pch = 16, cex = 1)
      ind = which(CC.dat[ ,13] == 1)
      points(CC.age[ind], rep(0.7*y.lim[2], length(ind)),
             col = "red", pch = 3, lwd = 2)

      mtext(("Age (cal yrs BP)"), side = 1, line = 2.5, cex = 1)
    } else {
      par(mfrow = c(4,1), mar = c(0.5, 5, 0.5, 5), oma = c(5,1,1,1))
      y.lim = c(min(CA.cpeak), 1.2*max(CA.cpeak))
      plot(0, 0, type = 'n', xlim = rev(range(CA.age)), ylim = y.lim,
           ylab = expression(CHAR[A]*~"residuals"), xaxt = "n", cex = 1)
      for (i in 1:n.peak) {
        polygon(c(overthresh.interval[i,], rev(overthresh.interval[i,])),
                rep(y.lim, each = 2), col = "mistyrose", border = NA)
      }
      for (i in 1:(n.peak + 1)) {
        polygon(c(underthresh.interval[i,],rev(underthresh.interval[i,])),
                rep(y.lim, each = 2), col = "lightcyan", border = NA)
      }
      lines(CA.age, CA.cpeak, type = 's', col = grDevices::grey(0.5))
      abline(h = 0, col = grDevices::grey(0.5))
      lines(CA.age, CA.thresh, col = 2)

      points(CA.age[peak.ind.screened] + CA.res/2,
             rep(0.9*y.lim[2], length(peak.ind.screened)),
             col = "red", pch = 3, lwd = 2)
      points(CA.age[peak.ind.notpass] + CA.res/2,
             rep(0.9*y.lim[2], length(peak.ind.notpass)),
             col = "gray", pch = 16, cex = 1)


      y.lim = c(min(Smpl$SmplCount), 1.2*max(Smpl$SmplCount))
      plot(0,0, type = 'n', xlim = rev(range(CA.age)), ylim = y.lim,
           ylab = expression(Pieces*~sample^{"-1"}), xaxt = "n", cex = 1)
      for (i in 1:nrow(overthresh.interval)) {
        polygon(c(overthresh.interval[i,], rev(overthresh.interval[i,])),
                rep(y.lim, each = 2), col = "mistyrose", border = NA)
      }
      for (i in 1:(n.peak + 1)) {
        polygon(c(underthresh.interval[i,], rev(underthresh.interval[i,])),
                rep(y.lim, each = 2), col = "lightcyan", border = NA)
      }
      lines(Smpl$Age_calBP, Smpl$SmplCount, type = 's',
            col = grDevices::grey(0.5))
      abline(h = 0, col = grDevices::grey(0.5))
      points(CA.age[peak.ind.screened] + CA.res/2,
             rep(0.9*y.lim[2], length(peak.ind.screened)),
             col = "red", pch = 3, lwd = 2)
      points(CA.age[peak.ind.notpass] + CA.res/2,
             rep(0.9*y.lim[2], length(peak.ind.notpass)),
             col = "gray", pch = 16, cex = 1)

      y.lim <- c(0, 1.2*max(Smpl$frag, na.rm = TRUE))
      plot(0, 0, type = 'n', xlim = rev(range(CA.age)), ylim = y.lim,
           ylab = expression(C["#"]/C[A] - ratio), xaxt = "n", cex = 1)
      for (i in 1:nrow(overthresh.interval)) {
        polygon(c(overthresh.interval[i,], rev(overthresh.interval[i,])),
                rep(y.lim, each = 2), col = "mistyrose", border = NA)
      }
      for (i in 1:(n.peak + 1)) {
        polygon(c(underthresh.interval[i,], rev(underthresh.interval[i,])),
                rep(y.lim, each = 2), col = "lightcyan", border = NA)
      }
      points(Smpl$Age_calBP, Smpl$frag,
             xlim = rev(range(CA.age)), pch = 16, cex = 1)
      abline(h = 0, col = grDevices::grey(0.5))


      y.lim <- c(min(CC.cpeak) + 0.01, 1.2*max(CC.cpeak))
      plot(0, 0, type = 'n', xlim = rev(range(CA.age)), ylim = y.lim,
           ylab = expression(CHAR["#"]*~"residuals"), cex = 1)

      lines(CC.age, CC.cpeak, type = 's', col = grDevices::grey(0.5))
      abline(h = 0, col = grDevices::grey(0.5))
      lines(CC.age, CC.thresh, col = 2)
      ind <- which(CC.dat[ ,12] == 1)
      points(CC.age[ind], rep(0.7*y.lim[2], length(ind)),
             col = "gray", pch = 16, cex = 1)
      ind <- which(CC.dat[ ,13] == 1)
      points(CC.age[ind], rep(0.7*y.lim[2], length(ind)),
             col = "red", pch = 3, lwd = 2)

      mtext(("Age (cal yrs BP)"), side = 1, line = 2.5, cex = 1)
    }


    # Plots comparison between fire events #####
    par(mfrow = c(1, 1))
    y.lim <- c(min(CA.cpeak), 1.7*max(CA.cpeak))
    plot(0, 0, type = 'n', xlim = rev(range(CA.age)), ylim = y.lim,
         ylab = expression(CHAR[A]*~mm^{2}), xlab = "Age cal BP")
    lines(CA.age, CA.cpeak, type = 's', col = grDevices::grey(0.5))
    abline(h = 0, col = grDevices::grey(0.5))
    lines(CA.age, CA.thresh, col = 2)
    ind <- which(CA.dat[ ,13] == 1)
    points(CA.age[ind] + CA.res/2, rep(0.8*y.lim[2], length(ind)), pch = 16,
           col = "grey", cex = 0.7)
    ind <- peak.ind.screened
    points(CA.age[ind] + CA.res/2, rep(0.8*y.lim[2], length(ind)),
           pch = 17, cex = 0.9)
    ind <- which(CC.dat[,12] == 1)
    points(CA.age[ind] + CA.res/2, rep(0.7*y.lim[2], length(ind)), pch = 16,
           col = "grey", cex = 0.8)
    ind <- which(CC.dat[,13] == 1)
    points(CA.age[ind] + CA.res/2, rep(0.7*y.lim[2], length(ind)),
           pch = 18, cex = 1.0)
    graphics::legend("topright", bty = 'n',
                     legend = c(expression(CHAR[A]*~screened),
                                expression(CHAR[C]*~screened), "Unscreened"),
                     pch = c(17,18,16), col = c("black","black","grey"),
                     ncol = 2, cex = .8)


    ## FRI plots - with boxplots -------
    if (all(is.na(FRI.CAs$FRI) == TRUE)) {
      ymax_FRI.CAs <- max(FRI.CA$FRI, na.rm = TRUE)
    } else {
      ymax_FRI.CAs <- max(FRI.CAs$FRI, na.rm = TRUE)
    }
    par(mfrow = c(2,2), mar = c(0.5, 4, 0.5, 0), oma = c(5,1,1,1), cex.lab = 1)
    y.lim <- c(min(FRI.CA$FRI, na.rm = T), 1.3*ymax_FRI.CAs)
    plot(0, 0, type = "n", ylim = y.lim, xlim = rev(range(CA.age)),
         ylab = expression(CHAR[A]*~FRI), xaxt = "n")
    points(FRI.CA[,2], FRI.CA[,3], pch = 16, col = "grey", cex = 0.7)
    points(FRI.CAs[,2], FRI.CAs[,3], pch = 17, cex = 0.9)
    graphics::legend("topleft", bty = 'n', legend = c("Screened", "Unscreened"),
                     pch = c(17,16), col = c("black","grey"), ncol = 2, cex = .8)
    graphics::boxplot(FRI.CA[,3], FRI.CAs[,3], axes = TRUE, varwidth = TRUE,
                      ylab = expression(CHAR[A]*~FRI))
    y.lim <- c(min(FRI.CC$FRI, na.rm = T), 1.3*max(FRI.CCs$FRI, na.rm = T))
    plot(0, 0, type = "n", ylim = y.lim, xlim = rev(range(CA.age)),
         ylab = expression(CHAR["#"]*~FRI))
    points(FRI.CC[,2], FRI.CC[,3], pch = 16, col = "grey", cex = 0.7)
    points(FRI.CCs[,2], FRI.CCs[,3], pch = 18, cex = 1.0)
    mtext("Age (cal yr BP)", side = 1, outer = TRUE, line = 2.5, cex = 1)
    graphics::legend("topleft", bty = 'n', legend = c("Screened", "Unscreened"),
                     pch = c(18,16), col = c("black","grey"), ncol = 2, cex = .8)
    graphics::boxplot(FRI.CC[ ,3], FRI.CCs[ ,3], axes = TRUE, varwidth = TRUE,
                      names = c("unscreened", "screened"),
                      ylab = expression(CHAR["#"]*~FRI))

    # Boxplot of charcoal-particle areas of samples that passed the screening test
    par(mfrow = c(1,1), mar = c(2, 5, 1, 1), cex = 1.2)
    graphics::boxplot(SdlArea ~ Age_calBP, data = valSdl.perc,
                      varwidth = TRUE, notch = FALSE, las = 2,
                      xlab = "",
                      ylab = expression(paste("Particle area (mm"^"2",")")))
    mtext("Age (cal yrs BP)", side = 1, line = 4, cex = 1.2)
  }

  ## Create new FireA.file with screened peaks
  fire_a_out <- FireA.file
  fire_a_out$thresh$RI_pos <- FRI.CAs$FRI
  fire_a_out$thresh$peaks.pos.age <- FRI.CAs$age_top_i
  fire_a_out$thresh$peaks.pos <- CA.dat.out$peaks_pos_sig
  fire_a_out$thresh$peaks.pos.insig <- CA.dat.out$peaks_pos_insig


  ## Gather output ------------------
  if (storedat == FALSE) {
    d_out <- fire_a_out
  } else {
    d_out <- structure(list(fire_a_out = fire_a_out, comp_perc = comp.perc))
  }
  return(d_out)
}
