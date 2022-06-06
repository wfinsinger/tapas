#' Performs peak-detection analysis of a time series.
#'
#' This wrapper function calls the following functions in correct order:
#' \code{\link{check_pretreat}()}, \code{\link{pretreatment_data}()},
#' \code{\link{SeriesDetrend}()}, \code{\link{global_thresh}()},
#' \code{\link{local_thresh}()}, \code{\link{Plot.Anomalies}()},
#' and \code{Plot_ReturnIntervals()}.
#'
#' @param series A \code{data.frame} with the following columns:
#'               \code{CmTop}, \code{CmBot}, \code{AgeTop},
#'               \code{AgeBot}, \code{Volume},
#'               and one or more columns with the data to be resampled.
#' @param out Desired return value: \describe{
#'            \item{"accI"}{the function returns resampled accumulation rates}
#'            \item{"conI"}{the function returns resampled concentrations}
#'            \item{"countI"}{the function returns resampled counts}
#'            }
#' @param series_name A character string defining typically the site name
#'                    (\code{NULL} by default).
#' @param first,last Age boundaries of the resampled time serie.
#'                   If unspecified (\code{first=NULL} and \code{last=NULL}),
#'                   then resampling is done over the entire sequence,
#'                   from \code{min(series$AgeTop)} to max(series$AgeBot).
#' @param yrInterp Resolution of the resampled timeseries.
#' @param smoothing_yr Smoothing-window width (in years) used for:\itemize{
#'                   \item{detrending the data}
#'                   \item{smoothing the signal-to-noise index (SNI)}
#'                   \item{the window width (in years)
#'                         for the local thresholds}
#'                   }
#' @param detr_type Smoothing function used: \describe{
#'              \item{"rob.loess"}{robust Loess}
#'              \item{"rob.lowess"}{robust Lowess with 4 iterations}
#'              \item{"mov.mode"}{moving median (aka. method #4
#'                                in Matlab CharAnalysis version)}
#'            }
#' @param proxy Select variable for the peak-detection analysis
#'              with \code{proxy = "VariableName"}.
#'              If the dataset includes only one variable,
#'              proxy does not need be specified.
#' @param t_lim Allows defining a portion of the time series.
#'              With \code{t_lim = NULL} (by default)
#'              the analysis is performed with the entire time series.
#' @param thresh_type Determines whether to use a single 'global' GMM,
#'                    or local GMMs to determine the threshold.
#'                    Defaults to \code{thresh_type = "local"}.
#' @param thresh_value Determines threshold as the nth-percentile
#'                     of the Gaussian Model of the noise component.
#'                     Defaults to \code{thresh_value = 0.95}.
#' @param noise_gmm Determines which of the two GMM components
#'                  should be considered as the noise component.
#'                  This is only needed if \code{thresh_type = "global"}.
#'                  Defaults to \code{noise_gmm = 1}.
#' @param keep_consecutive Logical. If \code{FALSE} (by default),
#'                         consecutive peak samples exceeding the threshold
#'                         are removed
#'                         and only the first (older) sample is retained.
#' @param min_CountP Probability that two resampled counts could arise
#'                   from the same Poisson distribution.
#'                   Defaults to \code{0.95}.
#'                   This is used to screen peak samples
#'                   and remove any that fail to pass the minimum-count test.
#'                   If \code{min_CountP = NULL},
#'                   the test will not be performed.
#' @param MinCountP_window Width (in years) of the search window
#'                         used for the minimum-count test.
#'                         Defaults to \code{MinCountP_window = 150}.
#' @param out_dir Path to the output directory
#'                where \code{*.pdf} figures files are written files.
#'                Defaults to \code{out_dir = NULL},
#'                so the default device is used.
#' @param plotit Logical. If \code{TRUE} (by default),
#'               results obtained with the selected settings
#'               are plotted in the display.
#' @param x_lim Age limits of the x-axis scale (the time scale).
#'              By default, \code{x_lim = t_lim}.
#' @param plot_x Logical. If \code{FALSE} (by default),
#'               the x-axis labels are omitted.
#' @param plot_neg Logical. If \code{FALSE} (by default), the return intervals
#'                 for both positive and negative anomalies are plotted.
#' @param plot_crosses Logical. If \code{TRUE} (by default),
#'                     crosses are added to indicate the location of the events.
#' @param plot_neg Logical. If \code{FALSE} (by default),
#'                 both positive and negative events
#'                 are marked with colored shaded areas.
#' @param sens Logical. Determines whether a sensitivity analysis is performed.
#'             This involves performing the peak-detection analysis
#'             with different smoothing-window widths (\code{smoothing_yr}),
#'             and plotting boxplots for the SNI distributions,
#'             the Return Intervals, and the number of detected events.
#' @param smoothing_yr_seq The smoothing-window widths used
#'                         for the sensitivity analysis.
#'                         If \code{smoothing_yr_seq = NULL} (by default),
#'                         the analysis is performed
#'                         for the following smoothing-window widths:
#'                         500, 600, 700, 800, 900, 1000 years.
#'
#' @importFrom graphics boxplot layout
#'
#' @export
peak_detection <- function(series = NULL, out = "accI", proxy = NULL,
                           first = NULL, last = NULL, yrInterp = NULL,
                           smoothing_yr = 500, detr_type = "mov.median",
                           thresh_type = "local", thresh_value = 0.95,
                           t_lim = NULL, noise_gmm = 1, keep_consecutive = F,
                           min_CountP = 0.05, MinCountP_window = 150,
                           series_name = NULL, out_dir = NULL,
                           plotit = TRUE, x_lim = t_lim,
                           plot_crosses = TRUE, plot_x = TRUE,
                           plot_neg = FALSE,
                           sens = TRUE, smoothing_yr_seq = NULL) {

  ## Initial check-ups ####

  ## Check if data is formatted correctly
  series <- tapas::check_pretreat(series)

  ## Check if arguments for detrending are coherent
  detr_options <- c("rob.loess", "rob.lowess", "mov.median")
  if (!(detr_type %in% detr_options)) {
    stop(
      "Unrecognized detrending type: '", detr_type, "'. ",
      "Accepted values: ", paste(detr_options, collapse = ", "), "."
    )
  }

  ## Check if the variable name (proxy) is coherent
  proxies_options <- colnames(series[ 6:dim(series)[2] ])

  if (is.null(proxy) & length(proxies_options) > 1) {
    # stop(
    #   "Unrecognized variable name: '", proxy, "'. ",
    #   "Accepted names: ", paste(proxies_options, collapse = ", "), "."
    # )
    stop(paste0("Fatal error: Unrecognized variable name (proxy = ). ",
                "Accepted names: ", paste(proxies_options, collapse = ", "), "."
    ))
  }

  if (is.null(proxy) & length(proxies_options) == 1) {
    proxy <- colnames(series) [6]
  }

  if (!(proxy %in% proxies_options)) {
    stop(
      "Unrecognized variable name: '", proxy, "'. ",
      "Accepted names: ", paste(proxies_options, collapse = ", "), "."
    )
  }


  ## Resample data ####
  d_i <- tapas::pretreatment_data(series = series, out = out,
                                  series.name = series_name,
                                  first = first, last = last, yrInterp = yrInterp)


  ## Detrend resampled data ####
  d_detr <- tapas::SeriesDetrend(series = d_i, smoothing.yr = smoothing_yr,
                                 detr.type = detr_type, out.dir = out_dir)


  ## Screen peaks using GMM-inferred threshold(s) and minimum-count test and ####
  ## evaluate suitability of the record with signal-to-noise index
  if (out == "accI") {
    proxy <- paste0(proxy,"AR")
  }

  if (thresh_type == "global") {
    d_thresh <- tapas::global_thresh(series = d_detr, proxy = proxy,
                                     t.lim = t_lim,
                                     thresh.value = thresh_value,
                                     noise.gmm = noise_gmm,
                                     smoothing.yr = smoothing_yr,
                                     keep_consecutive = keep_consecutive,
                                     minCountP = min_CountP,
                                     MinCountP_window = MinCountP_window,
                                     out.dir = out_dir,
                                     plot.global_thresh = plotit)
    d_thresh$thresh$thresh_type <- "global"
  } else {
    d_thresh <- tapas::local_thresh(series = d_detr, proxy = proxy,
                                    t.lim = t_lim,
                                    thresh.value = thresh_value,
                                    thresh.yr = smoothing_yr,
                                    smoothing.yr = smoothing_yr,
                                    keep_consecutive = keep_consecutive,
                                    minCountP = min_CountP,
                                    MinCountP_window = MinCountP_window,
                                    out.dir = out_dir,
                                    plot.local_thresh = plotit)
    d_thresh$thresh$thresh_type <- "local"
  }



  ## Perform sensitivity analysis ####
  if (sens == T) {

    ## Effect of different smoothing-window widths

    # Select smoothing windows
    if (is.null(smoothing_yr_seq)) {
      smoothing_yr_seq <- seq(from = 500, to = 1000, by = 100)
    }

    ## Prepare space to gather data generated in each loop
    sens_sni_pos <- as.data.frame(matrix(NA, ncol = length(smoothing_yr_seq),
                                         nrow = dim(d_thresh$int$series.int)[1]))
    names(sens_sni_pos) <- c(smoothing_yr_seq)


    sens_peaks_pos <- as.data.frame(matrix(NA, ncol = 2,
                                           nrow = length(smoothing_yr_seq)))
    sens_peaks_pos[ ,1] <- c(smoothing_yr_seq)

    sens_ri_pos <- vector('list', length(smoothing_yr_seq))
    names(sens_ri_pos) <- c(smoothing_yr_seq)


    ## Perform peak_detection() for each of the smoothing_yr_seq values:
    for (i in seq_along(smoothing_yr_seq)) {
      smoothing_yr <- smoothing_yr_seq[i]
      d_detr_i <- SeriesDetrend(series = d_i, smoothing.yr = smoothing_yr,
                                detr.type = detr_type, plot_pdf = F)


      ## Screen peaks using GMM-inferred threshold(s) and minimum-count test and
      ## evaluate suitability of the record with signal-to-noise index

      # sens_plots <- F # such that pdf files are not produced for the sensitivity runs

      if (thresh_type == "global") {
        d_thresh_i <- tapas::global_thresh(series = d_detr_i, proxy = proxy,
                                           t.lim = t_lim,
                                           thresh.value = thresh_value,
                                           noise.gmm = noise_gmm,
                                           smoothing.yr = smoothing_yr,
                                           keep_consecutive = keep_consecutive,
                                           minCountP = min_CountP,
                                           MinCountP_window = MinCountP_window,
                                           plot.global_thresh = F)
        d_thresh_i$thresh$thresh_type <- "global"
      } else {
        d_thresh_i <- tapas::local_thresh(series = d_detr_i, proxy = proxy,
                                          t.lim = t_lim,
                                          thresh.value = thresh_value,
                                          thresh.yr = smoothing_yr,
                                          smoothing.yr = smoothing_yr,
                                          keep_consecutive = keep_consecutive,
                                          minCountP = min_CountP,
                                          MinCountP_window = MinCountP_window,
                                          plot.local_thresh = F)
        d_thresh_i$thresh$thresh_type <- "local"
      }

      # Gather data
      sens_sni_pos[ ,i] <- d_thresh_i$thresh$SNI_pos$SNI_raw
      sens_peaks_pos[i,2] <- sum(d_thresh_i$thresh$peaks.pos)
      sens_ri_pos[[i]] <- d_thresh_i$thresh$RI_pos
    }

    ## Replace Inf values with NA values in SNI data frame
    sens_sni_pos <- do.call(data.frame,
                            lapply(sens_sni_pos,
                                   function(x) replace(x,
                                                       is.infinite(x),
                                                       NA)))

    ## Plot output of sensitivity analysis ####
    layout(1)
    par(mfrow = c(3,1), mar = c(2,4,2,2), oma = c(2,1,0,1), cex = 1)

    boxplot(sens_sni_pos, xlab = "",
            ylab = "SNI",
            notch = F, axes = F,
            ylim = c(0, max(sens_sni_pos[ ,1:length(smoothing_yr_seq)],
                            na.rm = T)),
            main = paste("Sensitivity for", thresh_type,
                         "threshold at", thresh_value))
    abline(h = 3, lty = 2)
    axis(1, at = 1:length(smoothing_yr_seq), labels = c(smoothing_yr_seq))
    axis(2)
    mtext("smoothing-window width (years)", side = 1, line = 2)

    boxplot(sens_ri_pos, xlab = "",
            ylab = "Return intervals", axes = F,
            ylim = c(0, 1.3*max(unlist(sens_ri_pos), na.rm = T)),
            na.rm = T)
    axis(1, at = 1:length(smoothing_yr_seq), labels = c(smoothing_yr_seq))
    axis(2)
    mtext("smoothing-window width (years)", side = 1, line = 2)


    plot(sens_peaks_pos[ ,1], sens_peaks_pos[ ,2],
         ylab = "# of Positive peaks",
         type = "o", col = "blue",
         lwd = 2, axes = F)
    axis(1, labels = T)
    axis(2)
    mtext("smoothing-window width (years)", side = 1, line = 2)
  }



  ## Plot results obtained with the selected parameters ####
  if (plotit == T) {
    layout(1)
    par(mfrow = c(2,1), mar = c(0,5,2,0.5), oma = c(4,1,0,0))
    Plot.Anomalies(series = d_thresh, x.lim = x_lim,
                   plot.crosses = plot_crosses, plot.x = F,
                   plot.neg = plot_neg)
    Plot_ReturnIntervals(series = d_thresh, x.lim = x_lim,
                         plot.x = plot_x, plot.neg = plot_neg)
    mtext("Age (cal yr BP)", side = 1, line = 2.5)
  }


  ## Return output to environment ####
  return(d_thresh)
}
