#' Plot return intervals determined based on the Trend And Peaks AnalysiS (tapas).
#'
#' Requires output from one of the following functions:
#' \code{\link{global_thresh}()}, \code{\link{local_thresh}()},
#' or \code{\link{peak_detection}()}.
#'
#' @param series The output of \code{\link{global_thresh}()},
#'               \code{\link{local_thresh}()},
#'               or \code{\link{peak_detection}()}.
#' @param x.lim Limit age for the x-axis scale (time scale)
#' @param y_lim A vector defining the limits of the y-axis scale.
#' @param plot.x Boolean. If \code{FALSE} (by default),
#'               The x-axis labels are omitted.
#' @param plot.neg Boolean. If \code{TRUE} (by default),
#'                 the return intervals for both positive and negative anomalies
#'                 are plotted.
#' @export
Plot_ReturnIntervals <- function(series = NULL, x.lim = NULL, y_lim = NULL,
                                 plot.x = F, plot.neg = F) {

  ## Gather data
  d <- series
  proxy <- d$thresh$proxy
  peaks_neg_ages <- d$thresh$peaks.neg.ages
  peaks_pos_ages <- d$thresh$peaks.pos.ages
  RI_neg <- d$thresh$RI_neg
  RI_pos <- d$thresh$RI_pos

  if (is.null(x.lim)) {
    x.lim <- rev(c(min(d$int$series.int$age), max(d$int$series.int$age)))
  }

  if (is.null(y_lim)) {
    if (plot.neg == T) {
      y_lim <- c(0, max(c(RI_neg, RI_pos), na.rm = T))
    } else {
      y_lim <- c(0, max(c(RI_pos), na.rm = T))
    }
  }

  # Plot data
  plot(peaks_pos_ages, RI_pos, type = "l", axes = F,
       ylab = paste(proxy, "events\nreturn interval (years)"), lwd = 2,
       xlab = "", xlim = x.lim, ylim = y_lim)
  if (plot.neg == T) {
    lines(peaks_neg_ages, RI_neg, col = "blue", lwd = 2)
  }
  axis(1, labels = plot.x)
  axis(2)
}
