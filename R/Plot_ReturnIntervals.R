#' Plot return intervals determined based on the Trend And Peaks AnalysiS (tapas).
#'
#' Requires output from the \code{\link{global_thresh}()} function.
#'
#' @param series The output of \code{\link{global_thresh}()}.
#' @param x.lim Limit age for the x-axis scale (time scale)
#' @param plot.x Boolean. If \code{FALSE} (by default),
#'               The x-axis labels are omitted.
#' @param plot.neg Boolean. If \code{TRUE} (by default),
#'                 the return intervals for both positive and negative anomalies
#'                 are plotted.
#' @export
Plot_ReturnIntervals <- function(series = NULL, x.lim = NULL, plot.x = F, plot.neg = F) {
  
  
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
  
  if (plot.neg == T) {
    y.lim <- c(0, max(c(RI_neg, RI_pos), na.rm = T))
  } else {
    y.lim <- c(0, max(c(RI_pos), na.rm = T))
  }
  
  
  # Plot data
  #par(mfrow = c(1,1), mar = c(5,4,2,0.5))
  plot(peaks_pos_ages, RI_pos, type = "l", axes = F,
       ylab = paste(proxy, "events\nreturn interval (years)"), lwd = 2,
       xlab = "", xlim = x.lim, ylim = y.lim)
  if (plot.neg == T) {
    lines(peaks_neg_ages, RI_neg, col = "blue", lwd = 2)
  }
  axis(1, labels = plot.x)
  axis(2)
  
  # mtext(paste0("thresh.value = ", d$thresh$thresh.value), side=2, las=0, line=1, cex=0.8)
}
