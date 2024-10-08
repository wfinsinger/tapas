#' Plots a summary figure of results obtained with 'tapas'.
#'
#' Requires output from one of the following functions:
#' \code{\link{global_thresh}()}, \code{\link{local_thresh}()},
#' or \code{\link{peak_detection}()}.
#'
#' @param series The output of \code{\link{global_thresh}()},
#'        \code{\link{local_thresh}()}, or \code{\link{peak_detection}()}.
#' @param x.lim A vector defining the limits for the x-axis scale (time scale).
#' @param y_lim A vector defining the limits for the y-axis scale.
#' @param plot.crosses Boolean. If \code{TRUE} (by default),
#'                     crosses are displayed to indicate
#'                     the location of anomalies.
#' @param plot.x Boolean. If \code{FALSE} (by default),
#'               the x-axis labels are omitted.
#' @param plot.neg Boolean. If \code{FALSE} (by default),
#'                 both positive and negative anomalies
#'                 are marked with colored shaded areas.
#'
#' @importFrom graphics polygon
#'
#' @author Walter Finsinger
#'
#' @export
Plot.Anomalies <- function(series = NULL, x.lim = NULL, y_lim = NULL,
                           plot.crosses = T, plot.x = F, plot.neg = T) {

  # Gather the data
  d <- series
  proxy <- d$thresh$proxy
  i <- which(colnames(d$detr$detr) == proxy)

  if (is.null(x.lim)) {
    x.lim <- rev(c(min(d$int$series.int$age), max(d$int$series.int$age)))
  }

  if (is.null(y_lim)) {
  y_lim <- c(min(d$int$series.int[[proxy]], na.rm = T),
             max(d$int$series.int[[proxy]], na.rm = T))
  }

  # Gather data from input list
  Peaks.pos.plot <- which(d$thresh$peaks.pos == 1)
  Peaks.neg.plot <- which(d$thresh$peaks.neg == 1)
  insig.peaks <- which(d$thresh$peaks.pos.insig == 1)

  d.int.proxy <- d$int$series.int[[proxy]]
  d.ages <- d$int$series.int$age
  d.detr.proxy <- d$detr$detr[[proxy]]

  thresh.plot <- data.frame(matrix(data = NA,
                                   nrow = length(d$thresh$ages.thresh)))
  thresh.plot[ ,1] <- intersect(d$detr$detr$age, d$thresh$ages.thresh)
  colnames(thresh.plot) [1] <- "age"
  thresh.plot$int <- d.int.proxy[which(is.na(match(d.ages,
                                                   thresh.plot$age)) == F)]
  thresh.plot$detr <- d.detr.proxy[which(is.na(match(d.ages,
                                                     thresh.plot$age)) == F)]
  thresh.plot$thresh.pos <- thresh.plot$int - thresh.plot$detr +
    d$thresh$thresh.pos.sm
  thresh.plot$thresh.neg <- thresh.plot$int - thresh.plot$detr +
    d$thresh$thresh.neg.sm

  peaks.pos.poly.start <- d$thresh$ages.thresh[Peaks.pos.plot] -
    d$int$yr.interp/2
  peaks.pos.poly.end   <- d$thresh$ages.thresh[Peaks.pos.plot] +
    d$int$yr.interp/2
  peaks.neg.poly.start <- d$thresh$ages.thresh[Peaks.neg.plot] -
    d$int$yr.interp/2
  peaks.neg.poly.end   <- d$thresh$ages.thresh[Peaks.neg.plot] +
    d$int$yr.interp/2

  # Plot data
  if (y_lim[2] < 0) {
    par(mar = c(3,5,2,0.5))
    plot(d$int$series.int$age, d$int$series.int[ ,i], type = "l", axes = F,
         ylab = paste0(proxy, "\n", 100*d$thresh$thresh.value, "nth"),
         xlab = "",
         xlim = x.lim, ylim = y_lim)
    lines(thresh.plot$age, thresh.plot$thresh.pos, col = "red")
    if (length(peaks.pos.poly.start > 0)) {
      for (j in 1:length(d$thresh$peaks.pos.ages)) {
        polygon(x = c(peaks.pos.poly.start[j], peaks.pos.poly.end[j],
                    peaks.pos.poly.end[j], peaks.pos.poly.start[j]),
                y = rep(y_lim, each = 2), col = "mistyrose", border = NA)
      }
    }
    if (plot.neg == T) {
      lines(thresh.plot$age, thresh.plot$thresh.neg, col = "red")
      if (length(peaks.neg.poly.start > 0)) {
        for (j in 1:length(d$thresh$peaks.neg.ages)) {
          polygon(x = c(peaks.neg.poly.start[j], peaks.neg.poly.end[j],
                      peaks.neg.poly.end[j], peaks.neg.poly.start[j]),
                  y = rep(y_lim, each = 2), col = "slategray1", border = NA)
        }
      }
    }
    par(new = T)
    plot(d$int$series.int$age, d$int$series.int[ ,i], type = "l", axes = F,
         ylab = paste0(proxy, "\n", 100*d$thresh$thresh.value, "nth"),
         xlab = "",
         xlim = x.lim, ylim = y_lim)
    axis(1, labels = plot.x)
    axis(2)
    lines(x = d$detr$detr$age,
          y = d$int$series.int[[proxy]] - d$detr$detr[[proxy]],
          col = "black", lwd = 2, xlim = x.lim, ylab = "")
    if (plot.crosses == T) {
      points(x = d$int$series.int$age[Peaks.pos.plot],
             y = rep(x = 1.1*y_lim[2], length(Peaks.pos.plot)),
             pch = 3, col = "red", lwd = 1.5)
      points(x = d$int$series.int$age[Peaks.neg.plot],
             y = rep(1.15*y_lim[2], length(Peaks.neg.plot)),
             pch = 3, col = "blue", lwd = 1.5)
    }
  } else {
    par(mar = c(3,5,2,0.5))
    plot(d$int$series.int$age, d$int$series.int[ ,i], type = "h", axes = F,
         ylab = paste(proxy, "\n", 100*d$thresh$thresh.value, "nth"),
         xlab = "", xlim = x.lim, ylim = y_lim)
    lines(thresh.plot$age, thresh.plot$thresh.pos, col = "red")
    if (length(peaks.pos.poly.start) > 0) {
      for (j in 1:length(peaks.pos.poly.start)) {
        polygon(x = c(peaks.pos.poly.start[j], peaks.pos.poly.end[j],
                    peaks.pos.poly.end[j], peaks.pos.poly.start[j]),
                y = rep(y_lim, each = 2), col = "mistyrose", border = NA)
      }
    }
    if (plot.neg == T) {
      lines(thresh.plot$age, thresh.plot$thresh.neg, col = "red")
      if (length(peaks.neg.poly.start > 0)) {
        for (j in 1:length(peaks.neg.poly.start)) {
          polygon(x = c(peaks.neg.poly.start[j], peaks.neg.poly.end[j],
                      peaks.neg.poly.end[j], peaks.neg.poly.start[j]),
                  y = rep(y_lim, each = 2), col = "lightgrey", border = NA)
        }
      }
    }
    par(new = T)
    plot(d$int$series.int$age, d$int$series.int[ ,i], type = "h", axes = F,
         ylab = paste(proxy, "\n", 100*d$thresh$thresh.value, "nth"),
         xlab = "",
         xlim = x.lim, ylim = y_lim)
    axis(1, labels = plot.x)
    axis(2)
    lines(x = d$detr$detr$age,
          y = d$int$series.int[[proxy]] - d$detr$detr[[proxy]],
          col = "grey70", lwd = 2, xlim = x.lim, ylab = "")
    if (plot.crosses == T) {
      points(x = d$int$series.int$age[Peaks.pos.plot],
             y = rep(x = 0.85 * y_lim[2], length(Peaks.pos.plot)),
             pch = 3, col = "red", lwd = 1.5)
      points(x = d$int$series.int$age[insig.peaks],
             y = rep(0.85 * y_lim[2], length(insig.peaks)),
             pch = 16, col = "darkgrey", lwd = 1.5)
      if (plot.neg == T) {
      points(x = d$int$series.int$age[Peaks.neg.plot],
             y = rep(0.80*y_lim[2], length(Peaks.neg.plot)),
             pch = 3, col = "blue", lwd = 1.5)
      }
    }
  }
}
