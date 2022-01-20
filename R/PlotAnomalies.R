#' Plots a summary figure of the PaleoDataAnomalies analysis.
#'
#' Requires output from the \code{\link{global_thresh}()} function.
#' 
#' @param series The output of \code{\link{global_thresh}()}.
#' @param x.lim Limit age for the x-axis scale (time scale).
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
#' @export
Plot.Anomalies <- function(series = NULL, x.lim = NULL, plot.crosses = T,
                           plot.x = F, plot.neg = T) {
  

  # Gather the data
  d <- series
  proxy <- d$thresh$proxy
  i <- which(colnames(d$detr$detr) == proxy)
  
  if (is.null(x.lim)) {
    x.lim <- rev(c(min(d$int$series.int$age), max(d$int$series.int$age)))
  }
  
  y.lim <- c(min(d$int$series.int[[proxy]], na.rm = T),
             max(d$int$series.int[[proxy]], na.rm = T))
  
  # Gather data from input list
  Peaks.pos.plot <- which(d$thresh$peaks.pos == 1)
  Peaks.neg.plot <- which(d$thresh$peaks.neg == 1)
  insig.peaks <- which(d$thresh$peaks.pos.insig == 1)
  
  d.int.proxy <- d$int$series.int[[proxy]]
  d.ages <- d$int$series.int$age
  d.detr.proxy <- d$detr$detr[[proxy]]
  
  thresh.plot <- data.frame(matrix(data = NA, nrow = length(d$thresh$ages.thresh)))
  thresh.plot[ ,1] <- intersect(d$detr$detr$age, d$thresh$ages.thresh)
  colnames(thresh.plot) [1] <- "age"
  thresh.plot$int <- d.int.proxy[which(is.na(match(d.ages, thresh.plot$age)) == F)]
  thresh.plot$detr <- d.detr.proxy[which(is.na(match(d.ages, thresh.plot$age)) == F)]
  thresh.plot$thresh.pos <- thresh.plot$int - thresh.plot$detr + d$thresh$thresh.pos.sm
  thresh.plot$thresh.neg <- thresh.plot$int - thresh.plot$detr + d$thresh$thresh.neg.sm
  
  peaks.pos.poly.start <- d$thresh$ages.thresh[Peaks.pos.plot] - d$int$yr.interp/2
  peaks.pos.poly.end   <- d$thresh$ages.thresh[Peaks.pos.plot] + d$int$yr.interp/2
  peaks.neg.poly.start <- d$thresh$ages.thresh[Peaks.neg.plot] - d$int$yr.interp/2
  peaks.neg.poly.end   <- d$thresh$ages.thresh[Peaks.neg.plot] + d$int$yr.interp/2
  
  # Plot data
  if (y.lim[2] < 0) {
    #par(mfrow=c(1,1), mar=c(5,4,2,0.5))
    plot(d$int$series.int$age, d$int$series.int[ ,i], type = "l", axes = F,
         ylab = paste0(proxy, "\n", 100*d$thresh$thresh.value, "nth"), xlab = "",
         xlim = x.lim, ylim = y.lim)
    # polygon(x = c(rev(thresh.plot$age), thresh.plot$age),
    #         y = (c(rev(thresh.plot$thresh.neg), thresh.plot$thresh.pos)),
    #         col = "lightgrey", border = NA)
    lines(thresh.plot$age, thresh.plot$thresh.pos, col = "red")
    if (length(peaks.pos.poly.start > 0)) {
      for (j in 1:length(d$thresh$peaks.pos.ages)) {
        polygon(x = c(peaks.pos.poly.start[j], peaks.pos.poly.end[j],
                    peaks.pos.poly.end[j], peaks.pos.poly.start[j]),
                y = rep(y.lim, each = 2), col = "mistyrose", border = NA)
      }
    }
    if (plot.neg == T) {
      lines(thresh.plot$age, thresh.plot$thresh.neg, col = "red")
      if (length(peaks.neg.poly.start > 0)) {
        for (j in 1:length(d$thresh$peaks.neg.ages)) {
          polygon(x = c(peaks.neg.poly.start[j], peaks.neg.poly.end[j],
                      peaks.neg.poly.end[j], peaks.neg.poly.start[j]),
                  y = rep(y.lim, each = 2), col = "slategray1", border = NA)
        }
      }
    }
    par(new = T)
    plot(d$int$series.int$age, d$int$series.int[ ,i], type = "l", axes = F,
         ylab = paste0(proxy, "\n", 100*d$thresh$thresh.value, "nth"), xlab = "",
         xlim = x.lim, ylim = y.lim)
    axis(1, labels = plot.x)
    axis(2)
    lines(x = d$detr$detr$age, y = d$int$series.int[[proxy]] - d$detr$detr[[proxy]],
          col = "black", lwd = 2, xlim = x.lim, ylab = "")
    # par(new = T)
    # plot(x = thresh.plot$age, y = thresh.plot$thresh.pos, type = "l", col = "red",
    #      xlim = x.lim, ylim = y.lim, axes = F, ylab = "")
    # par(new = T)
    # plot(x = thresh.plot$age, y = thresh.plot$thresh.neg, type = "l", col = "red",
    #      xlim = x.lim, ylim = y.lim, axes = F, ylab = "")
    if (plot.crosses == T) {
      points(x = d$int$series.int$age[Peaks.pos.plot],
             y = rep(x = 1.1*y.lim[2],length(Peaks.pos.plot)),
             pch = 3, col = "red", lwd = 1.5)
      points(x = d$int$series.int$age[Peaks.neg.plot],
             y = rep(1.15*y.lim[2], length(Peaks.neg.plot)),
             pch = 3, col = "blue", lwd = 1.5)
    }
    #mtext(paste0("thresh.value  =  ", d$thresh$thresh.value), side = 2, las = 0, line = 0.5, cex = 0.8)
  } else {
    #par(mfrow = c(1,1), mar = c(5,4,2,0.5))
    plot(d$int$series.int$age, d$int$series.int[ ,i], type = "h", axes = F,
         ylab = paste(proxy, "\n", 100*d$thresh$thresh.value, "nth"),
         xlab = "", xlim = x.lim, ylim = y.lim)
    # polygon(x = c(thresh.plot$age, rev(thresh.plot$age)),
    #         y = c(thresh.plot$thresh.pos, rev(thresh.plot$thresh.neg)),
    #         col = "grey90", border = NA)
    lines(thresh.plot$age, thresh.plot$thresh.pos, col = "red")
    if (length(peaks.pos.poly.start) > 0) {
      for (j in 1:length(peaks.pos.poly.start)) {
        polygon(x = c(peaks.pos.poly.start[j], peaks.pos.poly.end[j],
                    peaks.pos.poly.end[j], peaks.pos.poly.start[j]),
                y = rep(y.lim, each = 2), col = "mistyrose", border = NA)
      }
    }
    if (plot.neg == T) {
      lines(thresh.plot$age, thresh.plot$thresh.neg, col = "red")
      if (length(peaks.neg.poly.start > 0)) {
        for (j in 1:length(peaks.neg.poly.start)) {
          polygon(x = c(peaks.neg.poly.start[j], peaks.neg.poly.end[j],
                      peaks.neg.poly.end[j], peaks.neg.poly.start[j]),
                  y = rep(y.lim, each = 2), col = "lightgrey", border = NA)
        }
      }
    }
    par(new = T)
    plot(d$int$series.int$age, d$int$series.int[ ,i], type = "h", axes = F,
         ylab = paste(proxy, "\n", 100*d$thresh$thresh.value, "nth"), xlab = "",
         xlim = x.lim, ylim = y.lim)
    axis(1, labels = plot.x)
    axis(2)
    lines(x = d$detr$detr$age, y = d$int$series.int[[proxy]] - d$detr$detr[[proxy]],
          col = "grey70", lwd = 2, xlim = x.lim, ylab = "")
    # par(new = T)
    # plot(x = thresh.plot$age, y = thresh.plot$thresh.pos, type = "l", col = "red",
    #      xlim = x.lim, ylim = y.lim, axes = F, ylab = "")
    # par(new = T)
    # plot(x = thresh.plot$age, y = thresh.plot$thresh.neg, type = "l", col = "red",
    #       xlim = x.lim, ylim = y.lim, axes = F, ylab = "")
    if (plot.crosses == T) {
      points(x = d$int$series.int$age[Peaks.pos.plot],
             y = rep(x = 0.85 * y.lim[2], length(Peaks.pos.plot)),
             pch = 3, col = "red", lwd = 1.5)
      points(x = d$int$series.int$age[insig.peaks],
             y = rep(0.85 * y.lim[2], length(insig.peaks)),
             pch = 16, col = "darkgrey", lwd = 1.5)
      if (plot.neg == T) {
      points(x = d$int$series.int$age[Peaks.neg.plot],
             y = rep(0.80*y.lim[2], length(Peaks.neg.plot)),
             pch = 3, col = "blue", lwd = 1.5)
      }
    }
    # mtext(paste0("thresh.value = ", d$thresh$thresh.value), side=2, las=0, line=1, cex=0.8)
  }
}
