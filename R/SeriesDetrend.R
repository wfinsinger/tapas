#
#
#  ***************************************************************************

#' Detrend time series.
#'
#' Processes the output of pretreatment_data() function
#' to de-trend the time series.
#'
#' @param series Output of pretreatment_data().
#' @param smoothing.yr Smoothing window width.
#' @param detr.type Smoothing function used: \describe{
#'              \item{"rob.loess"}{robust Loess}
#'              \item{"rob.lowess"}{robust Lowess with 4 iterations}
#'              \item{"mov.mode"}{moving median (aka. method #4
#'                                in Matlab CharAnalysis version)}
#'            }
#' @param out.dir Folder to output the generated figures to.
#'                Defaults to NULL, and the figures are generated
#'                to default device instead.
#' @param plot_pdf Logical. If \code{TRUE}, a \code{*.pdf} file
#'                 will be written in the \code{out.dir} folder.
#'
#' @return A list similar to \code{series} with additional detrended data
#'         appended.
#'
#' @importFrom dplyr %>% select everything
#' @importFrom stats loess lowess predict
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis lines mtext par
#'
#' @export
SeriesDetrend <- function(series = NULL, smoothing.yr = NULL, detr.type = "rob.loess",
                          out.dir = NULL, plot_pdf = T) {
  
  # Determine path to Figure-output folder and create folder if it does not exist
  if (plot_pdf == T && !is.null(out.dir)) {
    out.path <- paste0("./", out.dir, "/")
    
    if (dir.exists(out.path) == F) {
      dir.create(out.path)
    }
  }
  
  # Extract data and parameters from input list
  a <- series$int$series.int
  s.name <- series$raw$series.name
  yr.interp <- series$int$yr.interp
  
  # Set Proportion of datapoints used for Loess and Lowess smoothing functions:
  n.smooth <- round(smoothing.yr/yr.interp)
  span <- n.smooth/dim(a)[1]
  
  # Get trends for every variable in dataset
  a.trend <- data.frame(a$age)
  colnames(a.trend) <- "age"
  
  if (detr.type == "rob.lowess") {
    for (i in 2:ncol(a)) {
      a.trend[ ,i] <- lowess(x = a[ ,i], f = span, iter = 4)$y
    }
  }
  
  if (detr.type == "rob.loess") {
    for (i in 2:ncol(a)) {
      loess.i <- loess(formula = a[ ,i] ~ age, data = a, span = span,
                       family = "symmetric")
      a.trend[ ,i] <- predict(object = loess.i, newdata = a.trend$age)
    }
  }
  
  # Moving median (Method #4 in Matlab version; translated from Matlab version)
  if (detr.type == "mov.median") {
    for (j in 2:ncol(a)) {
      for (i in 1:dim(a)[1]) {
        if (i <= round(n.smooth / 2)) { # if 1/2 n.smooth twards start
          #ai_t <- a[1:round(n.smooth / 2), j] # value for year t
          ai_t <- a[1:i + 1, j]
          a.trend[i,j] <- median(ai_t)
        }
        if (i >= dim(a)[1] - round(n.smooth / 2)) { # if 1/2 s twords end
          #ai_t <- a[dim(a)[1] - round(n.smooth / 2):dim(a)[1], j]
          ai_t <- a[i:dim(a)[1], j]
          a.trend[i,j] <- median(ai_t)
        }
        if (i > round(n.smooth/2) && i < dim(a)[1] - round(n.smooth / 2)) { # else, you're in the middle of the record
          ai_t <- a[round(i - 0.5 * n.smooth):round(i + 0.5 * n.smooth), j]
          a.trend[i,j] <- median(ai_t)
        }
      }
      a.trend[ ,j] <- lowess(x = a.trend[ ,j], f = span, iter = 0)$y
    }
  }
  
  colnames(a.trend) <- colnames(a)
  
  
  # Then detrend every variable in dataset
  a.detr <- as.data.frame(a[ ,2:ncol(a)] - a.trend[ ,2:ncol(a.trend)])
  colnames(a.detr) <- colnames(a) [2:ncol(a)]
  a.detr$age <- a.trend$age
  a.detr <- a.detr %>% select("age", everything()) # reorders the columns
  # a.detr <- a.detr[which(complete.cases(a.detr)), ] # removed this because in some datasets some variables are incomplete!
  
  
  # And plot the detrended dataseries
  x.lim <- rev(range(a.detr$age))
  a.names <- colnames(a)
  
  if (plot_pdf == T) {
    if (!is.null(out.dir)) {
      pdf(paste0(out.path, s.name, "_DetrendedSeries.pdf"))
    }
    for (i in 2:dim(a.detr)[2]) {
      i.name <- colnames(a)[i]
      par(mfrow = c(2,1), oma = c(2,0.5,0.5,0.5), mar = c(2,5,1,0.5))
      plot(a$age, a[ ,i], type = "s", axes = F, ylab = a.names[i], xlab = "", xlim = x.lim,
           main = paste(i.name, "trend with a", smoothing.yr,"yr", detr.type, "smoothing"))
      axis(side = 1, labels = F, tick = T)
      axis(2)
      lines(a.trend$age, a.trend[ ,i], col = "blue4", lwd = 2)
      
      plot(a.detr$age, a.detr[ ,i], type = "s", xlab = "Age", ylab = "Residuals\n(data-trend)",
           axes = F, xlim = x.lim)
      abline(h = 0, col = "red")
      axis(side = 1, labels = T, tick = T)
      axis(2)
      mtext(text = "Age (cal yr BP)", side = 1, line = 2.3)
    }
    if (!is.null(out.dir)) {
      dev.off()
    }
  }
  
  # Prepare output
  out1 <- structure(list(detr = a.detr, smoothing.yr = smoothing.yr,
                         detr.type = detr.type, series.name = s.name))
  out2 <- append(series, list(out1))
  names(out2) [4] <- "detr"
  
  # Get output
  return(out2)
  
}
