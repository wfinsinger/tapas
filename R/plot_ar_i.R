#' Plot to compare re-sampled and non-resampled accumulation rates
#'
#' Requires an object returned by the tapas::pretreatment() function.
#' Non-resampled accumulation rates are presented using grey bars. The re-sampled
#' accumulation rates at equal sampling intervals are presented by a black curve.
#' Modified from the plot.CHAR() function (paleofire v1.2.4 R package).
#'
#' @param x An list returned by the tapas::pretreatment() function.
#' @param xlim Optional. Determines the limits of the age scale.
#'              By default (\code{xlim = NULL}), the age-scale limits are in the
#'              form: \code{xlim = c(older limit, younger limit)}.
#' @param ylim Optional. Determines the limits of the y-axis scale.
#'              By default, \code{ylim = NULL}.
#' @param xlab Optional. A character string specifying the units of the age scale.
#'              By default \code{xlab = NULL} and the x-axis scale label is
#'              written as \code{age}.
#' @param ylab Optional. A character string specifying the units of the y-axis.
#'              By default \code{ylab = NULL} and the y-axis scale label is
#'              written as \code{accInit}.
#' @param frame Logical. Determines if the plot is framed. By default,
#'              \code{frame = TRUE}.
#' @param main Optional. A character string specifying the main plot title.
#' @param \dots \dots{}
#'
#' @author Olivier Blarquez
#' @author Walter Finsinger
#'
#' @examples
#' \dontrun{
#' # Here, we use the charcoal record from Code Lake (Higuera et al., 2009)
#' # Load raw charcoal data:
#' co <- tapas::co_char_data
#' c <- co[, 6] # charcoal counts
#' p <- co[, 1:5] # CmTop, CmBot, AgeTop, AgeBot, Volume
#'
#' # Calculate resampled charcoal accumulation rate (CHAR, as pieces cm-2 yr-1)
#' co_CHAR <- tapas::pretreatment(params = p, serie = c, Int = TRUE)
#' tapas::plot_ar_i(co_CHAR)
#' }
#' @importFrom stats na.omit
#' @export
plot_ar_i <- function(x, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
                      frame = TRUE, main = NULL, ...) {

  # Gather the values used for the plot
  age <- c(matrix(c(x$ageTop, x$ageBot), 2, byrow = TRUE))
  accInit <- rep(x$acc, each = 2)
  ageInt <- c(matrix(c(x$ybpI, x$ybpI + x$yrInterp), 2, byrow = TRUE))
  accInt <- rep(x$accI, each = 2)
  if (is.null(ylim)) ylim <- c(min(stats::na.omit(accInit)),
                               max(stats::na.omit(accInit)))
  if (is.null(xlim)) xlim <- c(max(age), min(age))

  # Make the plot
  plot(age, accInit,
       type = "l", col = "grey", ylim = ylim, xlim = xlim, yaxs = "i",
       xlab = xlab, ylab = ylab, frame = frame, main = main
  )
  polygon(c(age, age[length(age)]), c(accInit, -1e6), col = "grey")
  lines(age, accInit, type = "l", col = "grey")
  lines(ageInt, accInt, type = "l")
}
