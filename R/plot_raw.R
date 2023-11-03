#' Plots the raw input data as a polygon and optionally also bar plots against age
#'
#' The function produces a polygon plot and optionally overlays a bar plot of
#' the raw input data.
#'
#' @param series The input data frame. A matrix with the following columns:
#'    \code{CmTop, CmBot, AgeTop, AgeBot, Volume}, and one or more columns with
#'     the data which should be re-sampled (variables).
#' @param proxy a character string with the name of the selected variable. By
#'    default (i.e. \code{proxy = NULL}), the 6th column of the input data
#'    frame is selected.
#' @param my_col a character string defining the color for the polygon.
#' @param bars logical. By default \code{bars = FALSE}. If \code{bars = TRUE},
#'    a bar plot is added as an overlay, with sample ages estimated as
#'    mid-sample ages.
#' @param y_lim a numeric vector of length 2, giving the y coordinates range.
#' @param x_lab Optional. A character string specifying the units of the age
#'              scale. By default \code{xlab = NULL} and the x-axis scale label
#'              is written as \code{age}.
#' @param y_lab Optional. A character string specifying the y-axis label.
#' @param main Optional. A character string specifying the main plot title.
#' @param \dots \dots{}
#'
#' @return A plot.
#'
#' @export
#'
#' @author Walter Finsinger
#'
#' @examples
#' co <- tapas::co_char_data
#' tapas::plot_raw(co)

plot_raw <- function(series = NULL, proxy = NULL, my_col = "grey",
                     bars = FALSE, y_lim = NULL, x_lab = NULL, y_lab = NULL,
                     main = NULL, ...) {


  ## Remove rows with NA values
  series <- series[which(complete.cases(series)), ]

  ## Get vectors for sample ages (top and bottom sample ages)
  age_top <- series[ ,3]
  age_bot <- series[ ,4]

  ## Get data from variable
  if (is.null(proxy)) {
    v <- series[ ,6] # in case a variable was not selected
  } else {
    v <- series[[proxy]] # in case a variable was selected as an argument
  }

  ## Get values for polygon plot, barplot, and x_lim
  age <- c(matrix(c(age_top, age_bot), 2, byrow = TRUE))
  countInit <- rep(v, each = 2)
  midage <- (age_top + age_bot) / 2
  x_lim <- c(max(age_bot), min(age_top))
  if (is.null(y_lim)) {
    y_lim <- c(0, max(v, na.rm = TRUE))
  }
  if (is.null(y_lab)) {
    if(is.null(proxy)) {
      y_lab <- colnames(series[6])
    } else {
      y_lab <- proxy
    }
  } else {
    y_lab <- y_lab
  }


  ## Plot polygon and barplots
  plot(age, rep.int(v, times = 2), type = "h", col = "white",
       xlim = x_lim,
       ylim = y_lim,
       xlab = x_lab,
       ylab = y_lab,
       main = main)
  polygon(x = c(age, rev(age)),
          y = c(countInit, rep_len(0, length.out = length(countInit))),
          col = my_col, border = TRUE)
  if (bars == TRUE) {
    par(new = T)
    plot(midage, v, type = "h",
         xlim = x_lim, ylim = y_lim,
         xlab = "", ylab = "")
  }
}
