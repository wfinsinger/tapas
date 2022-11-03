#' Plots the raw input data as a polygon and bar plots against age
#'
#' This functions produces a polygon plot and optionally overlays a bar plot of
#' the raw input data.
#'
#' @param series The input data frame. A matrix with the following columns:
#'    \code{CmTop, CmBot, AgeTop, AgeBot, Volume},
#'    and one or more columns with the data
#'    which should be re-sampled (variables).
#' @param proxy a character string with the name of the selected variable. By
#'    default (proxy = NULL) the 6th column of the input data frame is selected.
#' @param my_col a character string defining the color for the polygon.
#' @param bars logical. By default bars = FALSE. If bars = TRUE, a bar plot is
#'    added as an overlay, with sample ages estimated as mid-sample ages.
#'
#' @return A plot.
#'
#' @export
#'
#' @author Walter Finsinger
#'

plot_raw <- function(series = NULL, proxy = NULL, my_col = "grey",
                     bars = FALSE) {

  # series = co
  # proxy = NULL
  # bars = TRUE


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
  x_lim = c(max(age_top), min(age_top))
  y_lim <- c(0, max(v, na.rm = TRUE))

  ## Plot polygon and barplots
  plot(age_top, v, type = "h", col = "white",
       xlim = x_lim,
       ylim = y_lim,
       xlab = "Age (cal BP)",
       ylab = "Number of pieces per sample")
  polygon(x = c(age, rev(age)),
          y = c(countInit,
                rep_len(0, length.out = length(countInit))),
          col = my_col, border = TRUE)
  if (bars == TRUE) {
    par(new = T)
    plot(midage, v, type = "h",
         xlim = x_lim, ylim = y_lim,
         xlab = "", ylab = "")
  }
}
