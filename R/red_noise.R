#' A random series with no peaks, only red noise
#'
#' A random series with no distinct charcoal peaks, only red noise
#' (random variation with 1st-order autocorrelation of ca 0.3).
#' Any real record should have a median signal-to-noise index greater than that
#' from the series of red noise.
#'
#' This record is similar to that displayed in Figure 2 ("RN") in
#' Kelly et al. (2011) A signal-to-noise index to quantify the potential for
#' peak detection in sediment-charcoal records. Quaternary Research 75:11-17.
#'
#' The original record was limited to samples having \code{age_bot < 10,000}
#' years to reduce its size.
#'
#' @format A data frame with 1000 rows and 6 variables:
#' \describe{
#'   \item{cm_top}{Sample top in centimeters}
#'   \item{cm_bot}{Sample bottom in centimeters}
#'   \item{age_top}{Age of sample Sample top}
#'   \item{age_bot}{Age of sample Sample bottom}
#'   \item{vol}{Sample volume in cm³}
#'   \item{char}{Sample count}
#' }
#'
#' @source
#' \url{https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/charanalysis/sample_records_2012.zip}
"red_noise"
