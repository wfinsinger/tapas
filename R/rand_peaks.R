#' A random series with peaks that are proportional to background CHAR
#'
#' A random series with charcoal peaks that are proportional to background
#' CHAR, i.e. peaks get bigger as background gets bigger, half way through the
#' record. The frequency of peaks does not change throughout this record. Thus,
#' peak-detection methods should identify no difference in peak frequency
#' between 1st and 2nd half of this record.
#'
#' This record is similar to that displayed in Figure 2 ("Scenario 2") in
#' Higuera et al (2010) Peak detection in sediment-charcoal records:
#' impacts of alternative data analysis methods on fire-history interpretations.
#' International Journal of Wildland Fire 19:996-1014.
#'
#' @format A data frame with 500 rows and 6 variables:
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
"rand_peaks"
