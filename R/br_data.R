#' Sedimentary charcoal data (Lake Brazi)
#'
#' Lake-sediment charcoal data from Lake Brazi, as in:
#' Finsinger W, Kelly R, Fevre J, Magyari EK (2014) A guide to screening
#' charcoal peaks in macrocharcoal-area records for fire-episodes
#' reconstructions. The Holocene 24: 1002–1008: doi:10.1177/0959683614534737.
#'
#' @format A data frame with 242 rows and 7 variables:
#' \describe{
#'   \item{top_depth}{Sample top in centimeters}
#'   \item{bot_depth}{Sample bottom in centimeters}
#'   \item{top_age}{Age of sample top (cal yr BP)}
#'   \item{bot_age}{Age of sample bottom (cal yr BP)}
#'   \item{vol}{Sample volume (in cm³)}
#'   \item{char_c}{number of charcoal particles (# pieces)}
#'   \item{char_a}{sum of charcoal-particle areas (mm2)}
#' }
#' @source
#' modified from Finsinger et al. (2009) \url{https://github.com/wfinsinger/ARCO/blob/master/data_in/arco_Smpl.csv}
"br_data"
