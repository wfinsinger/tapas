#' Extract data from a pretreatment_data() list to determine trend with a GAM.
#'
#' @description
#' Gathers data from a list produced with the \code{tapas::pretreatment_data()}
#' function, and returns a data frame that can be used to determine the trend
#' based on a Generalized Additive Model with the \code{mgcv::gam()} function.
#'
#' @details
#' The \code{tapas::pretreatment_data()} function sometimes may generate NA
#' values in the cm_top and cm_bot columns for the first and/or the last samples
#' in a temporal series. The \code{tapas::tapasa2mgcv()} function removes these
#' NA values automatically by taking the sediment-accumulation rate of adjacent
#' binned samples.
#'
#' @param series A list that was generated with the
#' \code{tapas::pretreatment_data()} function.
#'
#' @seealso [mgcv2tapas()]
#'
#' @return
#' The returned data frame includes the following columns: cm_top, cm_bot,
#' age_top, age_bot, volI, and in addition one column for each of the variables
#' that were included in the input data frame for the
#' \code{tapas::pretreatment_data()} function.
#' The returned age corresponds to the age_top of the samples.
#'
#' @author Walter Finsinger
#'
#' @importFrom dplyr mutate
#'
#' @export

tapas2mgcv <- function(series = NULL) {

  ## Gather data from the list that was produced by tapas::pretreatment_data()
  age_top <- series$int$series.int$age
  age <- as.data.frame(age_top)
  age <- age %>%
    dplyr::mutate(age_bot = age_top + series$int$yr.interp)
  data <- series$int$series.int %>%
    dplyr::select(-age)
  volI <- series$int$volI
  volI <- as.data.frame(volI)


  ## Get cm_top and cm_bot ####
  cm_top <- series$int$cmI

  ## If the pretreatment() function returns NA values for some cm_top samples:
  # Check if any of the first 10 cm_top sample(s) == NA
  if (is.na(cm_top[1])) {
    cm_top_na_id <- which(is.na(cm_top[1:10]))
    cm_top_na_id <- rev(cm_top_na_id) ## reverse vector
    # get the cm_top increment from next samples
    cm_top_first_diff <- abs(cm_top[cm_top_na_id[1] + 2] -
                               cm_top[cm_top_na_id[1] + 1])
    # replace NA values
    for (i in 1:length(cm_top_na_id)) {
      cm_top[cm_top_na_id[i]] <- cm_top[cm_top_na_id[i] + 1] -
        cm_top_first_diff
    }
  }

  # Check if any of the last 10 cm_top sample(s) == NA
  if (is.na(cm_top[length(cm_top)])) {
    cm_top_na_id <-
      which(is.na(cm_top[(length(cm_top) - 9):length(cm_top)]))
    cm_top_na_id <- cm_top_na_id - 10 + length(cm_top)
    # get the cm_top increment from previous samples
    cm_top_last_diff <- abs(cm_top[cm_top_na_id[1] - 2] -
                              cm_top[cm_top_na_id[1] - 1])
    # replace NA values
    for (i in 1:length(cm_top_na_id)) {
      cm_top[cm_top_na_id[i]] <- sum(cm_top[cm_top_na_id[i] - 1],
                                     cm_top_last_diff)
    }
  }

  # Get cm_bot for the binned samples except for the last one
  cm_bot <- cm_top[2:length(cm_top)]
  cm_bot_last_increment <- abs(cm_bot[length(cm_bot) - 1] -
                                 cm_bot[length(cm_bot)])
  # Add cm_bot for last binned sample
  cm_bot[length(cm_bot) + 1] <- sum(cm_bot[length(cm_bot) - 1],
                                    cm_bot_last_increment)

  ## Prepare output
  out <- cbind(cm_top, cm_bot, age, volI, data)

  ## Return output
  return(out)
}
