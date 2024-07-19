#' Check the input data
#'
#' This functions does some check-up on the data frame used for
#' the \code{tapas::pretreatment_data()} function,
#' which requires that the depth and age scales are continuous.
#' In other words, for every *i-th* row, it requires that \itemize{
#' \item{\code{CmBot[i] > CmTop[i]}}
#' \item{\code{AgeBot[i] > AgeTop[i]}}
#' \item{there shouldn't be duplicate values
#'       in the \code{CmTop, CmBot, AgeTop}, and \code{AgeBot} columns}.
#' }
#'
#'
#' @details
#'
#' If any of the following is true: \itemize{
#' \item{any \code{CmBot[i] < CmTop[i]}}
#' \item{\code{AgeBot[i] < AgeTop[i]}}
#' \item{there are duplicate values in one of the columns \code{CmTop, CmBot}}
#' }
#' ...the function returns a fatal error and stops.
#'
#'
#' The function also fixes a few things:
#'
#' Add missing samples:
#' If any \code{CmBot[i] < CmTop[i+1]}, then
#' the function adds a new row \code{[j = i+1]}
#' to fill in the missing \code{CmTop} and \code{CmBot} values.
#' Thus, the added samples will have \code{CmTop[j] == CmBot[i]}
#' and \code{CmBot[j] == CmTop[i+1]}.
#' The new age scale values will be:
#' \code{AgeTop[j] == AgeBot[i]} and \code{AgeBot[j] == AgeTop[i+1]}.
#'
#' Check for overlapping sample depths:
#' If any \code{CmBot[i] > CmTop[i+1]}, the function flags a warning due to the
#' presence of overlapping sample depths. To get a data.frame of the
#' overlapping sample depths, use the \code{tapas::get_overlaps()} function.
#'
#' Check presence of slumps:
#' If after these checks any \code{AgeTop[i] == AgeBot[i]}, then
#' the function *removes the flagged rows*,
#' and *creates new* corrected \code{CmTop} and \code{CmBot} scales
#' that exclude slumps,
#' such that \code{CmBot[i] > CmTop[i]}
#' and \code{CmBot[i] == CmTop[i+1]}.
#'
#' @param series The input data frame. A matrix with the following columns:
#'    \code{CmTop, CmBot, AgeTop, AgeBot, Volume},
#'    and one or more columns with the data
#'    which should be resampled (variables).
#'
#' @return A data frame.
#'
#' @seealso [get_overlap_depths]
#'
#' @export
#'
#' @author Walter Finsinger
#'
#' @examples
#' co <- tapas::co_char_data
#' co_check <- tapas::check_pretreat(co)
#'
check_pretreat <- function(series) {

  ## Gather data --------------------------------------------------------------
  A <- series
  cm <- A[ ,1]
  cmB <- A[ ,2]
  ybp <- A[ ,3]
  ybpB <- A[ ,4]

  n_var <- dim(A)[2] - 5


  # Initial check ups ---------------------------------------------------------

  ### Check that all CmTop < CmBot --------------------------------------------
  if (any((cmB < cm))) {
    A[which(cmB < cm), ]
    stop('Fatal error: CmTop > CmBot in some of the samples')
  }

  ### Check that all AgeTop < AgeBot ------------------------------------------
  if (any((ybpB < ybp))) {
    A[which(ybpB < ybp), ]
    stop('Fatal error: AgeTop > AgeBot in some of the samples')
  }

  ### Check for duplicate sample depths ---------------------------------------
  if (any(duplicated(cm) == T)) {
    print('Fatal error: the following CmTop samples depths are duplicated')
    return(A[duplicated(cm) == T, ])
  }

  if (any(duplicated(cmB) == T)) {
    print('Fatal error: the following CmBot samples depths are duplicated')
    return(A[duplicated(cmB) == T, ])
  }


  ## Check if there are slump samples ----------------------------------------

  ## This could arise if there's a slump or a tephra, thus for sediment
  ## sections that were deposited in a very short time interval.
  ## The result will be that AgeBot[i] == AgeTop[i]

  # Check for presence of samples having AgeBot[i] == AgeTop[i]
  slumps <- ifelse(ybpB == ybp, yes = 1, no = 0)
  slumps_index <- which(slumps == 1)
  slumps_n <- sum(slumps)

  #### If slumps were found ---------------------------------------------------
  ## Create a new depth scale based on the sample thicknesses of samples that
  ## are not identified as slump samples
  if (!length(slumps_index) > 0) {
    print('No slump samples detected.')
    print('')
  } else {
    slumps_cm <- vector()

    # Quantify total thickness of slumps
    for (i in 1:length(slumps_index)) {
      slumps_cm[i] <- cmB[ slumps_index[i] ] - cm[ slumps_index[i] ]
    }
    slumps_tot_cm <- sum(slumps_cm)

    print(paste0('Warning: AgeTop = AgeBot for ',slumps_n,' samples, totaling ',slumps_tot_cm,'cm'))
    print('')

    # Remove samples identified as slumps (those with AgeBot[i] == AgeTop[i])
    A <- A[-c(slumps_index), ]

    # Build a new, corrected depth scale based on the sample thicknesses
    smpl_thickn <- A[ ,2] - A[ ,1]

    cm_corr <- cumsum(smpl_thickn)
    cmB_corr <- cm_corr + smpl_thickn

    A[ ,1] <- cm_corr
    A[ ,2] <- cmB_corr
  }


  ## Check if the depth scale is continuous (cmBot[i] = cmTop[i+1]) ----------

  ## Get difference between cmBot of sample[i] and cmTop of sample[i+1].
  ## The idea being that:
  ## - negative differences point to gaps, whereas
  ## - positive differences point to partially overlapping samples.

  ## First, redefine the depth scale as it may have changed in the previous
  ## loop:
  cm <- A[ ,1]
  cmB <- A[ ,2]
  gaps_index <- cmB[1:length(cmB) - 1] - cm[2:length(cm)]

  # Get indices for samples above missing samples
  gaps <- which(gaps_index < 0)

  if (length(gaps) > 0) { # if gaps were found...
    print('Warning: there are gaps in the depth scale,')
    print('a new depth scale will be created')
    print('')
  } else {
    print('No gaps found in the depth scale')
    print('')
  }

  # Get indices for overlapping samples
  overlaps <- which(gaps_index > 0)

  if (length(overlaps) > 0) { #if overlaps were found...
    print('Warning: overlapping depths found,')
    print('you may get them with the "tapas::get_overlap_samples() function')
    print('')
  } else {
    print('No overlapping depths found')
    print('')
  }




  ## Check if the age scale is continuous (AgeBot[i] = AgeTop[i+1]) ----------

  ## Difference between cmBot of sample[i] - cmTop of sample[i+1]
  gaps_age_index <- ybpB[1:length(ybpB) - 1] - ybp[2:length(ybp)]

  if (any(gaps_age_index != 0) == T) {
    print('Warning: the age scale is not continuous.')
  } else {
    print('The age scale is continuous')
    print('')
  }

  # Indices for samples above missing samples
  gaps_age <- which(gaps_age_index < 0)

  if (length(gaps_age) > 0) { # if gaps were found...
    print('Gaps found in the age scale. Missing samples will be added')
  }


  ## If age scale is continuous & depth scale is not --------------------------

  if (!length(gaps_age) > 0 && length(gaps) > 0) {

    # Build a new, corrected depth scale based on the sample thicknesses
    smpl_thickn <- A[ ,2] - A[ ,1]

    cm_corr <- cumsum(smpl_thickn) - smpl_thickn[1]
    cmB_corr <- c(cm_corr[2:length(cm_corr)],
                  cm_corr[length(cm_corr)] + smpl_thickn[length(smpl_thickn)])

    A[ ,1] <- cm_corr
    A[ ,2] <- cmB_corr
  }


  ## If both the age scale & the depth scale are not continuous ---------------

  if (length(gaps_age) > 0 && length(gaps) > 0) {

    ## Add rows for missing samples
    gap_cm <- vector()

    # Get gaps
    for (i in 1:length(gaps)) {
      gap_cm[i] <- cm[ gaps[i] + 1] - cmB[ gaps[i] ]
    }
    tot_gap <- sum(gap_cm)

    # print warning
    print(paste0("Warning: added ",length(gaps)," missing samples, totaling ",
                 tot_gap," cm."))


    # get params and series data for each of the detected gaps
    A_gaps <- A[0, ]

    for (i in 1:length(gaps)) {
      A_gaps[i, 1] <- cmB[ gaps[i] ]
      A_gaps[i, 2] <- cm[ gaps[i] + 1 ]
      A_gaps[i, 3] <- ybpB[ gaps[i] ]
      A_gaps[i, 4] <- ybp[ gaps[i] + 1 ]
      A_gaps[i, 5] <- 0
      A_gaps[i, c(6:(5 + n_var) )] <- NA
    }

    # add gaps to input data frame
    A <- rbind(A, A_gaps)
    A <- A[order(A[ ,1]), ]

    # If there are still any A_gaps$TopAge[i] == A_gaps$BotAge[i],
    # slightly modify the sample ages
    p <- which(A$TopAge == A$BotAge)
    A$BotAge[p] <- A$BotAge[p] - 0.01
    A$TopAge[p + 1] <- A$BotAge[p]
  }

  ## Prepare output ####
  output <- A

  return(output)
}
