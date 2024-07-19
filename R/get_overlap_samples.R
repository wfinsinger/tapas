#' Get overlapping sample depths
#'
#' This function makes it easy to see which samples in the input file have
#' overlapping sample depths. It may be used if the
#' \code{tapas::check_pretreat()} function flags a warning for the presence of
#' overlapping sample depths.
#'
#' @details
#' The function checks if any \code{CmBot[i] > CmTop[i+1]}. If \code{TRUE},
#' a data.frame with the overalpping samples is returned.
#'
#' @param series The input data frame. A matrix with the following columns:
#'    \code{CmTop, CmBot, AgeTop, AgeBot, Volume},
#'    and one or more columns with the data
#'    that is being analysed (variables).
#'
#' @return A data frame including only the rows of the overlapping samples
#' sorted in ascending order by depth. The rownames are inherited from the
#' input data.frame, and mya thus help finding the flagged samples in raw data
#' files.
#'
#' @seealso [check_pretreat()]
#'
#' @importFrom dplyr distinct arrange pick
#'
#' @export
#'
#' @author Walter Finsinger
#'
get_overlap_depths <- function(series) {


  ## Gather data --------------------------------------------------------------
  A <- series
  cm <- A[ ,1]
  cmB <- A[ ,2]


  ## Check if the depth scale is continuous (cmBot[i] > cmTop[i+1]) ----------

  ## Get difference between cmBot of sample[i] and cmTop of sample[i+1].
  overlaps_index <- cmB[1:length(cmB) - 1] - cm[2:length(cm)]

  # Get indices for overlapping samples
  overlaps <- which(overlaps_index > 0)

  # Get overlapping samples from input data.frame:
  all_overlap_samples <- rbind(A[overlaps, ], A[overlaps + 1, ])

  # Remove duplicates and sort in ascending order by top depth
  overlap_samples <- dplyr::distinct(all_overlap_samples) %>%
    dplyr::arrange(dplyr::pick(1))

  # Return output
  return(overlap_samples)
}
