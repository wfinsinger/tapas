#' Check data frame used for \code{paleofire::pretreatment()} function.
#'
#' This functions does some check-up on the data frame used for
#' the \code{paleofire::pretreatment()} function,
#' which requires that the depth and age scales are continuous.
#' In other words, for every i-th row, it requires that \itemize{
#' \item{\code{CmBot[i] > CmTop[i]}}
#' \item{\code{AgeBot[i] > AgeTop[i]}}
#' \item{there shouldn't be duplicate values
#'       in the \code{CmTop, CmBot, AgeTop}, and \code{AgeBot} columns}
#' \item{\code{CmBot[i] == CmTop[i+1]}}
#' \item{\code{AgeBot[i] == AgeTop[i+1]}}
#' }
#'
#' If any of the following is true: \itemize{
#' \item{any \code{CmBot[i] < CmTop[i]}}
#' \item{\code{AgeBot[i] < AgeTop[i]}}
#' \item{there are duplicate values in one of the columns \code{CmTop, CmBot}}
#' }
#' .. then the function returns a fatal error and stops.
#'
#' The function also fixes a few things:
#'
#' If any \code{CmBot[i] < CmTop[i+1]}, then
#' the function adds a new row \code{[j = i+1]}
#' to fill in the missing \code{CmTop} and \code{CmBot} values.
#' Thus, the added samples will have \code{CmTop[j] == CmBot[i]}
#' and \code{CmBot[j] == CmTop[i+1]}.
#' The age scale values will be:
#' \code{AgeTop[j] == AgeBot[i]}
#' and \code{AgeBot[j] == AgeTop[i+1]}.
#'
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
#' @export
check_pretreat <- function(series) {
  
  
  ## Gather data ####
  A <- series
  cm <- A[ ,1]
  cmB <- A[ ,2]
  ybp <- A[ ,3]
  ybpB <- A[ ,4]
  
  n_var <- dim(A)[2] - 5
  
  
  # Initial check ups ####
  
  # Check that CmTop < CmBot
  if (any((cmB < cm))) {
    A[which(cmB < cm), ]
    stop('Fatal error: CmTop > CmBot in some of the samples')
  }
  
  # Check that AgeTop < AgeBot
  if (any((ybpB < ybp))) {
    A[which(ybpB < ybp), ]
    stop('Fatal error: AgeTop > AgeBot in some of the samples')
  }
  
  # Check for duplicate sample depths
  if (any(duplicated(cm) == T)) {
    print('Fatal error: the following CmTop samples depths are duplicated')
    return(A[duplicated(cm) == T, ])
  }
  
  if (any(duplicated(cmB) == T)) {
    print('Fatal error: the following CmBot samples depths are duplicated')
    return(A[duplicated(cmB) == T, ])
  }
  
  
  ## Check if the depth scale is continuous (cmBot[i] = cmTop[i+1]) ####
  
  ## Difference between cmBot of sample[i] - cmTop of sample[i+1]
  gaps_index <- cmB[1:length(cmB) - 1] - cm[2:length(cm)]
  
  # Indices for samples above missing samples
  gaps <- which(gaps_index < 0) 
  gap_cm <- vector()
  
  if (length(gaps) > 0) { # if gaps were found...
    
    # Get gaps
    for (i in 1:length(gaps)) {
      gap_cm[i] <- cm[ gaps[i] + 1] - cmB[ gaps[i] ]
    }
    tot_gap <- sum(gap_cm)
    
    # print warning
    print(paste0("Warning: added ",length(gaps)," missing samples, totaling ",tot_gap," cm."))
    
    
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
    A <- A[order(A$cmTop), ]
  } else {
    print("No missing samples detected.")
  }
  
  
  ## Check if the age scale is continuous ####
  ## 
  ## This could arise if there's a slump or a tephra, thus for sediment sections that
  ## were deposited in a very short time interval.
  ## The result will be that AgeBot[i] == AgeTop[i]
  
  # Check for presence of samples having AgeBot[i] == AgeTop[i]
  slumps <- ifelse(ybpB == ybp, yes = 1, no = 0)
  slumps_index <- which(slumps == 1)
  slumps_n <- sum(slumps)
  
  # If slumps were found:
  if (length(slumps_index) > 0) {
    slumps_cm <- vector()
    
    # Quantify total thickness of slumps
    for (i in 1:length(slumps_index)) {
    slumps_cm[i] <- cmB[ slumps_index[i] ] - cm[ slumps_index[i] ]
    }
    slumps_tot_cm <- sum(slumps_cm)
    
    print(paste0('Warning: AgeTop = AgeBot for ',slumps_n,' samples, totaling ',slumps_tot_cm,'cm'))
    
    # Remove samples identified as slumps (those with AgeBot[i] == AgeTop[i])
    A <- A[-c(slumps_index), ]
    
    # Build a new, corrected depth scale based on the sample thicknesses
    smpl_thickn <- A[ ,2] - A[ ,1]
    
    cm_corr <- cumsum(smpl_thickn)
    cmB_corr <- cm_corr + smpl_thickn
    
    A$CmTop <- cm_corr
    A$CmBot <- cmB_corr
    
  } else {
    print('No slumps detected; the age scale is continuous.')
  }
  
  
  ## Prepare output ####
  output <- A
  
  return(output)
  
}
