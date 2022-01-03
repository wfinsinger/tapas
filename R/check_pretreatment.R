#  ***************************************************************************
#   check_pretreat.R
#  ---------------------
#   Date                 : January 2022
#   Copyright            : (C) 2022 by Walter Finsinger
#   Email                : walter.finsinger@umontpellier.fr
#  ---------------------
#
#  ***************************************************************************
#  
#### This functions does some check-up on the data frame used for
#### the paleofire::pretreatment() function, which requires that the depth and age scales
#### are continuous. In other words, for every i-th row, it requires that
#### - CmBot[i] > CmTop[i],
#### - AgeBot[i] > AgeTop[i],
#### - there shouldn't be duplicate values in the CmTop, CmBot, AgeTop, and AgeBot columns,
#### - CmBot[i] == CmTop[i+1],
#### - AgeBot[i] == AgeTop[i+1].
#### 
#### 
#### If:
#### - any CmBot[i] < CmTop[i], or
#### - any AgeBot[i] < AgeTop[i], or
#### - there are duplicate values in one of the columns CmTop, CmBot,
####    the function returns a fatal error and stops 
#### 
#### If any CmBot[i] < CmTop[i+1]:
#### - the function adds a new row [j = i+1] to fill in the missing CmTop and CmBot values.
####    Thus, the added samples will have CmTop[j] == CmBot[i] and CmBot[j] == CmTop[i+1].
####    The age scale values will be: AgeTop[j] == AgeBot[i] and AgeBot[j] == AgeTop[i+1].
####  
#### 
#### If after these checks any AgeTop[i] == AgeBot[i]:
#### - the function removes the flagged rows, and
#### - creates new corrected CmTop and CmBot scales that exclude slumps, such that
####    CmBot[i] > CmTop[i] and CmBot[i] == CmTop[i+1].
########################################################################################## 


check_pretreat <- function(series) {
  
  
  ## Gather data ####
  A <- series
  cm <- A[ ,1]
  cmB <- A[ ,2]
  ybp <- A[ ,3]
  ybpB <- A[ ,4]
  
  n_var <- dim(A)[2] - 5
  
  
  # Initial check ups ####
  
  # Check if CmTop < CmBot
  if (any((cmB < cm))) {
    A[which(cmB < cm), ]
    print('Fatal error: CmTop > CmBot in some of the samples')
    return()
  }
  
  # Check if AgeTop < AgeBot
  if (any((ybpB < ybp))) {
    A[which(ybpB < ybpB), ]
    print('Fatal error: AgeTop > AgeBot in some of the samples')
    return()
  }
  
  # Check if there are duplicate sample depths
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
      #i=1
      gap_cm[i] <- cm[ gaps[i] + 1] - cmB[ gaps[i] ]
    }
    tot_gap <- sum(gap_cm)
    
    # print warning
    print(paste0("Warning: added ",length(gaps)," missing samples, totaling ",tot_gap," cm."))
    
    
    # get params and series data for each of the detected gaps
    A_gaps <- A[0, ]
    
    for (i in 1:length(gaps)) {
      #i = 1
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
  
  rm(gap_cm, gaps, gaps_index)
  
  
  
  
  ## Check if the age scale is continuous ####
  ## 
  ## These could arise if there's a slump or a tephra, thus for sediment sections that
  ## were deposited in a very short time interval.
  ## The result will be that AgeBot[i] == AgeTop[i]
  
  slumps <- ifelse(ybpB == ybp, yes = 1, no = 0) # check for slumps 
  slumps_index <- which(slumps == 1)
  slumps_n <- sum(slumps)
  
  
  if (length(slumps_index) > 0) {
    slumps_cm <- vector()
    
    for (i in 1:length(slumps_index)) {
    slumps_cm[i] <- cmB[ slumps_index[i] ] - cm[ slumps_index[i] ]
    }
    slumps_tot_cm <- sum(slumps_cm)
    
    print(paste0('Warning: AgeTop = AgeBot for ',slumps_n,' samples, totaling ',slumps_tot_cm,'cm'))
    
    A <- A[-c(slumps_index), ]
    
    smpl_thickn <- A[ ,2] - A[ ,1]
    
    cm_corr <- cumsum(smpl_thickn)
    cmB_corr <- cm_corr + smpl_thickn
    
    A$cmTop <- cm_corr
    A$cmBot <- cmB_corr
    
  } else {
    print('No slumps detected. The age scale is continuous.')
  }
  
  
  ## Prepare output ####
  output <- A
  
  return(output)
  
}
