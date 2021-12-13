#  ***************************************************************************
#   pretreatment_data.R
#  ---------------------
#   Date                 : July 2020
#   Copyright            : (C) 2020 by Walter Finsinger
#   Email                : walter.finsinger@umontpellier.fr
#  ---------------------
#
#  ***************************************************************************
#
#   The script uses the pretreatment_full() function to virtually resample timeseries
#     at equal sampling intervals.
#
#   Requires a matrix as input with the following columns 
#   CmTop CmBot AgeTop AgeBot Volume 
#   and one or more columns with the data which should be resampled.
#
#   The user-defined parameters are as follows:
#     series      ->  the name of the input matrix
#     out         ->  with out = "accI" the function returns resampled accumulation rates,
#                     with out = "conI" the function returns resampled concentrations,
#                     with out = "countI" the function returns resampled counts.
#     series.name ->  a character string defining the name of the input matrix. Is NA by default.
#     first, last ->  determine the age boundaries of the resampled timeserie.
#                     If they are not specified (first & last == NULL),
#                     the resampling is done over the entire sequence,
#                     from min(series$AgeTop) to max(series$AgeBot).
#     yrInterp    ->  determines the resolution of the resampled timeseries.
#
#  ***************************************************************************

pretreatment_data <- function(series=NULL, out=NULL, series.name=NA, first=NULL, last=NULL,
                              yrInterp=yr.interp) {
  
  
  ## Initial check up ####
  if (is.null(first)) first <- min(series$AgeTop)
  if (is.null(last)) last <- max(series$AgeBot)
  
  ybpI <- seq(first, last, yrInterp) # [yr BP] Years to resample record to.
  
  
  ## If output == resampled accumulation rates ####
  if (out == "accI") {
    raw <- data.frame(series$AgeTop)
    colnames(raw) <- "age"
    
    int <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4, nrow = length(ybpI)))
    colnames(int) [1] <- "age"
    int$age <- ybpI
    
    conI <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4, nrow = length(ybpI)))
    colnames(conI) [1] <- "age"
    conI$age <- ybpI
    
    
    j = 2
    for (i in 6:ncol(series)) {
      #i=6
      pre.i <- pretreatment_full(params = series[ ,1:5], serie = series[ ,i], Int = F,
                                 first = first, last = last, yrInterp = yr.interp)
      
      raw[ ,j] <- pre.i$acc
      colnames(raw) [j] <- colnames(series) [i]
      
      int[ ,j] <- pre.i$accI
      colnames(int) [j] <- paste0(colnames(series) [i], "AR")
      
      conI[ ,j] <- pre.i$conI
      colnames(conI) [j] <- colnames(series) [i]
      
      volI <- pre.i$volI
      
      j = j + 1
    }
  }
  
  
  ## If output = resampled concentration values ####
  if (out == "conI") {
    raw <- data.frame(series$AgeTop)
    colnames(raw) <- "age"
    
    int <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4, nrow = length(ybpI)))
    colnames(int) [1] <- "age"
    int$age <- ybpI
    
    conI <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4, nrow = length(ybpI)))
    colnames(conI) [1] <- "age"
    conI$age <- ybpI
    
    j = 2
    for (i in 6:ncol(series)) {
      #i = 6
      pre.i <- pretreatment_full(params = series[ ,1:5], serie = series[ ,i], Int = F,
                                 first = first, last = last, yrInterp = yr.interp)
      
      raw[ ,j] <- series[ ,i] / series[ ,5] # non-resampled concentration values
      colnames(raw) [j] <- colnames(series) [i]
      
      int[ ,j] <- pre.i$conI # resampled concentration values
      colnames(int) [j] <- colnames(series) [i]
      
      conI[ ,j] <- pre.i$conI # resampled concentration values
      colnames(conI) [j] <- colnames(series) [i]
      
      volI <- pre.i$volI
      
      j = j + 1
    }
  }
  
  
  ## If output = resampled counts ####
  if (out == "countI") {
    raw <- data.frame(series$AgeTop)
    colnames(raw) <- "age"
    
    int <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4, nrow = length(ybpI)))
    colnames(int) [1] <- "age"
    int$age <- ybpI
    
    conI <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4, nrow = length(ybpI)))
    colnames(conI) [1] <- "age"
    conI$age <- ybpI
    
    j = 2
    for (i in 6:ncol(series)) {
      #i = 6
      pre.i <- pretreatment_full(params = series[ ,1:5], serie = series[ ,i], Int = F,
                                 first = first, last = last, yrInterp = yr.interp)
      
      raw[ ,j] <- series[ ,i]  # non-resampled count values
      colnames(raw) [j] <- colnames(series) [i]
      
      int[ ,j] <- pre.i$countI
      colnames(int) [j] <- colnames(series) [i]
      
      conI[ ,j] <- pre.i$conI
      colnames(conI) [j] <- colnames(series) [i]
      
      volI <- pre.i$volI
      
      j = j + 1
    }
  }
  
  d.out <- structure(list(list(series = raw, series.name = series.name),
                          list(series.int = int, series.conI = conI, volI = volI,
                               yr.interp = yrInterp,
                               span.l = NA, type = "pretreatment"), out = out))
  names(d.out) [1] <- "raw"
  names(d.out) [2] <- "int"
  
  return(d.out)
}
