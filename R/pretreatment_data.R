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
#     series      ->  the matrix
#     out         ->  with out="accI" the function returns resampled accumulation rates,
#                     with out="conI" the function returns resampled concentrations.
#     series.name ->  a character string defining the name of the input matrix. Is NA by default.
#     first, last ->  determine the age boundaries of the resampled timeserie.
#                     If they are not specified (first & last == NULL), the resampling is done over the entire
#                     sequence (from min(series$AgeTop) to max(series$AgeBot)).
#     yrInterp    ->  determines the resolution of the resampled timeseries.
#
#  ***************************************************************************

pretreatment_data <- function(series=NULL, out=NULL, series.name=NA, first=NULL, last=NULL,
                              yrInterp=yr.interp) {
  
  
  if (is.null(first)) first <- min(series$AgeTop)
  if (is.null(last)) last <- max(series$AgeBot)
  
  ybpI <- seq(first, last, yrInterp) # [yr BP] Years to resample record to.
  
  
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
      #i=7
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
  
  if (out == "conI") {
    raw <- data.frame(series$AgeTop)
    colnames(raw) <- "age"
    
    int <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4, nrow = length(ybpI)))
    colnames(int) [1] <- "age"
    int$age <- ybpI
    
    j = 2
    for (i in 6:ncol(series)) {
      #i = 6
      pre.i <- pretreatment_full(params = series[ ,1:5], serie = series[ ,i], Int = F,
                                 first = first, last = last, yrInterp = yr.interp)
      
      raw[ ,j] <- series[ ,i]
      colnames(raw) [j] <- colnames(series) [i]
      
      int[ ,j] <- pre.i$countI
      colnames(int) [j] <- colnames(series) [i]
      
      volI <- pre.i$volI
      
      j = j + 1
    }
  }
  
  d.out <- structure(list(list(series = raw, series.name = series.name),
                          list(series.int = int, series.conI = conI, volI = volI,
                               yr.interp = yrInterp,
                               span.l = NA, type = "pretreatment")))
  names(d.out) [1] <- "raw"
  names(d.out) [2] <- "int"
  
  return(d.out)
}
