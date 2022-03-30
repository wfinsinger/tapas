#' Pre-process data by resampling all columns.
#'
#' Uses paleofire::pretreatment() function
#' to virtually resample timeseries at equal sampling intervals.
#'
#' @param series A matrix with the following first columns:
#'               \code{c("CmTop", "CmBot", "AgeTop", "AgeBot", "Volume")}
#'               and additional columns with the data which should be resampled.
#' @param out Desired return value: \describe{
#'             \item{"accI"}{the function returns resampled accumulation rates}
#'             \item{"conI"}{the function returns resampled concentrations}
#'             \item{"countI"}{the function returns resampled counts}
#'            }
#' @param series.name A string defining the name of the input matrix.
#'                    Defaults to \code{NA}.
#' @param first,last Age boundaries of the resampled time serie.
#'                   If unspecified (\code{first=NULL} and \code{last=NULL}),
#'                   then resampling is done over the entire sequence,
#'                   from \code{min(series$AgeTop)} to max(series$AgeBot).
#' @param yrInterp Resolution of the resampled timeseries.
#'
#' @return A list with the resampled data according to \code{out} parameter.
#'
#' @importFrom stats median
#'
#' @export
pretreatment_data <- function(series=NULL, out="accI", series.name=NA,
                              first=NULL, last=NULL, yrInterp=NULL) {

  # series=co; out="accI"; series.name=NA;
  # first=NULL; last=NULL; yrInterp=NULL
  
  ## Gather data ####
  ybp <- series[ ,3]  # sample AgeTop
  ybpB <- series[ ,4] # sample AgeBot
  
  ## Initial check up ####
  if (is.null(first)) first <- min(ybp)
  if (is.null(last)) last <- max(ybpB)
  
  ## Calculate yrInterp ####
  if (is.null(yrInterp)) {
    yrInterp <- round(median(diff(ybp)))
  }
  
  
  ## [yr BP] Years to resample record to ####
  ybpI <- seq(from = first, to = last, by = yrInterp) 
  
  
  
  ## If output == resampled accumulation rates ####
  if (out == "accI") {
    raw <- data.frame(ybp)
    colnames(raw) <- "age"
    
    int <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4,
                             nrow = length(ybpI)))
    colnames(int) [1] <- "age"
    int$age <- ybpI
    
    conI <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4,
                              nrow = length(ybpI)))
    colnames(conI) [1] <- "age"
    conI$age <- ybpI
    
    countI <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4,
                                nrow = length(ybpI)))
    colnames(countI) [1] <- "age"
    countI$age <- ybpI
    
    j = 2
    for (i in 6:ncol(series)) {
      pre.i <- paleofire::pretreatment(params = series[ ,1:5],
                                       serie = series[ ,i], Int = T,
                                       first = first, last = last,
                                       yrInterp = yrInterp)
      
      raw[ ,j] <- pre.i$acc
      colnames(raw) [j] <- colnames(series) [i]
      
      int[ ,j] <- pre.i$accI
      colnames(int) [j] <- paste0(colnames(series) [i], "AR")
      
      conI[ ,j] <- pre.i$conI
      colnames(conI) [j] <- colnames(series) [i]
      
      countI[ ,j] <- pre.i$countI
      colnames(countI) [j] <- colnames(series) [i]
      
      volI <- pre.i$volI
      cmI <- pre.i$cmI
      
      j = j + 1
    }
  }
  
  
  ## If output = resampled concentration values ####
  if (out == "conI") {
    raw <- data.frame(ybp)
    colnames(raw) <- "age"
    
    int <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4,
                             nrow = length(ybpI)))
    colnames(int) [1] <- "age"
    int$age <- ybpI
    
    conI <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4,
                              nrow = length(ybpI)))
    colnames(conI) [1] <- "age"
    conI$age <- ybpI
    
    countI <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4,
                                nrow = length(ybpI)))
    colnames(countI) [1] <- "age"
    countI$age <- ybpI
    
    j = 2
    for (i in 6:ncol(series)) {
      pre.i <- paleofire::pretreatment(params = series[ ,1:5],
                                       serie = series[ ,i], Int = T,
                                       first = first, last = last,
                                       yrInterp = yrInterp)
      
      raw[ ,j] <- series[ ,i] / series[ ,5] # non-resampled concentration values
      colnames(raw) [j] <- colnames(series) [i]
      
      int[ ,j] <- pre.i$conI # resampled concentration values
      colnames(int) [j] <- colnames(series) [i]
      
      conI[ ,j] <- pre.i$conI # resampled concentration values
      colnames(conI) [j] <- colnames(series) [i]
      
      countI[ ,j] <- pre.i$countI
      colnames(countI) [j] <- colnames(series) [i]
      
      volI <- pre.i$volI
      cmI <- pre.i$cmI
      
      j = j + 1
    }
  }
  
  
  ## If output = resampled counts ####
  if (out == "countI") {
    raw <- data.frame(ybp)
    colnames(raw) <- "age"
    
    int <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4,
                             nrow = length(ybpI)))
    colnames(int) [1] <- "age"
    int$age <- ybpI
    
    conI <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4,
                              nrow = length(ybpI)))
    colnames(conI) [1] <- "age"
    conI$age <- ybpI
    
    countI <- data.frame(matrix(data = NA, ncol = dim(series)[2] - 4,
                                nrow = length(ybpI)))
    colnames(countI) [1] <- "age"
    countI$age <- ybpI
    
    j = 2
    for (i in 6:ncol(series)) {
      pre.i <- paleofire::pretreatment(params = series[ ,1:5],
                                       serie = series[ ,i], Int = T,
                                       first = first, last = last,
                                       yrInterp = yrInterp)
      
      raw[ ,j] <- series[ ,i]  # non-resampled count values
      colnames(raw) [j] <- colnames(series) [i]
      
      int[ ,j] <- pre.i$countI
      colnames(int) [j] <- colnames(series) [i]
      
      conI[ ,j] <- pre.i$conI
      colnames(conI) [j] <- colnames(series) [i]
      
      countI[ ,j] <- pre.i$countI
      colnames(countI) [j] <- colnames(series) [i]
      
      volI <- pre.i$volI
      cmI <- pre.i$cmI
      
      j = j + 1
    }
  }
  
  d.out <- structure(list(list(series = raw, series.name = series.name),
                          list(series.int = int, series.conI = conI,
                               series.countI = countI, cmI = cmI,
                               volI = volI, yr.interp = yrInterp,
                               type = "pretreatment"), out = out))
  names(d.out) [1] <- "raw"
  names(d.out) [2] <- "int"
  
  return(d.out)
}