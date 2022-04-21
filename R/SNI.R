#' Calculate signal-to-noise index (SNI) for paleoecological timeseries records.
#'
#' This function computes SNI as described in Kelly et al. 2011.
#' Note that your data must be interpolated
#' to constant sample resolution (yr/sample)
#' before input to the function.
#' The function makes no assumption about prior analysis on the input series,
#' i.e. any background and threshold methods may be used.
#' However, input data should still be consistent
#' with the interpretation that a value (column 2)
#' greater than the corresponding threshold value (column 3)
#' is a *signal* sample,
#' whereas a value below the threshold is *noise*.
#' Refer to Kelly et al. 2010 for details and discussion.
#'
#' @author Ryan Kelly (University of Illinois, USA)
#' @references Supplementary materials to the publication:
#'             Kelly, R., Higuera, P., Barrett, C., & Hu, F. (2011).
#'             A signal-to-noise index to quantify
#'             the potential for peak detection in sediment–charcoal records.
#'             Quaternary Research, 75(1), 11-17.
#'             \url{https://doi.org/10.1016/j.yqres.2010.07.011}
#' @seealso Supplementary materials:
#'          \url{https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0033589400006785/resource/name/S0033589400006785sup001.txt}
#'
#' @param ProxyData Matrix of input data with one row per sample,
#'                  containing: \describe{
#'        \item{Column 1}{age associated with the sample (yr)}
#'        \item{Column 2}{proxy accumulation rate (e.g. CHAR)
#'                        of the sample (pieces/cm^2/yr)}
#'        \item{Column 3}{threshold value (pieces/cm^2/yr)}
#'        }
#' @param BandWidth Width of moving window for computing SNI.
#'
#' @return The function was modified by Walter Finsinger
#'         to return only a subset
#'         of the original \code{SNI_output} list: \describe{
#'           \item{\code{$SNI_sm}}{Lowess-smoothed SNI values.}
#'           \item{\code{$SNI_raw}}{Raw SNI values.}
#'         }
#'         In the original script,
#'         the returned data list contained the computed SNI
#'         and related data,
#'         with one row for each row in the input variable (ProxyData):
#'         \describe{
#'  \item{\code{$SNI_sm}}{The smoothed SNI computed for each sample.}
#'  \item{\code{$winInd}}{Indexes of the first and last samples
#'                        included in each moving window.
#'                        E.g. \code{SNI_output$winInd(X) == [A, B]}
#'                        indicates that the moving window used
#'                        to calculate SNI for the \code{X}th sample
#'                        contained all samples between \code{A} and \code{B},
#'                        inclusive.}
#'  \item{\code{$popN}}{The CHAR values
#'                      of all samples in the noise (N) population
#'                      (samples in the moving window
#'                      with CHAR below threshold)}
#'  \item{\code{$popS}}{The CHAR values
#'                      of all samples in the signal (S) population
#'                      (samples in the moving window
#'                      with CHAR at or above threshold)}
#'  \item{\code{$meanN}}{Mean CHAR of the samples in popN}
#'  \item{\code{$stdN}}{standard deviation of the samples in popN}
#'  \item{\code{$CF}}{The *correction factor*
#'                    used in computing SNI.
#'                    Equal to (v - 2)/v,
#'                    where v is the number of samples in popN}
#'  }
#'
#' @importFrom stats sd
#'
#' @export
SNI <- function(ProxyData, BandWidth) {

  # Data setup
  ages <- ProxyData[ ,1]
  CHAR <- ProxyData[ ,2]
  thresh <- ProxyData[ ,3]

  r <- mean(diff(ages)) # r = mean sampling resolution

  if (!sd(diff(ages)) == 0) { # If resolution (yr/sample) is not constant
    print(paste0('Warning: Sampling resolution (mean: ',r,
                ' yr/sample) is not constant (sd = ', sd(diff(ages)),').'))
  }


  # Preallocate some space
  SNI_output 				<- list()
  #SNI_output$SNI 		<- NA * numeric(length(ages))
  SNI_output$winInd	<- matrix(data = NA, nrow = length(ages), ncol = 2)
  SNI_output$popN 	<- vector(mode = "list", length = length(ages))
  SNI_output$popS 	<- vector(mode = "list", length = length(ages))
  SNI_output$meanN 	<- NA * numeric(length(ages))
  SNI_output$stdN 	<- NA * numeric(length(ages))
  SNI_output$CF 		<- NA * numeric(length(ages))
  rawSNI 						<- NA * numeric(length(ages))

  for (i in 1:length(ages)) { # Perform calculations at each sample age
    # Calculate window indexes
    if ( i < round(0.5*(BandWidth/r)) + 1 ) {  # Samples near beginning (moving window truncated)
      SNI_output$winInd[i, ] <- c(1, (i + round(0.5*(BandWidth/r))))
    } else {
      if ( i > length(ages) - round(0.5*(BandWidth/r)) ) { # Samples near end
        SNI_output$winInd[i, ] <- c((i - round(0.5*(BandWidth/r))),
                                    length(ages) )
      } else {
        SNI_output$winInd[i, ] <- c((i - round(0.5*(BandWidth/r))),
                                    (i + round(0.5*(BandWidth/r))))
      }
    }

    # In each moving window...

    # Get CHAR for samples in window (X is entire window population)
    X <- CHAR[SNI_output$winInd[i,1]:SNI_output$winInd[i,2]]

    # inS, inN: boolean indexes of samples counted as S & N, respectively
    inS <- (X >= thresh[SNI_output$winInd[i,1]:SNI_output$winInd[i,2]])
    inN <- !inS

    # Fill in S & N populations, means, stds
    SNI_output$popS[[i]] <- X[inS]
    SNI_output$popN[[i]] <- X[inN]
    SNI_output$meanN[i] <- mean(X[inN])
    SNI_output$stdN[i] <- sd(X[inN])

    # Calculate correction factor
    v <- length(SNI_output$popN[[i]]) # d.f. of N
    SNI_output$CF[i] <- (v - 2) / v # SNI = Z * (v-2)/v

    # Calculate raw SNI (unsmoothed)
    if ( length(SNI_output$popS[[i]]) > 0) {
      rawSNI[i] <- (mean( (SNI_output$popS[[i]] - SNI_output$meanN[i]) )
                    / SNI_output$stdN[i]) * SNI_output$CF[i]
    } else {
      rawSNI[i] <- 0 # SNI = 0 by definition when no samples exceed threshold
    }

  }

  ## Replace NA and Inf values with zero (added for tapas::)
  rawSNI[is.na(rawSNI)] <- 0
  rawSNI[is.infinite(rawSNI)] <- 0

  # Smooth raw values to obtain final SNI
  SNI_sm <- lowess(ages, rawSNI, f = (BandWidth/r)/length(ages), iter = 0)$y
  SNI_raw <- rawSNI

  # Prepare output
  SNI_output <- as.data.frame(cbind(SNI_raw, SNI_sm))

  return(SNI_output)

}
