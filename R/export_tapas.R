#' Prepares an output data frame to summarize the results of the peak-detection analysis.
#'
#' @description
#' Requires an output from one of the following functions:
#' \code{\link{global_thresh}()}, \code{\link{local_thresh}()},
#' or \code{\link{peak_detection}()}.
#'
#' @param series The output of \code{\link{global_thresh}()},
#'               \code{\link{local_thresh}()},
#'               or \code{\link{peak_detection}()}.
#'
#' @return A data frame.
#'
#' @examples
#' \dontrun{
#' co <- tapas::co_char_data
#' co_loc <- tapas::peak_detection(co, proxy = "char")
#' co_glob_exp <- tapas::tapas_export(co_loc)
#' }
#'
#' @importFrom stringr str_remove
#'
#' @importFrom dplyr select pull left_join bind_cols all_of
#'
#' @export
#'
#' @author Walter Finsinger
#'
tapas_export <- function(series = NULL) {

  ## extract data from the list
  cm_i <- series$int$cmI
  age_top_i <- series$int$series.int$age
  vol_i <- series$int$volI

  want_proxy <- series$thresh$proxy

  if (series$out == "accI") {
    want_proxy_conc <- stringr::str_remove(want_proxy, "AR")
    con_i <- dplyr::pull(series$int$series.conI, all_of(want_proxy_conc))
    count_i <- dplyr::pull(series$int$series.countI, all_of(want_proxy_conc))
  } else {
    con_i <- dplyr::pull(series$int$series.conI, all_of(want_proxy))
    count_i <- dplyr::pull(series$int$series.countI, all_of(want_proxy))
  }

  acc_i <- dplyr::pull(series$int$series.int, want_proxy)
  bkg_trend <- dplyr::pull(series$detr$detr, want_proxy)
  sni_raw <- series$thresh$SNI_pos$SNI_raw
  sni_smooth <- series$thresh$SNI_pos$SNI_sm
  peaks_pos_sig <- series$thresh$peaks.pos[ ,1]
  peaks_pos_insig <- series$thresh$peaks.pos.insig


  if (length(series$thresh$thresh.pos) == 1) {
    thresh_final_pos <- rep_len(series$thresh$thresh.pos,
                                length.out = length(age_top_i))
  } else {
    thresh_final_pos <- series$thresh$thresh.pos[ ,1]
  }

  thresh_final_pos_sm <- series$thresh$thresh.pos.sm

  if (length(thresh_final_pos_sm) == 1) {
    thresh_final_pos_sm <- rep_len(thresh_final_pos_sm,
                                   length.out = length(age_top_i))
  }

  if (length(series$thresh$thresh.neg) == 1) {
    thresh_final_neg <- rep_len(series$thresh$thresh.neg,
                                length.out = length(age_top_i))
  } else {
    thresh_final_neg <- series$thresh$thresh.neg[ ,1]
  }

  thresh_final_neg_sm <- series$thresh$thresh.neg.sm
  if (length(thresh_final_neg_sm) == 1) {
    thresh_final_neg_sm <- rep_len(thresh_final_neg_sm,
                                   length.out = length(age_top_i))
  }

  ## Gather output so far
  output <- cbind(cm_i, age_top_i, vol_i, count_i, con_i,
                  acc_i, bkg_trend,
                  thresh_final_pos, thresh_final_pos_sm,
                  sni_raw, sni_smooth,
                  peaks_pos_insig, peaks_pos_sig)
  output <- as.data.frame(output)

  ## Add PeakMag and Return Intervals
  peak_mag <- series$thresh$PeakMag_pos
  colnames(peak_mag) <- c("age_top_i", "peak_mag")
  output <- dplyr::left_join(output, peak_mag,
                             by = "age_top_i")
  output$peak_mag[is.na(output$peak_mag)] <- 0


  peaks_pos_ages <- series$thresh$peaks.pos.ages
  RI_pos <- series$thresh$RI_pos

  RI_pos <- as.data.frame(cbind(peaks_pos_ages, RI_pos))
  output <- dplyr::left_join(output, RI_pos,
                             by = c("age_top_i" = "peaks_pos_ages"))

  return(output)
}
