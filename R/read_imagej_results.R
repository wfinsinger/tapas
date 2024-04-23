#' Reads results of charcoal analyses made with ImageJ
#'
#' @param folder_in Sets the path to the directory where the *.csv files are
#'                  located (e.g. "./data_in/" if the files are in a subfolder
#'                  relatively to the working directory).
#' @param metadata A vector that describes the metadata stored in the file names.
#'                 By default, metadata = c("core_id", "depth", "image_id"),
#'                 for file names formatted as "CoreID_Depth_ImageID.csv".
#'                 However, other file name formats are possible. See the
#'                 Description here below.
#' @param image_id_default Specifies the value or character that shall be
#'                          added if sometimes the file name does not include
#'                          the information for the ImageID. By default,
#'                          \code{image_id_default = "a"}.
#' @param na_replace Logical. If \code{na_replace = FALSE} (default), \code{NA}
#'                    values will not be replaced with \code{0} values for all
#'                    numeric columns.
#'
#' @description
#' The function loads spreadsheets obtained with ImageJ (https://imagej.net/)
#' that were stored as *.csv files (comma-separated) in a common folder whose
#' path is specified by the \code{folder_in} argument.
#'
#' @details
#' The function expects that the file names include a set of information that
#' identify the samples unequivocally, and that the information is given as
#' strings separated by a dash (-) or an underscore (_) (do not use empty
#' spaces (space bar) to separate the metadata).
#' It is recommended to use a hierarchical file-name format that includes the
#' ImageID information last.
#'
#' By default the following set of information is expected in the file names:
#'  - a site identifier;
#'  - a sample identifier, for instance the master core depth;
#'  - an image identifier (a character string) in case several images were
#'     taken for a sample. If only one image was taken from a sample, the
#'     image identifier is optional (an "a" will be added if the file name does
#'     not include this information).
#'
#' However, users can choose a different set and number of identifiers.
#' Possible alternatives may include:
#' - "CoreID_SectionID_Depth_ImageID.csv" or
#' - "CoreID_SectionID_DepthInDrive_ImageID.csv".
#'
#' If the file names include a "depth" value, make sure to take note as to
#' which depth values are mentioned (top depth, mid depth, bottom depth), as
#' well as to which depth scale the values are referring to (e.g. the master
#' core depth, depth in drive, or other).
#' Further, make sure the file names include useful information to link
#' the data with other data types from the same sediment archive,
#' such as the chronology or other proxy data.
#'
#' @returns
#' A data.frame with the specified \code{metadata} columns, and all columns
#' that appear in any of the input files.
#'
#' @author Walter Finsinger
#'
#' @examples
#' \dontrun{
#' folder_in = "./data_in/",
#' metadata = c("site_id", "depth", "image_id"),
#' image_id_default = "a",
#' na_replace = FALSE
#' }
#'
#' @export
#'
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#' @importFrom utils read.csv
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyr replace_na
#' @importFrom tidyselect where

read_imagej_results <- function(folder_in = NULL,
                                metadata = c("site_id", "depth", "image_id"),
                                image_id_default = "a", na_replace = FALSE) {

  ## Read file names in the 'folder_in'
  file_names <- list.files(folder_in)
  n_files <- length(file_names)

  ## Loop to read the *.csv files, split the file names, and add columns with
  ## metadata written in the file names:
  df_list <- list()

  for (i in 1:n_files) {
    # i = 1
    ## Get the metadata from the file name:
    f_ids_i <- unlist(strsplit(file_names[i], split = '[. _ -]+'))
    f_ids_i <- f_ids_i[!f_ids_i == "csv"]

    # Check if file name is compatible with the desired metadata:
    if (!length(f_ids_i) == length(metadata) &&
        !length(f_ids_i) == length(metadata) - 1) {
      stop(print(
        paste("The filename", file_names[i],
              "is incompatible with the desired set of metadata")))
    }

    # If the file name passes the check-up, write metadata in metadata_i:
    if (length(f_ids_i) == length(metadata) - 1) {
      f_ids_i[length(f_ids_i) + 1] <- image_id_default
    }

    metadata_i <- matrix(NA, ncol = length(metadata))
    colnames(metadata_i) <- metadata

    for (j in 1:length(metadata)) {
      if (metadata[j] == "depth") {
        metadata_i[j] <- as.numeric(f_ids_i[j])
      } else {
        metadata_i[j] <- f_ids_i[j]
      }
    }

    ## Get data from the .csv file:
    f_data_i <- read.csv(paste0(folder_in, file_names[i]))
    colnames(f_data_i)[1] <- "char_id"

    ## If a csv has no data (for samples without charcoal pieces),
    ## add an NA row:
    if (dim(f_data_i)[1] == 0) {
      f_data_i[1, ] <- NA
    }

    ## Bind metadata to the loaded data.frame, and move it to the list:
    f_data_i <- dplyr::bind_cols(metadata_i, f_data_i)
    df_list[[i]] <- f_data_i
  }

  ## Gather data from list of data.frames
  df <- dplyr::bind_rows(df_list)

  ## If one wants to replace the NA values
  if (na_replace == TRUE) {
    df <- df %>% dplyr::mutate(dplyr::across(tidyselect::where(is.numeric),
                                      ~tidyr::replace_na(.x, 0)))
  }

  ## Return output:
  return(df)
}
