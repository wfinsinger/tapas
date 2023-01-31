rm(list = ls())

library(tidyverse)
library(usethis)


## Load raw data --------------------------------------------------------------

# Source: https://github.com/wfinsinger/ARCO/blob/master/data_in/arco_Seedle.csv
if (!file.exists("./br_sdl.csv")) {
  download.file(
    "https://raw.githubusercontent.com/wfinsinger/ARCO/master/data_in/arco_Seedle.csv",
    "./br_sdl.csv")
}
br_sdl <- read_csv("./br_Sdl.csv", col_names = TRUE, col_types = "n")

## NOT RUN #####
# # Format data.frame for peak-identification analysis  -------------------------
# # Source: https://github.com/wfinsinger/ARCO/blob/master/data_in/arco_Smpl.csv
# download.file(
#   "https://raw.githubusercontent.com/wfinsinger/ARCO/master/data_in/arco_Smpl.csv",
#   "./br_smpl.csv")
# br_smpl <- read_csv("./br_smpl.csv", col_names = TRUE, col_types = "n")
#
# top_depth <- br_smpl$Depth
# bot_depth <- c(top_depth[2:length(top_depth)],
#                top_depth[length(top_depth)] + 1)
# top_age <- br_smpl$Age_calBP
# bot_age <- c(top_age[2:length(top_age)], top_age[length(top_age)] +
#                top_age[length(top_age)] - top_age[length(top_age) - 1])
# vol <- rep_len(1, length.out = length(top_depth))
# char_c <- br_smpl$SmplCount
# char_a <- br_smpl$SmplArea
#
# br_data <- as.data.frame(cbind(top_depth, bot_depth,
#                                top_age, bot_age, vol,
#                                char_c, char_a))
# write_csv(br_data, "./br_data.csv", col_names = TRUE)

br_data <- read_csv("./br_data.csv", col_names = TRUE, col_types = "n")

# Write output ----------------------------------------------------------------
usethis::use_data(br_sdl, overwrite = TRUE)
usethis::use_data(br_data, overwrite = TRUE)
