

## Load raw data files
# From https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/charanalysis/sample_records_2012.zip
co_char_data <- read.csv("CO_charData.csv", header = T)
rand_peaks <- read.csv("randPeaks_charData.csv", header = T)
red_noise <- read.csv("redNoise_charData.csv", header = T)



## Change colnames

co_colnames <- c("cmTop", "cmBot", "AgeTop", "AgeBot", "charvol", "char")
col_names <- c("cm_top", "cm_bot", "age_top", "age_bot", "vol", "char")

colnames(co_char_data) <- co_colnames
colnames(rand_peaks) <- col_names
colnames(red_noise) <- col_names


## Add 'dummy_chararea' to data.frame 'co'
co_char_data$dummy_chararea <- co_char_data$char / 4.1


## Reduce length of the redNoise record
red_noise <- dplyr::filter(red_noise, age_bot <= 10000)


## Check records
# tapas::plot_raw(co_char_data)
# tapas::plot_raw(rand_peaks)
# tapas::plot_raw(red_noise)


# save(co, file = "co_char_data.rda")
# save(rand_peaks, file = "rand_peaks.rda")
# save(red_noise, file = "red_noise.rda")


usethis::use_data(co_char_data, overwrite = TRUE)
usethis::use_data(rand_peaks, overwrite = TRUE)
usethis::use_data(red_noise, overwrite = TRUE)


