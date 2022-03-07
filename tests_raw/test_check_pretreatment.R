

## Load source files ####
source("../R/check_pretreatment.r")
source("./R/pretreatment_full.r")

## Load data ####
load("../data/co_char_data.rda")


## Shorten the data, such that plots are clearer ####
co <- co_char_data[1:20, ]


### Test with contiguous samples ####

co_i <- paleofire::pretreatment(params = co[ ,1:5],
                                serie = co[ ,6],
                                Int = T)
paleofire::plot.CHAR(co_i)



### Test when gaps are introduced ####
co_gaps <- co[-c(2:3, 12), ]

# The pretreatment() function fails because
# the sample depths and the sample ages are not continuous
co_gaps_i <- paleofire::pretreatment(params = co_gaps[ ,1:5],
                                     serie = co_gaps[ ,6],
                                     Int = T)
paleofire::plot.CHAR(co_gaps_i)

# Add missing samples, and use pretreatment
co_gaps_checked <- check_pretreat(series = co_gaps)
co_gaps_checked_i <- paleofire::pretreatment(params = co_gaps_checked[ ,1:5],
                                             serie = co_gaps_checked[ ,6],
                                             Int = T)

# Compare the output
par(mfrow = c(3,1), mar = c(3,5,1,1), oma = c(2,0,0,0))
paleofire::plot.CHAR(co_i)
paleofire::plot.CHAR(co_gaps_i)
paleofire::plot.CHAR(co_gaps_checked_i)



### Test with a 'slump' sample ####
### (AgeTop[i] == AgeBot[i])
co_slump <- co

## Introduce a slump in the nth sample
slump_sample <- 4
co_slump$AgeBot[slump_sample] <- co_slump$AgeTop[slump_sample]
co_slump$AgeTop[slump_sample + 1] <- co_slump$AgeTop[slump_sample]

### The pretreatment() function fails because AgeTop[i] == AgeBot[i]
co_slump_i <- paleofire::pretreatment(params = co_slump[ ,1:5],
                                      serie = co_slump[ ,6],
                                      Int = T)


co_slump_checked <- check_pretreat(series = co_slump)
co_slump_checked_i <- paleofire::pretreatment(params = co_slump_checked[ ,1:5],
                                              serie = co_slump_checked[ ,6],
                                              Int = T)

# Compare the output
par(mfrow = c(2,1), mar = c(3,5,1,1), oma = c(2,0,0,0))
paleofire::plot.CHAR(co_i)
paleofire::plot.CHAR(co_slump_checked_i)

