##### Setup #####

## load packages
library(plyr)
library(ggplot2)
library(dplyr)
library(ggridges)
library(colorspace)

## Define standard error function

se <- function(x) {
  sd(x)/sqrt(length(x))
}

#### Combine all data ####

## Combine lab and field data with isotopes data

foodweb1 <- 
  left_join(animal_info, 
            isotopes,
            by = c('ID', 'site', 'sample.type'))

## Combine all particle counts data into one dataframe

full_spec_data <- 
  rbind.data.frame(
    clams_polymer2,
    MU_polymer2,
    CR_polymer3[1:12],
    CU_polymer3[1:12],
    SS_polymer3[1:12],
    FF_polymer3[1:12],
    SP_polymer3[1:12],
    RF_polymer3[1:12]
  )

## Combine particle counts data with lab, field, and isotopes data

foodweb2 <- left_join(full_spec_data,
                      foodweb1,
                      by = c('ID', 'site', 'sample.type'))

foodweb2$ID <- as.factor(foodweb2$ID)
foodweb2$site <- as.factor(foodweb2$site)
foodweb2$sample.type <- as.factor(foodweb2$sample.type)

