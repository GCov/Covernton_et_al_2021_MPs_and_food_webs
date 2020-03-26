##### Project Setup #####

## load packages
library(plyr)
library(ggplot2)
library(dplyr)
library(ggridges)

## define standard error function

se <- function(x) {
  sd(x)/sqrt(length(x))
}



#### Field Data ####

## load field data

field_data <- read.csv("NCAG_field_data.csv", header = TRUE)

## Clean up data and merge

## Clean up field data

names(field_data)
summary(field_data$site)
summary(field_data$animal.type)
summary(field_data$species)
summary(field_data$sex)

field_data$sex <- mapvalues(field_data$sex, 
                            from = c('f',
                                     'F',
                                     'm',
                                     'M'),
                            to = c('Female',
                                   'Female',
                                   'Male',
                                   'Male'))

head(field_data)
