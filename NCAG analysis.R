## load packages
library(plyr)
library(ggplot2)
library(dplyr)

## load data

field_data <- read.csv("NCAG_field_data.csv", header = TRUE)
plankton_tows <- read.csv("plankton_tow_counts.csv", header = TRUE)

## Clean up data and merge

## Clean up field data

names(field_data)
summary(field_data$site)
summary(field_data$animal.type)
summary(field_data$species)

field_data$species <- mapvalues(field_data$species,
                                from = c('C. miniata',
                                         'M. magister',
                                         'P. californicus'),
                                to = c('Cucumaria miniata',
                                       'Metacarcinus magister',
                                       'Parastichopus californicus'))
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

## Clean up plankton tow data

names(plankton_tows)
summary(plankton_tows$size.fraction)
summary(plankton_tows$shape)
summary(plankton_tows$colour)

plankton_tows <- subset(plankton_tows, 
                        shape != 'fragment' |
                          colour != 'blue')

summary(plankton_tows$shape)

plankton_tows$shape <- mapvalues(plankton_tows$shape,
                                 from = c('fibre',
                                          'fibre ',
                                          'fragment'),
                                 to = c('Fibre',
                                        'Fibre',
                                        'Fragment'))
summary(plankton_tows$shape)

summary(plankton_tows$raman.ID)