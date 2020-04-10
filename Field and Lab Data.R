##### Setup #####

## load packages
library(plyr)
library(ggplot2)
library(dplyr)
library(ggridges)

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

summary(field_data)

#### Lab Data ####

## Import data

lab_data <- read.csv('NCAG_lab_data.csv', header = TRUE)

## Check

names(lab_data)

summary(lab_data)

#### Combine ####

animal_info <- left_join(lab_data,
                         field_data,
                         by = c('ID', 'site'))

animal_info$ID <- as.factor(animal_info$ID)

summary(animal_info)

summary(animal_info$sample.type)

animal_info$sample.type <- 
  mapvalues(animal_info$sample.type,
            from = c('Flatfish Guts',
                     'Rockfish Guts',
                     'Surfperch Guts'),
            to = c('Flatfish',
                   'Rockfish',
                   'Surfperch'))

ggplot(animal_info) +  # plot size distribution by species
  geom_density_ridges(aes(x = log(total.body.wet.weight),
                          y = reorder(species, total.body.wet.weight, mean)),
                      fill = 'lightblue') +
  labs(x = 'ln(Body Weight) (g)',
       y = 'Species') +
  theme_classic()
