## load packages
library(plyr)
library(ggplot2)
library(dplyr)

## load data

field_data <- read.csv("NCAG_field_data.csv", header = TRUE)
plankton_tows <- read.csv("plankton_tow_counts.csv", header = TRUE)
PT_field <- read.csv("PT_field_data.csv", header = TRUE)

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

plankton_tows$particle.type <- 
  mapvalues(plankton_tows$raman.ID,
            from = levels(plankton_tows$raman.ID),
            to = c('Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Unknown Anthropogenic',
                   'Natural',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Natural',
                   'Natural',
                   'Unknown'))
summary(plankton_tows$particle.type)

summary(plankton_tows$colour)

## Clean up plankton tow field data 

summary(PT_field$sample.volume)

PT <- left_join(plankton_tows, PT_field, by = 'ID')

PT$ID <- as.factor(PT$ID)

head(PT)

PT$num <- ifelse(is.na(PT$length), 0, 1)

## Separate blanks data

PT_blanks <- subset(PT, sample.type == 'Blanks')
PT <- subset(PT, sample.type == 'Plankton tows')

PT_polymer <- 
  PT %>% 
  group_by(ID, site, shape, particle.type, raman.ID, sample.volume) %>% 
  summarize(count = sum(num))

PT_polymer$concentration <- with(PT_polymer, count/sample.volume)*1000

PT_particle_type <- 
  PT %>% 
  group_by(ID, site, shape, sample.volume, particle.type) %>% 
  summarize(count = sum(num))

PT_particle_type$concentration <- 
  with(PT_particle_type, 
       count/sample.volume)*1000

PT_polymer_summary <-
  PT_polymer %>%
  group_by(site, shape, raman.ID) %>%
  summarize(mean = mean(concentration),
            sd = sd(concentration))

PT_particle_type_summary <-
  PT_particle_type %>% 
  group_by(site, shape, particle.type) %>% 
  summarize(mean = mean(concentration),
            sd = sd(concentration))
  
ggplot(PT_particle_type_summary) + 
  geom_col(aes(x = site,
               y = mean),
           colour = 'black',
           size = 1) + 
  labs(x = 'Particles/1000 L',
       y = 'Site') +
  facet_grid(shape ~ particle.type)
