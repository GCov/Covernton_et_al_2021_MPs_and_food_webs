## load packages
library(plyr)
library(ggplot2)
library(dplyr)

## load field data

field_data <- read.csv("NCAG_field_data.csv", header = TRUE)

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

## Load plankton tow data

plankton_tows <- read.csv("plankton_tows.csv", header = TRUE)
PT_field <- read.csv("PT_field_data.csv", header = TRUE)


## Clean up plankton tow data

names(plankton_tows)
summary(plankton_tows$size.fraction)
summary(plankton_tows$shape)
summary(plankton_tows$colour)

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

## Summarize blanks data

summary(PT_blanks)
PT_blanks_particle_type <- 
  PT_blanks %>% 
  group_by(ID, shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.count = sum(num))

PT_blanks_means <- 
  PT_blanks_particle_type %>% 
  group_by(shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.mean = mean(blank.count))

## Summarize PT data

PT_polymer <- 
  PT %>% 
  group_by(ID, site, shape, colour, blank.match, particle.type, raman.ID, 
           sample.volume) %>% 
  summarize(count = sum(num))

## Blank subtract

PT_polymer2 <-
  left_join(PT_polymer, 
            PT_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'raman.ID', 
                   'particle.type'))

PT_polymer2$blank.mean[is.na(PT_polymer2$blank.mean)] <- 0

PT_polymer2$adj.count <- with(PT_polymer2, count - blank.mean)

PT_particle_type <-
  PT_polymer2 %>%
  group_by(ID, site, particle.type, sample.volume) %>%
  summarize(count = sum(adj.count))

PT_synthetic <- subset(PT_particle_type,
                       particle.type == 'Synthetic Polymer')

PT_synthetic_summary <-
  PT_synthetic %>% 
  group_by(site) %>% 
  summarize(mean = mean(count/sample.volume),
            sd = sd(count/sample.volume))
  
ggplot(PT_synthetic_summary) + 
  geom_point(aes(x = site,
               y = mean),
           colour = 'black',
           size = 2) + 
  geom_errorbar(aes(x = site,
                    ymin = mean - sd,
                    ymax = mean + sd),
                size = 1) +
  labs(x = 'Particles/L',
       y = 'Site')

## Do all of the above for the plankon jar samples

## Load plankton tow data

plankton_jars <- read.csv("plankton_jars.csv", header = TRUE)

## Clean up plankton tow data

names(plankton_jars)
summary(plankton_jars$size.fraction)
summary(plankton_jars$shape)
summary(plankton_jars$colour)

summary(plankton_jars$shape)

plankton_jars$shape <- mapvalues(plankton_jars$shape,
                                 from = c('fibre',
                                          'fragment'),
                                 to = c('Fibre',
                                        'Fragment'))
summary(plankton_jars$shape)

summary(plankton_jars$raman.ID)

plankton_jars$particle.type <- 
  mapvalues(plankton_jars$raman.ID,
            from = levels(plankton_jars$raman.ID),
            to = c('Unknown',
                   'Natural Anthropogenic',
                   'Unknown Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(plankton_jars$particle.type)

summary(plankton_jars$colour)

plankton_jars$num <- with(plankton_jars,
                          ifelse(length == 'na', 0, 1))

## Separate blanks data

PJ_blanks <- subset(plankton_jars, sample.type == 'Blanks')
PJ <- subset(plankton_jars, sample.type == 'Plankton Jars')

## Summarize blanks data

summary(PJ_blanks)
PJ_blanks_particle_type <- 
  PJ_blanks %>% 
  group_by(ID, shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.count = sum(num))

PJ_blanks_means <- 
  PJ_blanks_particle_type %>% 
  group_by(shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.mean = mean(blank.count))

## Summarize PJ data

PJ_polymer <- 
  PJ %>% 
  group_by(ID, site, shape, colour, blank.match, particle.type, raman.ID) %>% 
  summarize(count = sum(num))

## Blank subtract

PJ_polymer2 <-
  left_join(PJ_polymer, 
            PJ_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'raman.ID', 
                   'particle.type'))

PJ_polymer2$blank.mean[is.na(PJ_polymer2$blank.mean)] <- 0

PJ_polymer2$adj.count <- with(PJ_polymer2, count - blank.mean)

PJ_particle_type <-
  PJ_polymer2 %>%
  group_by(ID, site, particle.type) %>%
  summarize(count = sum(adj.count))

PJ_synthetic <- subset(PJ_particle_type,
                       particle.type == 'Synthetic Polymer')

ID <- c('CBPJ1', 'CBPJ3', 'CBPJ4', 
        'EBPJ1', 'EBPJ2', 'EBPJ3', 
        'HPPJ2', 'HPPJ3', 'HPPJ5')
site <- as.factor(c(
  rep("Cole's Bay", 3),
  rep('Elliot Bay', 3),
  rep('Victoria Harbour', 3)
))
PJ_IDs <- data.frame(ID, site)

PJ_synthetic <- left_join(PJ_IDs, PJ_synthetic, by = c('ID', 'site'))

PJ_synthetic$count[is.na(PJ_synthetic$count)] <- 0

PJ_synthetic_summary <-
  PJ_synthetic %>% 
  group_by(site) %>% 
  summarize(mean = mean(count),
            sd = sd(count))
            
ggplot(PJ_synthetic_summary) + 
  geom_point(aes(x = site,
                 y = mean),
             colour = 'black',
             size = 2) + 
  geom_errorbar(aes(x = site,
                    ymin = mean - sd,
                    ymax = mean + sd),
                size = 1) +
  labs(x = 'Site',
       y = 'Particles/L')
