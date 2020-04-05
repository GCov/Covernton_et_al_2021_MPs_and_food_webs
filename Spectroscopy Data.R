##### Setup #####

## load packages
library(plyr)
library(ggplot2)
library(dplyr)

#### Plankton Tows ####

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

plankton_tows$particle.type[plankton_tows$shape == 'Fibre' &
                              plankton_tows$colour == 'clear' &
                              plankton_tows$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

## Clean up plankton tow field data 

summary(PT_field$sample.volume)

PT <- left_join(plankton_tows, PT_field, by = 'ID')

PT$ID <- as.factor(PT$ID)

head(PT)

PT$num <- 
  ifelse(is.na(PT$length), 0, 1)  # Separate out samples that had 0 counts

## Separate blanks data

PT_blanks <- subset(PT, sample.type == 'Blanks')
PT <- subset(PT, sample.type == 'Plankton Tows')

## Summarize blanks data

summary(PT_blanks)
PT_blanks_particle_type <- 
  PT_blanks %>% 
  group_by(ID, shape, colour, size.fraction, blank.match, raman.ID, 
           particle.type) %>% 
  summarize(blank.count = sum(num))

PT_blanks_means <- 
  PT_blanks_particle_type %>% 
  group_by(shape, colour, size.fraction, blank.match, raman.ID, 
           particle.type) %>% 
  summarize(blank.mean = mean(blank.count))

## Summarize PT data

PT_polymer <- 
  PT %>% 
  group_by(ID, site, sample.type, size.fraction, shape, colour, blank.match, 
           particle.type, raman.ID, sample.volume) %>% 
  summarize(count = sum(num))

## Blank subtract

PT_polymer2 <-
  left_join(PT_polymer, 
            PT_blanks_means, 
            by = c('shape', 'colour', 'size.fraction', 'blank.match', 
                   'raman.ID', 'particle.type'))

PT_polymer2$blank.mean[is.na(PT_polymer2$blank.mean)] <- 0

PT_polymer2$adj.count <- with(PT_polymer2, count - blank.mean)

PT_polymer2 %>% 
  group_by(site) %>% 
  summarize(num.samples = length(unique(ID)))  ## 5 samples each every site

#### Plankton Jars ####

## Load plankton jar data

plankton_jars <- read.csv("plankton_jars.csv", header = TRUE)

## Clean up plankton jar data

names(plankton_jars)
summary(plankton_jars$size.fraction)
summary(plankton_jars$shape)
summary(plankton_jars$colour)

summary(plankton_jars$shape)

plankton_jars$shape <- mapvalues(plankton_jars$shape,
                                 from = c('fiber',
                                          'fibre',
                                          'fragment'),
                                 to = c('Fibre',
                                        'Fibre',
                                        'Fragment'))
summary(plankton_jars$shape)

summary(plankton_jars$raman.ID)
plankton_jars$raman.ID <-
  mapvalues(plankton_jars$raman.ID,
            from = c('Cellulose\n',
                     'Plastic Dye'),
            to = c('Cellulose',
                   'Plastics Dye'))

plankton_jars$particle.type <- 
  mapvalues(plankton_jars$raman.ID,
            from = levels(plankton_jars$raman.ID),
            to = c('Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-Synthetic',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(plankton_jars$particle.type)

plankton_jars$particle.type <- as.character(plankton_jars$particle.type)

plankton_jars$particle.type[plankton_jars$colour == 'clear' &
                              plankton_jars$shape == 'Fibre' &
                              plankton_jars$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

plankton_jars$particle.type <- as.factor(plankton_jars$particle.type)

plankton_jars$num <- 
  with(plankton_jars,
       ifelse(is.na(length), 0, 1))  # Separate out samples that had 0 counts

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
  group_by(ID, site, sample.type, shape, colour, blank.match, particle.type, 
           raman.ID) %>% 
  summarize(count = sum(num))

## Blank subtract

PJ_polymer2 <-
  left_join(PJ_polymer, 
            PJ_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'raman.ID', 
                   'particle.type'))

PJ_polymer2$blank.mean[is.na(PJ_polymer2$blank.mean)] <- 0

PJ_polymer2$adj.count <- ceiling(with(PJ_polymer2, count - blank.mean))
PJ_polymer2$adj.count[PJ_polymer2$adj.count < 0] <- 0

PJ_polymer2 %>% 
  group_by(site) %>% 
  summarize(num.samples = length(unique(ID)))  ## 5 samples each every site

#### Mussels ####

## Load mussels data

mussels <- read.csv("mussels.csv", header = TRUE)

## Clean up mussels data

names(mussels)
summary(mussels$size.fraction)
summary(mussels$shape)
summary(mussels$colour)
summary(mussels$site)

mussels$shape <- mapvalues(mussels$shape,
                                 from = levels(mussels$shape),
                                 to = c('Fibre',
                                        'Fibre',
                                        'Film',
                                        'Fragment'))
summary(mussels$shape)

mussels$colour <- mapvalues(mussels$colour,
                            from = 'clear\n',
                            to = 'clear')

summary(mussels$raman.ID)

mussels$raman.ID <- mapvalues(mussels$raman.ID,
                              from = c('Acrylic\n',
                                       'Acyrlic',
                                       'Cellulosee'),
                              to = c('Acyrlic',
                                     'Acrylic',
                                     'Cellulose'))

mussels$particle.type <- 
  mapvalues(mussels$raman.ID,
            from = levels(mussels$raman.ID),
            to = c('Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic'))
summary(mussels$particle.type)

mussels$particle.type <- as.character(mussels$particle.type)

mussels$particle.type[mussels$colour == 'clear' &
                              mussels$shape == 'Fibre' &
                              mussels$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

mussels$particle.type <- as.factor(mussels$particle.type)

mussels$num <- with(mussels,
                          ifelse(is.na(length), 0, 1))

## Separate blanks data

MU_blanks <- subset(mussels, sample.type == 'Blanks')
MU <- subset(mussels, sample.type == 'Mussels')
MU$ID <- as.character(MU$ID)
MU$ID <- as.factor(MU$ID)

## Summarize blanks data

summary(MU_blanks)
MU_blanks_particle_type <- 
  MU_blanks %>% 
  group_by(ID, shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.count = sum(num))

MU_blanks_means <- 
  MU_blanks_particle_type %>% 
  group_by(shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.mean = mean(blank.count))

## Summarize MU data

MU_polymer <- 
  MU %>% 
  group_by(ID, site, sample.type, size.fraction, shape, colour, blank.match, 
           particle.type, raman.ID) %>% 
  summarize(count = sum(num))

## Blank subtract

MU_polymer2 <-
  left_join(MU_polymer, 
            MU_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'raman.ID', 
                   'particle.type'))

MU_polymer2$blank.mean[is.na(MU_polymer2$blank.mean)] <- 0

MU_polymer2$adj.count <- ceiling(with(MU_polymer2, count - blank.mean))
MU_polymer2$adj.count[MU_polymer2$adj.count < 0] <- 0

MU_polymer2 %>% 
  group_by(site) %>% 
  summarize(num.samples = 
              length(unique(ID)))  ## 19 samples CB, 15 EB, 15 VH

#### Clams ####

# Load clams data

clams <- read.csv("clams.csv", header = TRUE)

## Clean up clams data

names(clams)
summary(clams$size.fraction)
summary(clams$shape)
summary(clams$colour)

clams$shape <- mapvalues(clams$shape,
                           from = levels(clams$shape),
                           to = c('Fibre',
                                  'Fibre'))
summary(clams$shape)

summary(clams$raman.ID)

clams$raman.ID <- mapvalues(clams$raman.ID,
                            from = 'Cellulose\n',
                            to = 'Cellulose')

clams$particle.type <- 
  mapvalues(clams$raman.ID,
            from = levels(clams$raman.ID),
            to = c('Unknown',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymerr',
                   'Semi-synthetic',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(clams$particle.type)

clams$particle.type <- as.character(clams$particle.type)

clams$particle.type[clams$colour == 'clear' &
                        clams$shape == 'Fibre' &
                        clams$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

clams$particle.type <- as.factor(clams$particle.type)

clams$num <- with(clams,
                    ifelse(is.na(length), 0, 1))

## Summarize clams data

clams_polymer <- 
  clams %>% 
  group_by(ID, site, sample.type, size.fraction, shape, colour, blank.match, 
           particle.type, raman.ID) %>% 
  summarize(count = sum(num))

## Blank subtract

clams_polymer2 <-
  left_join(clams_polymer, 
            MU_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'raman.ID', 
                   'particle.type'))

clams_polymer2$blank.mean[is.na(clams_polymer2$blank.mean)] <- 0

clams_polymer2$adj.count <- ceiling(with(clams_polymer2, count - blank.mean))
clams_polymer2$adj.count[clams_polymer2$adj.count < 0] <- 0

clams_polymer2 %>% 
  group_by(site) %>% 
  summarize(num.samples = length(unique(ID)))  ## 13 samples CB, 15 EB

#### Sea Cucumbers #####

# Load sea_cucumbers data

sea_cucumbers <- read.csv("sea_cucumbers.csv", header = TRUE)

## Clean up sea cucumbers data

names(sea_cucumbers)
summary(sea_cucumbers$size.fraction)
summary(sea_cucumbers$shape)
summary(sea_cucumbers$colour)

sea_cucumbers$shape <- mapvalues(sea_cucumbers$shape,
                                 from = c('fibre',
                                          'fibre ',
                                          'fragment'),
                                 to = c('Fibre',
                                        'Fibre',
                                        'Fragment'))
summary(sea_cucumbers$shape)

levels(sea_cucumbers$raman.ID)

sea_cucumbers$particle.type <- 
  mapvalues(sea_cucumbers$raman.ID,
            from = levels(sea_cucumbers$raman.ID),
            to = c('Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Natural',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Unknown'))
summary(sea_cucumbers$particle.type)

sea_cucumbers$particle.type <- as.character(sea_cucumbers$particle.type)

sea_cucumbers$particle.type[sea_cucumbers$colour == 'clear' &
                      sea_cucumbers$shape == 'Fibre' &
                      sea_cucumbers$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

sea_cucumbers$particle.type <- as.factor(sea_cucumbers$particle.type)

sea_cucumbers$num <- with(sea_cucumbers,
                          ifelse(is.na(length), 0, 1))

## Separate blanks data

CU_blanks <- subset(sea_cucumbers, sample.type == 'Blanks')
CU <- subset(sea_cucumbers, sample.type == 'Sea Cucumbers')
CU$ID <- as.character(CU$ID)
CU$ID <- as.factor(CU$ID)

## Summarize blanks data

summary(CU_blanks)
CU_blanks_particle_type <- 
  CU_blanks %>% 
  group_by(ID, shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.count = sum(num))

CU_blanks_means <- 
  CU_blanks_particle_type %>% 
  group_by(shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.mean = mean(blank.count))

## Summarize CU data

CU_polymer <- 
  CU %>% 
  group_by(ID, site, sample.type, size.fraction, shape, colour, blank.match, 
           particle.type, raman.ID) %>% 
  summarize(count = sum(num))

## Blank subtract

CU_polymer2 <-
  left_join(CU_polymer, 
            CU_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'raman.ID', 
                   'particle.type'))

CU_polymer2$blank.mean[is.na(CU_polymer2$blank.mean)] <- 0

CU_polymer2$adj.count <- ceiling(with(CU_polymer2, count - blank.mean))
CU_polymer2$adj.count[CU_polymer2$adj.count < 0] <- 0

## Remove samples where all size fractions aren't counted yet

CU_incomplete <- 
  CU_polymer2 %>% 
  group_by(ID) %>%
  summarize(count = length(unique(size.fraction))) %>% 
  filter(count < 2)

CU_polymer2$incomplete <- CU_polymer2$ID %in% CU_incomplete$ID
  
CU_polymer3 <-
  CU_polymer2 %>% 
  filter(incomplete == 'FALSE')

CU_polymer3 %>% 
  group_by(site) %>% 
  summarize(num.samples = length(unique(ID)))  ## 16 samples CB, 15 EB, 15 VH

#### Sea Stars ####

# Load sea_stars data

sea_stars <- read.csv("sea_stars.csv", header = TRUE)

## Clean up sea_stars data

names(sea_stars)
summary(sea_stars$size.fraction)
summary(sea_stars$shape)
summary(sea_stars$colour)

sea_stars$shape <- mapvalues(sea_stars$shape,
                         from = c('fiber',
                                  'fibre',
                                  'fibre ',
                                  'fragment'),
                         to = c('Fibre',
                                'Fibre',
                                'Fibre',
                                'Fragment'))
sea_stars$shape[sea_stars$shape == ''] <- NA
summary(sea_stars$shape)

summary(sea_stars$raman.ID)

sea_stars$particle.type <- 
  mapvalues(sea_stars$raman.ID,
            from = levels(sea_stars$raman.ID),
            to = c('Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Natural',
                   'Synthetic Polymer',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(sea_stars$particle.type)

sea_stars$particle.type <- as.character(sea_stars$particle.type)

sea_stars$particle.type[sea_stars$colour == 'clear' &
                              sea_stars$shape == 'Fibre' &
                              sea_stars$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

sea_stars$particle.type <- as.factor(sea_stars$particle.type)

sea_stars$num <- with(sea_stars,
                  ifelse(is.na(length), 0, 1))

## Separate blanks data

SS_blanks <- subset(sea_stars, sample.type == 'Blanks')
SS <- subset(sea_stars, sample.type == 'Sea Stars')
SS$ID <- as.character(SS$ID)
SS$ID <- as.factor(SS$ID)

## Summarize blanks data

summary(SS_blanks)
SS_blanks_particle_type <- 
  SS_blanks %>% 
  group_by(ID, shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.count = sum(num))

SS_blanks_means <- 
  SS_blanks_particle_type %>% 
  group_by(shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.mean = mean(blank.count))

## Summarize SS data

SS_polymer <- 
  SS %>% 
  group_by(ID, site, sample.type, size.fraction, shape, colour, blank.match,
           particle.type, raman.ID) %>% 
  summarize(count = sum(num))

## Blank subtract

SS_polymer2 <-
  left_join(SS_polymer, 
            SS_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'raman.ID', 
                   'particle.type'))

SS_polymer2$blank.mean[is.na(SS_polymer2$blank.mean)] <- 0

SS_polymer2$adj.count <- ceiling(with(SS_polymer2, count - blank.mean))
SS_polymer2$adj.count[SS_polymer2$adj.count < 0] <- 0

## Add together samples that were divided

SS_polymer2$ID <- mapvalues(SS_polymer2$ID,
                            from = c('CBSS3 1/2',
                                     'CBSS3 2/2',
                                     'CBSS8 (1/4)', 
                                     'CBSS8 (2/4)',
                                     'CBSS8 (3/4)',
                                     'CBSS8 (4/4)',
                                     'EBSS10-1',
                                     'EBSS10-2',
                                     'EBSS10-3',
                                     'EBSS11 (1/3)',
                                     'EBSS11 (2/3)',
                                     'EBSS11 (3/3)'),
                            to = c('CBSS3', 'CBSS3',
                                   'CBSS8', 'CBSS8', 'CBSS8', 'CBSS8',
                                   'EBSS10', 'EBSS10', 'EBSS10',
                                   'EBSS11', 'EBSS11', 'EBSS11'))

## Remove samples where all size fractions aren't counted yet


SS_incomplete <- 
  SS_polymer2 %>% 
  group_by(ID) %>%
  summarize(count = length(unique(size.fraction))) %>% 
  filter(count < 2)

SS_polymer2$incomplete <- SS_polymer2$ID %in% SS_incomplete$ID

SS_polymer3 <-
  SS_polymer2 %>% 
  filter(incomplete == 'FALSE')

SS_polymer3 %>% 
  group_by(site) %>% 
  summarize(num.samples = length(unique(ID)))  ## 15 samples CB, 15 EB


#### Crabs ####

# Load crabs data

crabs <- read.csv("crabs.csv", header = TRUE)

## Clean up crabs data

names(crabs)
summary(crabs$size.fraction)
summary(crabs$shape)
summary(crabs$colour)
crabs$colour <- mapvalues(crabs$colour,
                          from = c('Black',
                                   'Blue',
                                   'Brown',
                                   'Clear',
                                   'Purple'),
                          to = c('black',
                                 'blue',
                                 'brown',
                                 'clear',
                                 'purple'))

crabs$shape <- mapvalues(crabs$shape,
                                 from = c('fibre',
                                          'fibre ',
                                          'fragment',
                                          'fiber'),
                                 to = c('Fibre',
                                        'Fibre',
                                        'Fragment',
                                        'Fibre'))
summary(crabs$shape)

levels(crabs$raman.ID)

crabs$particle.type <- 
  mapvalues(crabs$raman.ID,
            from = levels(crabs$raman.ID),
            to = c('Unknown',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(crabs$particle.type)

crabs$particle.type <- as.character(crabs$particle.type)

crabs$particle.type[crabs$colour == 'clear' &
                              crabs$shape == 'Fibre' &
                              crabs$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

crabs$particle.type <- as.factor(crabs$particle.type)

crabs$num <- with(crabs, ifelse(is.na(length), 0, 1))

## Separate blanks data

CR_blanks <- subset(crabs, sample.type == 'Blanks')
CR <- subset(crabs, sample.type == 'Crabs')
CR$ID <- as.character(CR$ID)
CR$ID <- as.factor(CR$ID)

## Summarize blanks data

summary(CR_blanks)
CR_blanks_particle_type <- 
  CR_blanks %>% 
  group_by(ID, shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.count = sum(num))

CR_blanks_means <- 
  CR_blanks_particle_type %>% 
  group_by(shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.mean = mean(blank.count))

## Summarize CR data

CR_polymer <- 
  CR %>% 
  group_by(ID, site, sample.type, size.fraction, shape, colour, blank.match, 
           particle.type, 
           raman.ID) %>% 
  summarize(count = sum(num))

## Blank subtract

CR_polymer2 <-
  left_join(CR_polymer, 
            CR_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'raman.ID', 
                   'particle.type'))

CR_polymer2$blank.mean[is.na(CR_polymer2$blank.mean)] <- 0

CR_polymer2$adj.count <- ceiling(with(CR_polymer2, count - blank.mean))
CR_polymer2$adj.count[CR_polymer2$adj.count < 0] <- 0

## Remove samples where all size fractions aren't counted yet

CR_incomplete <- 
  CR_polymer2 %>% 
  group_by(ID) %>%
  summarize(count = length(unique(size.fraction))) %>% 
  filter(count < 2)

CR_polymer2$incomplete <- CR_polymer2$ID %in% CR_incomplete$ID

CR_polymer3 <-
  CR_polymer2 %>% 
  filter(incomplete == 'FALSE')

CR_polymer3 %>% 
  group_by(site) %>% 
  summarize(num.samples = length(unique(ID)))  ## 15 samples CB, 15 EB, 13 VH

#### Surfperch ####

# Load surfperch data

surfperch <- read.csv("surfperch.csv", header = TRUE)

## Clean up surfperch data

names(surfperch)
summary(surfperch$size.fraction)
surfperch$size.fraction <- mapvalues(surfperch$size.fraction,
                                     from = '>1',
                                     to = '1-150')

summary(surfperch$shape)
summary(surfperch$colour)

surfperch$shape <- mapvalues(surfperch$shape,
                         from = levels(surfperch$shape),
                         to = c('Fibre',
                                'Fibre',
                                'Film',
                                'Fragment',
                                'Sphere'))
summary(surfperch$shape)

levels(surfperch$raman.ID)

surfperch$raman.ID <-
  mapvalues(surfperch$raman.ID,
            from = c('Cellulose\n',
                     'Nylon\n',
                     'Polyester\n',
                     'Polyster',
                     'Unknown\n'),
            to = c('Cellulose',
                   'Nylon',
                   'Polyester',
                   'Polyester',
                   'Unknown'))

surfperch$particle.type <- 
  mapvalues(surfperch$raman.ID,
            from = levels(surfperch$raman.ID),
            to = c('Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Natural',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(surfperch$particle.type)

surfperch$particle.type <- as.character(surfperch$particle.type)

surfperch$particle.type[surfperch$colour == 'clear' &
                      surfperch$shape == 'Fibre' &
                      surfperch$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

surfperch$particle.type <- as.factor(surfperch$particle.type)

surfperch$num <- with(surfperch,
                  ifelse(is.na(length), 0, 1))

## Separate blanks data

SP_blanks <- subset(surfperch, sample.type == 'Blanks')
SP <- subset(surfperch, 
             sample.type == 'Surfperch Guts' |
             sample.type == 'Surfperch Livers'  , )
SP$ID <- as.character(SP$ID)
SP$ID <- as.factor(SP$ID)

## Summarize blanks data

summary(SP_blanks)
SP_blanks_particle_type <- 
  SP_blanks %>% 
  group_by(ID, shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.count = sum(num))

SP_blanks_means <- 
  SP_blanks_particle_type %>% 
  group_by(shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.mean = mean(blank.count))

## Summarize SP data

SP_polymer <- 
  SP %>% 
  group_by(ID, site, sample.type, size.fraction, shape, colour, blank.match, 
           particle.type, raman.ID) %>% 
  summarize(count = sum(num))

## Blank subtract

SP_polymer2 <-
  left_join(SP_polymer, 
            SP_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'raman.ID', 
                   'particle.type'))

SP_polymer2$blank.mean[is.na(SP_polymer2$blank.mean)] <- 0

SP_polymer2$adj.count <- ceiling(with(SP_polymer2, count - blank.mean))
SP_polymer2$adj.count[SP_polymer2$adj.count < 0] <- 0

## Remove samples where all size fractions aren't counted yet

SP_incomplete <- 
  SP_polymer2 %>% 
  group_by(ID) %>%
  summarize(count = length(unique(size.fraction)),
            sample.type = length(unique(sample.type))) %>% 
  filter(count < 2 | sample.type < 2)

SP_polymer2$incomplete <- SP_polymer2$ID %in% SP_incomplete$ID

SP_polymer3 <-
  SP_polymer2 %>% 
  filter(incomplete == 'FALSE')

SP_polymer3 %>% 
  group_by(site) %>% 
  summarize(num.samples = length(unique(ID)))  ## 11 samples CB, 10 EB, 12 VH

#### Flatfish ####

# Load surfperch data

flatfish <- read.csv("flatfish.csv", header = TRUE)

## Clean up flatfish data

names(flatfish)
summary(flatfish$size.fraction)
flatfish$size.fraction <- mapvalues(flatfish$size.fraction,
                                     from = '>1',
                                     to = '1-150')

summary(flatfish$shape)
summary(flatfish$colour)

flatfish$shape <- mapvalues(flatfish$shape,
                             from = levels(flatfish$shape),
                             to = c('Fibre',
                                    'Fibre',
                                    'Fragment'))
summary(flatfish$shape)

levels(flatfish$raman.ID)

flatfish$raman.ID <-
  mapvalues(flatfish$raman.ID,
            from = c('Paint Chip',
                     'Wool\n'),
            to = c('Paint',
                   'Wool'))

flatfish$particle.type <- 
  mapvalues(flatfish$raman.ID,
            from = levels(flatfish$raman.ID),
            to = c('Unknown',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(flatfish$particle.type)

flatfish$particle.type <- as.character(flatfish$particle.type)

flatfish$particle.type[flatfish$colour == 'clear' &
                          flatfish$shape == 'Fibre' &
                          flatfish$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

flatfish$particle.type <- as.factor(flatfish$particle.type)

flatfish$num <- with(flatfish,
                      ifelse(is.na(length), 0, 1))

## Separate blanks data

FF_blanks <- subset(flatfish, sample.type == 'Blanks')
FF <- subset(flatfish, 
             sample.type == 'Flatfish Guts' |
               sample.type == 'Flatfish Livers'  , )
FF$ID <- as.character(FF$ID)
FF$ID <- as.factor(FF$ID)

## Summarize blanks data

summary(FF_blanks)
FF_blanks_particle_type <- 
  FF_blanks %>% 
  group_by(ID, shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.count = sum(num))

FF_blanks_means <- 
  FF_blanks_particle_type %>% 
  group_by(shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.mean = mean(blank.count))

## Summarize FF data

FF_polymer <- 
  FF %>% 
  group_by(ID, site, sample.type, size.fraction, shape, colour, blank.match, 
           particle.type, raman.ID) %>% 
  summarize(count = sum(num))

## Blank subtract

FF_polymer2 <-
  left_join(FF_polymer, 
            FF_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'raman.ID', 
                   'particle.type'))

FF_polymer2$blank.mean[is.na(FF_polymer2$blank.mean)] <- 0

FF_polymer2$adj.count <- ceiling(with(FF_polymer2, count - blank.mean))
FF_polymer2$adj.count[FF_polymer2$adj.count < 0] <- 0

## Remove samples where all size fractions aren't counted yet

FF_incomplete <- 
  FF_polymer2 %>% 
  group_by(ID) %>%
  summarize(count = length(unique(size.fraction)),
            sample.type = length(unique(sample.type))) %>% 
  filter(count < 2 | sample.type < 2)

FF_polymer2$incomplete <- FF_polymer2$ID %in% FF_incomplete$ID

FF_polymer3 <-
  FF_polymer2 %>% 
  filter(incomplete == 'FALSE')

FF_polymer3 %>% 
  group_by(site) %>% 
  summarize(num.samples = length(unique(ID)))  ## 10 samples CB


#### Rockfish ####

# Load rockfish data

rockfish <- read.csv("rockfish.csv", header = TRUE)

## Clean up crabs data

names(rockfish)
summary(rockfish$size.fraction)
rockfish$size.fraction <- mapvalues(rockfish$size.fraction,
                                    from = '1-150',
                                    to = '>1')

summary(rockfish$shape)
summary(rockfish$colour)

rockfish$shape <- mapvalues(rockfish$shape,
                            from = levels(rockfish$shape),
                            to = c('Fibre',
                                   'Fragment'))
summary(rockfish$shape)

levels(rockfish$raman.ID)

rockfish$raman.ID <-
  mapvalues(rockfish$raman.ID,
            from = c('Cellulose\n'),
            to = c('Cellulose'))

rockfish$particle.type <- 
  mapvalues(rockfish$raman.ID,
            from = levels(rockfish$raman.ID),
            to = c('Unknown',
                   'Synthetic Polymer',
                   'Natural',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(rockfish$particle.type)

summary(rockfish$colour)

rockfish$particle.type <- as.character(rockfish$particle.type)

rockfish$particle.type[rockfish$colour == 'clear' &
                         rockfish$shape == 'Fibre' &
                         rockfish$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

rockfish$particle.type <- as.factor(rockfish$particle.type)

rockfish$num <- with(rockfish,
                     ifelse(is.na(length), 0, 1))

## Separate blanks data

RF_blanks <- subset(rockfish, sample.type == 'Blanks')
RF <- subset(rockfish, 
             sample.type == 'Rockfish Guts' |
               sample.type == 'Rockfish Livers' |
               sample.type == 'Rockfish: Ingested Animals')
RF$ID <- as.character(RF$ID)
RF$ID <- as.factor(RF$ID)

## Summarize blanks data

summary(RF_blanks)
RF_blanks_particle_type <- 
  RF_blanks %>% 
  group_by(ID, shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.count = sum(num))

RF_blanks_means <- 
  RF_blanks_particle_type %>% 
  group_by(shape, colour, blank.match, raman.ID, particle.type) %>% 
  summarize(blank.mean = mean(blank.count))

## Summarize RF data

RF_polymer <- 
  RF %>% 
  group_by(ID, site, sample.type, size.fraction, shape, colour, blank.match, 
           particle.type, raman.ID) %>% 
  summarize(count = sum(num))

## Blank subtract

RF_polymer2 <-
  left_join(RF_polymer, 
            RF_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'raman.ID', 
                   'particle.type'))

RF_polymer2$blank.mean[is.na(RF_polymer2$blank.mean)] <- 0

RF_polymer2$adj.count <- ceiling(with(RF_polymer2, count - blank.mean))
RF_polymer2$adj.count[RF_polymer2$adj.count < 0] <- 0

## Add together samples that were divided

RF_polymer2$ID <- mapvalues(RF_polymer2$ID,
                            from = c('EBRF16A', 
                                     'EBRF16B',
                                     'HPRF21A',
                                     'HPRF21B',
                                     'HPRF4B'),
                            to = c('EBRF16',
                                   'EBRF16',
                                   'HPRF21',
                                   'HPRF21',
                                   'HPRF4'))

## Remove samples where all size fractions aren't counted yet

RF_incomplete <- 
  RF_polymer2 %>% 
  group_by(ID) %>%
  summarize(count = length(unique(size.fraction)),
            sample.type = length(unique(sample.type))) %>% 
  filter(count < 3 | sample.type < 2)

RF_polymer2$incomplete <- RF_polymer2$ID %in% RF_incomplete$ID

RF_polymer3 <-
  RF_polymer2 %>% 
  filter(incomplete == 'FALSE')

RF_polymer2 %>% 
  group_by(site) %>% 
  summarize(num.samples = 
              length(unique(ID)))  ## 10 samples CB, 21 EB, 22 VH
