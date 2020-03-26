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
plankton_tows$particle.type[is.na(plankton_tows$particle.type)] <- 'Unknown'

summary(plankton_tows$colour)

## Clean up plankton tow field data 

summary(PT_field$sample.volume)

PT <- left_join(plankton_tows, PT_field, by = 'ID')

PT$ID <- as.factor(PT$ID)

head(PT)

PT$num <- ifelse(is.na(PT$length), 0, 1)  # Separate out samples that had 0 counts

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
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Natural Anthropogenic',
                   'Unknown Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-Synthetic',
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

#### Mussels ####

## Load mussels data

mussels <- read.csv("mussels.csv", header = TRUE)

## Clean up mussels data

names(mussels)
summary(mussels$size.fraction)
summary(mussels$shape)
summary(mussels$colour)

mussels$shape <- mapvalues(mussels$shape,
                                 from = c('fibre',
                                          'fragment'),
                                 to = c('Fibre',
                                        'Fragment'))
summary(mussels$shape)

summary(mussels$raman.ID)

mussels$particle.type <- 
  mapvalues(mussels$raman.ID,
            from = levels(mussels$raman.ID),
            to = c('Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Natural Anthropogenic',
                   'Natural Anthropogenic',
                   'Unknown Anthropogenic',
                   'Unknown Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic'))
summary(mussels$particle.type)

summary(mussels$colour)

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

#### Clams ####

# Load clams data

clams <- read.csv("clams.csv", header = TRUE)

## Clean up clams data

names(clams)
summary(clams$size.fraction)
summary(clams$shape)
summary(clams$colour)

clams$shape <- mapvalues(clams$shape,
                           from = c('fibre',
                                    'fragment'),
                           to = c('Fibre',
                                  'Fragment'))
summary(clams$shape)

summary(clams$raman.ID)

clams$particle.type <- 
  mapvalues(clams$raman.ID,
            from = levels(clams$raman.ID),
            to = c('Unknown',
                   'Natural Anthropogenic',
                   'Natural Anthropogenic',
                   'Unknown Anthropogenic',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(clams$particle.type)

summary(clams$colour)

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
                   'Unknown Anthropogenic',
                   'Synthetic Polymer',
                   'Natural',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Unknown'))
summary(sea_cucumbers$particle.type)

summary(sea_cucumbers$colour)

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

#### Sea Stars ####

# Load sea_stars data

sea_stars <- read.csv("sea_stars.csv", header = TRUE)

## Clean up sea_stars data

names(sea_stars)
summary(sea_stars$size.fraction)
summary(sea_stars$shape)
summary(sea_stars$colour)

sea_stars$shape <- mapvalues(sea_stars$shape,
                         from = c('fibre',
                                  'fibre ',
                                  'fragment'),
                         to = c('Fibre',
                                'Fibre',
                                'Fragment'))
summary(sea_stars$shape)

summary(sea_stars$raman.ID)

sea_stars$particle.type <- 
  mapvalues(sea_stars$raman.ID,
            from = levels(sea_stars$raman.ID),
            to = c('Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Unknown Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Natural',
                   'Synthetic Polymer',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(sea_stars$particle.type)

summary(sea_stars$colour)

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
                            from = c('CBSS8 (1/4)', 
                                     'CBSS8 (2/4)',
                                     'CBSS8 (3/4)',
                                     'CBSS8 (4/4)',
                                     'EBSS11 (1/3)',
                                     'EBSS11 (2/3)',
                                     'EBSS11 (3/3)'),
                            to = c('CBSS8', 'CBSS8', 'CBSS8', 'CBSS8',
                                   'EBSS11', 'EBSS11', 'EBSS11'))

#### Crabs ####

# Load crabs data

crabs <- read.csv("crabs.csv", header = TRUE)

## Clean up crabs data

names(crabs)
summary(crabs$size.fraction)
summary(crabs$shape)
summary(crabs$colour)

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
                   'Unknown Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(crabs$particle.type)

summary(crabs$colour)

crabs$num <- with(crabs,
                          ifelse(is.na(length), 0, 1))

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

#### Surfperch ####

# Load surfperch data

surfperch <- read.csv("surfperch.csv", header = TRUE)

## Clean up crabs data

names(surfperch)
summary(surfperch$size.fraction)
summary(surfperch$shape)
summary(surfperch$colour)

surfperch$shape <- mapvalues(surfperch$shape,
                         from = levels(surfperch$shape),
                         to = c('Fibre',
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
                     'Unknown\n'),
            to = c('Cellulose',
                   'Nylon',
                   'Polyester',
                   'Unknown'))

surfperch$particle.type <- 
  mapvalues(surfperch$raman.ID,
            from = levels(surfperch$raman.ID),
            to = c('Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Unknown Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Natural',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(surfperch$particle.type)

summary(surfperch$colour)

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

#### Flatfish ####

# Load flatfish data

flatfish <- read.csv("flatfish.csv", header = TRUE)

## Clean up crabs data

names(flatfish)
summary(flatfish$size.fraction)
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
            from = c('Wool\n'),
            to = c('Wool'))

flatfish$particle.type <- 
  mapvalues(flatfish$raman.ID,
            from = levels(flatfish$raman.ID),
            to = c('Unknown',
                   'Natural Anthropogenic',
                   'Unknown Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Unknown',
                   'Natural Anthropogenic'))
summary(flatfish$particle.type)

summary(flatfish$colour)

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

#### Rockfish ####