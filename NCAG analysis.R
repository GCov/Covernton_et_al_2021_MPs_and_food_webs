## load packages
library(plyr)
library(ggplot2)
library(dplyr)
library(ggridges)

## define standard error function

se <- function(x) {
  sd(x)/sqrt(length(x))
}

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

PT_particle_type <-
  PT_polymer2 %>%
  group_by(ID, site, sample.type, particle.type, sample.volume) %>%
  summarize(count = sum(adj.count))

PT_summary <-
  PT_particle_type %>% 
  group_by(ID, site, sample.type, sample.volume) %>% 
  summarize(count = sum(count))
  
ggplot(PT_summary) + 
  geom_boxplot(aes(x = site,
               y = count/sample.volume),
           colour = 'black',
           size = 1) + 
  labs(x = 'Particles/L',
       y = 'Site')

with(PT_summary,tapply(count/sample.volume, site, mean))
with(PT_summary,tapply(count/sample.volume, site, se))

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

PJ_particle_type <-
  PJ_polymer2 %>%
  group_by(ID, site, sample.type, particle.type) %>%
  summarize(count = sum(adj.count))

ID <- c('CBPJ1', 'CBPJ2', 'CBPJ3', 'CBPJ4', 'CBPJ5', 
        'EBPJ1', 'EBPJ2', 'EBPJ3', 'EBPJ4', 'EBPJ5', 
        'HPPJ1', 'HPPJ2', 'HPPJ3', 'HPPJ4', 'HPPJ5')
site <- as.factor(c(
  rep("Cole's Bay", 5),
  rep('Elliot Bay', 5),
  rep('Victoria Harbour', 5)
))
PJ_IDs <- data.frame(ID, site)

PJ_particle_type <- left_join(PJ_IDs, PJ_particle_type, by = c('ID', 'site'))

PJ_particle_type$count[is.na(PJ_synthetic$count)] <- 0

PJ_particle_type$sample.type[is.na(PJ_particle_type$sample.type)] <- 
  'Plankton Jars'

PJ_summary <-
  PJ_particle_type %>% 
  group_by(ID, site, sample.type) %>% 
  summarize(count = sum(count))

PJ_summary$count[is.na(PJ_summary$count)] <- 0
            
ggplot(PJ_summary) + 
  geom_boxplot(aes(x = site,
                   y = count),
               colour = 'black',
               size = 1) + 
  labs(x = 'Site',
       y = 'Particles/L')

with(PJ_summary, tapply(count, site, mean))
with(PJ_summary, tapply(count, site, se))

## Do all of the above for the mussels data

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
  group_by(ID, site, sample.type, shape, colour, blank.match, particle.type, 
           raman.ID) %>% 
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

MU_particle_type <-
  MU_polymer2 %>%
  group_by(ID, site, sample.type, particle.type) %>%
  summarize(count = sum(adj.count))

MU_synthetic <- subset(MU_particle_type,
                       particle.type == 'Synthetic Polymer')

ID <- c(levels(MU$ID))

site <- as.factor(c(
  rep("Cole's Bay", 16),
  rep('Elliot Bay', 10),
  rep('Victoria Harbour', 10)
))
MU_IDs <- data.frame(ID, site)

MU_synthetic <- left_join(MU_IDs, MU_synthetic, by = c('ID', 'site'))

MU_synthetic$count[is.na(MU_synthetic$count)] <- 0

MU_synthetic$sample.type[is.na(MU_synthetic$sample.type)] <- 'Mussels'

MU_synthetic_summary <-
  MU_synthetic %>% 
  group_by(site) %>% 
  summarize(mean = mean(count),
            sd = sd(count))

ggplot(MU_synthetic) + 
  geom_boxplot(aes(x = site,
                   y = count),
               colour = 'black',
               size = 1) +
  labs(x = 'Site',
       y = 'Particles/Ind')


## Do all of the above for clams

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
  group_by(ID, site, sample.type, shape, colour, blank.match, particle.type, 
           raman.ID) %>% 
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

clams_particle_type <-
  clams_polymer2 %>%
  group_by(ID, site, sample.type, particle.type) %>%
  summarize(count = sum(adj.count))

clams_synthetic <- subset(clams_particle_type,
                       particle.type == 'Synthetic Polymer')

ID <- c(levels(clams$ID))

site <- as.factor(c(
  rep("Cole's Bay", 8),
  rep('Elliot Bay', 10)))
  
clams_IDs <- data.frame(ID, site)

clams_synthetic <- left_join(clams_IDs, clams_synthetic, by = c('ID', 'site'))

clams_synthetic$count[is.na(clams_synthetic$count)] <- 0

clams_synthetic$sample.type[is.na(clams_synthetic$sample.type)] <- 'Clams'

clams_synthetic_summary <-
  clams_synthetic %>% 
  group_by(site) %>% 
  summarize(mean = mean(count),
            sd = sd(count))

ggplot(clams_synthetic) + 
  geom_boxplot(aes(x = site,
                   y = count),
               colour = 'black',
               size = 1) +
  labs(x = 'Site',
       y = 'Particles/Ind')


## Do all of the above for sea stars

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
  group_by(ID, site, sample.type, shape, colour, blank.match, particle.type, 
           raman.ID) %>% 
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

SS_particle_type <-
  SS_polymer2 %>%
  group_by(ID, site, sample.type, particle.type) %>%
  summarize(count = sum(adj.count))

SS_synthetic <- subset(SS_particle_type,
                       particle.type == 'Synthetic Polymer')

SS2 <- subset(SS, size.fraction == '1-150')
SS2$ID <- as.character(SS2$ID)
SS2$ID <- as.factor(SS2$ID)

ID <- levels(SS2$ID)

ID <- mapvalues(ID,
                from = c('CBSS8 (1/4)', 
                         'CBSS8 (2/4)',
                         'CBSS8 (3/4)',
                         'CBSS8 (4/4)',
                         'EBSS11 (1/3)',
                         'EBSS11 (2/3)',
                         'EBSS11 (3/3)'),
                to = c('CBSS8', 'CBSS8', 'CBSS8', 'CBSS8',
                       'EBSS11', 'EBSS11', 'EBSS11'))

ID <- unique(ID)

site <- as.factor(c(
  rep("Cole's Bay", 11),
  rep('Elliot Bay', 10)))

SS_IDs <- data.frame(ID, site)

SS_synthetic <- left_join(SS_IDs, SS_synthetic, by = c('ID', 'site'))

SS_synthetic$count[is.na(SS_synthetic$count)] <- 0

SS_synthetic$sample.type[is.na(SS_synthetic$sample.type)] <- 'Sea Stars'

SS_synthetic_summary <-
  SS_synthetic %>% 
  group_by(site) %>% 
  summarize(mean = mean(count),
            sd = sd(count))

ggplot(SS_synthetic) + 
  geom_boxplot(aes(x = site,
                   y = count),
               colour = 'black',
               size = 1) +
  labs(x = 'Site',
       y = 'Particles/Ind')

## Do all of the above for sea cucumbers

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
  group_by(ID, site, sample.type, shape, colour, blank.match, particle.type, 
           raman.ID) %>% 
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

CU_particle_type <-
  CU_polymer2 %>%
  group_by(ID, site, sample.type, particle.type) %>%
  summarize(count = sum(adj.count))

CU_synthetic <- subset(CU_particle_type,
                       particle.type == 'Synthetic Polymer')

CU2 <- subset(CU, size.fraction == '1-150')
CU2$ID <- as.character(CU2$ID)
CU2$ID <- as.factor(CU2$ID)

ID <- levels(CU2$ID)

site <- as.factor(c(
  rep("Cole's Bay", 13),
  rep('Elliot Bay', 11),
  rep('Victoria Harbour', 11)))

CU_IDs <- data.frame(ID, site)

CU_synthetic <- left_join(CU_IDs, CU_synthetic, by = c('ID', 'site'))

CU_synthetic$count[is.na(CU_synthetic$count)] <- 0

CU_synthetic$sample.type[is.na(CU_synthetic$sample.type)] <- 'Sea Cucumbers'

CU_synthetic_summary <-
  CU_synthetic %>% 
  group_by(site) %>% 
  summarize(mean = mean(count),
            sd = sd(count))

ggplot(CU_synthetic) + 
  geom_boxplot(aes(x = site,
                   y = count),
               colour = 'black',
               size = 1) +
  labs(x = 'Site',
       y = 'Particles/Ind')

## do all of the above for crabs

# Load crabs data

crabs <- read.csv("crabs.csv", header = TRUE)

## Clean up sea cucumbers data

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
  group_by(ID, site, sample.type, shape, colour, blank.match, particle.type, 
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

CR_particle_type <-
  CR_polymer2 %>%
  group_by(ID, site, sample.type, particle.type) %>%
  summarize(count = sum(adj.count))

CR_synthetic <- subset(CR_particle_type,
                       particle.type == 'Synthetic Polymer')

CR2 <- subset(CR, size.fraction == '<1000')
CR2$ID <- as.character(CR2$ID)
CR2$ID <- as.factor(CR2$ID)

ID <- levels(CR2$ID)

site <- as.factor(c(
  rep("Cole's Bay", 11),
  rep('Elliot Bay', 11),
  rep('Victoria Harbour', 10)))

CR_IDs <- data.frame(ID, site)

CR_synthetic <- left_join(CR_IDs, CR_synthetic, by = c('ID', 'site'))

CR_synthetic$count[is.na(CR_synthetic$count)] <- 0

CR_synthetic$sample.type[is.na(CR_synthetic$sample.type)] <- 'Crabs'

CR_synthetic_summary <-
  CR_synthetic %>% 
  group_by(site) %>% 
  summarize(mean = mean(count),
            sd = sd(count))

ggplot(CR_synthetic) + 
  geom_boxplot(aes(x = site,
                   y = count),
               colour = 'black',
               size = 1) +
  labs(x = 'Site',
       y = 'Particles/Ind')

## Plot everything together

PT_synthetic2 <- PT_synthetic
PT_synthetic2$count <- with(PT_synthetic2, count/sample.volume)

allcounts <- rbind.fill(PT_synthetic2[c(1:4, 6)],
                        PJ_synthetic,
                        MU_synthetic,
                        clams_synthetic,
                        SS_synthetic,
                        CU_synthetic,
                        CR_synthetic)

palette1 <- c('#A5A57D', '#E80D6D', '#0E6B1D', '#AC3E7E', '#91721E', '#84D6D8', 
              '#8F2C40')

png('NCAG Preliminary Counts.png', width = 20, height = 15, units = 'cm',
    pointsize = 16, res = 300)

ggplot(allcounts) +
  geom_violin(aes(x = site, 
                   y = count,
                   fill = sample.type)) +
  labs(x = '',
       y = 'MP Count per ind or L') +
  scale_fill_manual(values = palette1) +
  theme_bw() +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    text = element_text(size = 16),
    panel.spacing = unit(0.5, "lines"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()

## Ridgeline plot

ggplot(allcounts) +
  geom_density_ridges2(aes(x = count, 
                          y = sample.type,
                          fill = site),
                      alpha = 0.4,
                      panel_scaling = FALSE,
                      scale = 1.1,
                      rel_min_height = 0.01) +
  labs(x = 'MP Count per ind or L',
       y = '') +
  theme_ridges() +
  scale_fill_manual(values = palette1[5:7]) +
  scale_x_continuous(limits = c(0, 5)) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    text = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 16)
  )


## Summarize by polymer type

PT_polymer2 %>% group_by(sample.type, raman.ID) %>% summarize(count = sum(count))
PJ_polymer2 %>% group_by(sample.type, raman.ID) %>% summarize(count = sum(count))
MU_polymer2 %>% group_by(sample.type, raman.ID) %>% summarize(count = sum(count))
clams_polymer2 %>% group_by(sample.type, raman.ID) %>% summarize(count = sum(count))
SS_polymer2 %>% group_by(sample.type, raman.ID) %>% summarize(count = sum(count))
CU_polymer2 %>% group_by(sample.type, raman.ID) %>% summarize(count = sum(count))
