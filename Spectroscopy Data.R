##### Setup #####

## load packages
library(plyr)
library(ggplot2)
library(dplyr)
library(nnet)
library(MuMIn)
library(boot)

#### Load spectroscopy data ####

## Load plankton tow data

plankton_tows <- read.csv("plankton_tows.csv", header = TRUE)
PT_field <- read.csv("PT_field_data.csv", header = TRUE)

## Clean up and combine plankton tow field data 

summary(PT_field$sample.volume)

PT <- left_join(plankton_tows, PT_field, by = 'ID')

PT$ID <- as.factor(PT$ID)

head(PT)

## Clean up plankton tow data

names(PT)
summary(PT$size.fraction)
summary(PT$shape)
summary(PT$colour)

summary(PT$shape)

PT$shape <- mapvalues(PT$shape,
                                 from = c('fibre',
                                          'fibre ',
                                          'fragment'),
                                 to = c('Fibre',
                                        'Fibre',
                                        'Fragment'))
summary(PT$shape)

summary(PT$raman.ID)

PT$particle.type <- 
  mapvalues(PT$raman.ID,
            from = levels(PT$raman.ID),
            to = c(NA,
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

PT$particle.type[is.na(PT$particle.type)] <-
  with(subset(PT, is.na(particle.type)), 
       ifelse(is.na(length),
              NA,
              'Unknown'))

summary(PT$particle.type)

PT$particle.type[PT$shape == 'Fibre' &
                              PT$colour == 'clear' &
                              PT$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

PT$num <- 
  ifelse(is.na(PT$length), 0, 1)  # Separate out samples that had 0 counts

## Load plankton jars data

PJ <- read.csv("plankton_jars.csv", header = TRUE)  # plankton jars

## Clean up plankton jars data

names(PJ)
summary(PJ$size.fraction)
summary(PJ$shape)
summary(PJ$colour)

summary(PJ$shape)

PJ$shape <- mapvalues(PJ$shape,
                                 from = c('fiber',
                                          'fibre',
                                          'fragment'),
                                 to = c('Fibre',
                                        'Fibre',
                                        'Fragment'))
summary(PJ$shape)

summary(PJ$raman.ID)

PJ$raman.ID <-
  mapvalues(PJ$raman.ID,
            from = c('Cellulose\n',
                     'Plastic Dye'),
            to = c('Cellulose',
                   'Plastics Dye'))

PJ$particle.type <- 
  mapvalues(PJ$raman.ID,
            from = levels(PJ$raman.ID),
            to = c(NA,
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Unknown',
                   'Natural Anthropogenic'))

PJ$particle.type[is.na(PJ$particle.type)] <-
  with(subset(PJ, is.na(particle.type)), 
       ifelse(length == 'NA',
              NA,
              'Unknown'))

summary(PJ$particle.type)

PJ$particle.type <- as.character(PJ$particle.type)

PJ$particle.type[PJ$shape == 'Fibre' &
                              PJ$colour == 'clear' &
                              PJ$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

PJ$particle.type <- as.factor(PJ$particle.type)

summary(PJ$particle.type)

PJ$num <- 
  ifelse(is.na(PJ$length), 0, 1)  # Separate out samples that had 0 counts



## Load all other data

MU <- read.csv("mussels.csv", header = TRUE)  # mussels
C <- read.csv("clams.csv", header = TRUE)  # clams
CU <- read.csv("sea_cucumbers.csv", header = TRUE)  # sea cucumbers
CR <- read.csv("crabs.csv", header = TRUE)  # crabs
SS <- read.csv("sea_stars.csv", header = TRUE)  # sea stars
SP <- read.csv("surfperch.csv", header = TRUE)  # surfperch
FF <- read.csv("flatfish.csv", header = TRUE)  # flatfish
RF <- read.csv("rockfish.csv", header = TRUE)  # rockfish

#### Combine all data and clean up ####

## Combine all animal particle counts data into one dataframe

full_spec_data <-
  rbind.data.frame(MU,
                   C,
                   CU,
                   CR,
                   SS,
                   SP,
                   FF,
                   RF)

## Clean up

names(full_spec_data)
full_spec_data$blank.match <- as.factor(full_spec_data$blank.match)
summary(full_spec_data$blank.match)
summary(full_spec_data$size.fraction)
summary(full_spec_data$shape)
summary(full_spec_data$colour)

full_spec_data$colour <-
  mapvalues(full_spec_data$colour,
            from = levels(full_spec_data$colour),
            to = c('Black',
                   'Blue',
                   'Clear',
                   'Green',
                   'Orange',
                   'Pink',
                   'Purple',
                   'Red',
                   'Yellow',
                   'Brown',
                   'olive',
                   'yellow green',
                   'Black',
                   'Blue',
                   'Brown',
                   'Clear',
                   'Purple',
                   'Multi-colour',
                   'Multi-colour',
                   'White'))

summary(full_spec_data$colour)

summary(full_spec_data$shape)

full_spec_data$shape <- mapvalues(full_spec_data$shape,
                                 from = levels(full_spec_data$shape),
                                 to = c('Fibre',
                                        'Fibre',
                                        'Film',
                                        'Fragment',
                                        'Fibre'))
summary(full_spec_data$shape)

summary(full_spec_data$raman.ID)

full_spec_data$raman.ID <- 
  mapvalues(full_spec_data$raman.ID,
            from = c('',
                     'Acrylic\n',
                     'Cellulosee',
                     'Cellulose\n',
                     'Unknown\n',
                     'Nylon\n',
                     'Polyester\n',
                     'Wool\n',
                     'Paint Chip'),
            to = c(NA,
                   'Acrylic',
                   'Cellulose',
                   'Cellulose',
                   'Unknown',
                   'Nylon',
                   'Polyester',
                   'Wool',
                   'Paint'))

full_spec_data$particle.type <- 
  mapvalues(full_spec_data$raman.ID,
            from = levels(full_spec_data$raman.ID),
            to = c('Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Natural',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Natural',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Natural'))

full_spec_data$particle.type[is.na(full_spec_data$particle.type)] <-
  with(subset(full_spec_data, is.na(particle.type)), 
       ifelse(length == 'NA',
              NA,
              'Unknown'))

summary(full_spec_data$particle.type)

full_spec_data$particle.type[full_spec_data$shape == 'Fibre' &
                              full_spec_data$colour == 'Clear' &
                              full_spec_data$raman.ID == 'Cellulose'] <-
  'Natural'  # call clear cellulosic fibres natural to be safe

summary(full_spec_data$particle.type)


full_spec_data$num <-
  ifelse(is.na(full_spec_data$length),
         0, 1)  # Separate out samples that had 0 counts



#### Estimate the particle type for unknown particles ####

## Construct a multinomial regression model for known particle types

moddata1 <- rbind(full_spec_data, PT[c(1:14, 22:23)], PJ)

moddata2 <- subset(moddata,
                   particle.type != 'Unknown' &
                     particle.type != 'NA')

moddata2$particle.type <- as.character(moddata2$particle.type)
moddata2$particle.type <- as.factor(moddata2$particle.type)

## Split data into training/test data (70/30)

moddata2$row.number <- as.numeric(rownames(moddata2))
moddata2 <- ungroup(moddata2)

set.seed(123)
train <- sample_frac(moddata2, size = 0.7)
moddata2$subset <- moddata2$row.number %in% train$row.number
test <- subset(moddata2, subset != 'TRUE')

## Fit intial model

particle.mod1 <-
  multinom(particle.type ~
             colour*shape*user.id + sample.type,
           data = train)

## Test model accuracy by building a classification table

# Predicting the values for train dataset
train$predicted <- predict(particle.mod1, newdata = train, "class")

# Building classification table
ctable1 <- table(train$particle.type, train$predicted)

# Calculating accuracy - sum of diagonal elements divided by total obs
round((sum(diag(ctable1))/sum(ctable1))*100,2)

## Repeat for test data

test$predicted <- predict(particle.mod1, newdata = test, "class")

# Building classification table
ctable2 <- table(test$particle.type, test$predicted)

# Calculating accuracy - sum of diagonal elements divided by total obs
round((sum(diag(ctable2))/sum(ctable2))*100,2)

## Refit to all data

particle.mod2 <- multinom(particle.type ~
                            colour*shape*user.id + sample.type,
                          data = moddata2)

predPT_data <- subset(PT, particle.type == 'Unknown')

predPJ_data <- subset(PJ, particle.type == 'Unknown')

pred_data <- subset(full_spec_data,
                    particle.type == 'Unknown')

## predict

predictPT.mlr <- predict(particle.mod2, newdata = predPT_data, type = 'class')

predictPJ.mlr <- predict(particle.mod2, newdata = predPJ_data, type = 'class')

predict.mlr <- predict(particle.mod2, newdata = pred_data, type = 'class')


## add predicitons back to original data

PT[PT$particle.type == 'Unknown' &
     !is.na(PT$particle.type), ]$particle.type <- predictPT.mlr

PJ[PJ$particle.type == 'Unknown' &
     !is.na(PJ$particle.type),]$particle.type <- predictPJ.mlr

full_spec_data$particle.type[full_spec_data$particle.type == 'Unknown' &
                               !is.na(full_spec_data$particle.type)] <-
  predict.mlr

PT$particle.type <- as.character(PT$particle.type)
PT$particle.type <- as.factor(PT$particle.type)

PJ$particle.type <- as.character(PJ$particle.type)
PJ$particle.type <- as.factor(PJ$particle.type)

full_spec_data$particle.type <- as.character(full_spec_data$particle.type)
full_spec_data$particle.type <- as.factor(full_spec_data$particle.type)

#### Blank subtract ####

## Summarize by particle type

PT_particle_type <-
  PT %>%
  group_by(ID,
           site,
           shape,
           sample.type,
           size.fraction,
           colour,
           blank.match,
           particle.type,
           sample.volume) %>%
  summarize(count = sum(num))  # plankton tows

PJ_particle_type <-
  PJ %>%
  group_by(ID,
           site,
           shape,
           sample.type,
           size.fraction,
           colour,
           blank.match,
           particle.type) %>%
  summarize(count = sum(num))  # plankton jars

animals_particle_type <-
  full_spec_data %>%
  group_by(ID,
           site,
           shape,
           sample.type,
           size.fraction,
           colour,
           blank.match,
           particle.type) %>%
  summarize(count = sum(num))  # animal data

## Separate blanks data

PT_blanks <- subset(PT_particle_type, sample.type == 'Blanks')
PT_data <- subset(PT_particle_type, sample.type != 'Blanks')

PJ_blanks <- subset(PJ_particle_type, sample.type == 'Blanks')
PJ_data <- subset(PJ_particle_type, sample.type != 'Blanks')

animal_blanks <- subset(animals_particle_type, sample.type == 'Blanks')
animal_data <- subset(animals_particle_type, sample.type != 'Blanks')

## Summarize blanks data

PT_blanks_means <- 
  PT_blanks %>% 
  group_by(shape, colour, blank.match, particle.type) %>% 
  summarize(blank.mean = mean(count))

PJ_blanks_means <- 
  PJ_blanks %>% 
  group_by(shape, colour, blank.match, particle.type) %>% 
  summarize(blank.mean = mean(count))

animal_blanks_means <- 
  animal_blanks %>% 
  group_by(shape, colour, blank.match, particle.type) %>% 
  summarize(blank.mean = mean(count))

## Blank subtract

# plankton tows

PT_data2 <-
  left_join(PT_data, 
            PT_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'particle.type'))

PT_data2$blank.mean[is.na(PT_data2$blank.mean)] <- 0

PT_data2$adj.count <- ceiling(with(PT_data2, count - blank.mean))
PT_data2$adj.count[PT_data2$adj.count < 0] <- 0

# plankton jars

PJ_data2 <-
  left_join(PJ_data, 
            PJ_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'particle.type'))

PJ_data2$blank.mean[is.na(PJ_data2$blank.mean)] <- 0

PJ_data2$adj.count <- ceiling(with(PJ_data2, count - blank.mean))
PJ_data2$adj.count[PJ_data2$adj.count < 0] <- 0

# animal data

animal_data2 <-
  left_join(animal_data, 
            animal_blanks_means, 
            by = c('shape', 'colour', 'blank.match', 'particle.type'))

animal_data2$blank.mean[is.na(animal_data2$blank.mean)] <- 0

animal_data2$adj.count <- ceiling(with(animal_data2, count - blank.mean))
animal_data2$adj.count[animal_data2$adj.count < 0] <- 0

## Combine samples that were split

animal_data2$ID <- mapvalues(
  animal_data2$ID,
  from = c(
    'CBSS3 1/2',
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
    'EBSS11 (3/3)',
    'EBRF16A',
    'EBRF16B',
    'HPRF21A',
    'HPRF21B',
    'HPRF4B'
  ),
  to = c(
    'CBSS3',
    'CBSS3',
    'CBSS8',
    'CBSS8',
    'CBSS8',
    'CBSS8',
    'EBSS10',
    'EBSS10',
    'EBSS10',
    'EBSS11',
    'EBSS11',
    'EBSS11',
    'EBRF16',
    'EBRF16',
    'HPRF21',
    'HPRF21',
    'HPRF4'
  )
)

## Remove incomplete samples

animal_incomplete <- 
  animal_data2 %>% 
  group_by(ID, sample.type) %>%
  filter(sample.type != 'Mussels' &
           sample.type != 'Clams'&
           sample.type != 'Surfperch Livers' &
           sample.type != 'Flatfish Livers' &
           sample.type != 'Rockfish Livers') %>% 
  summarize(size.fractions = length(unique(size.fraction))) %>% 
  filter(size.fractions < 2)

animal_data2$incomplete <- animal_data2$ID %in% animal_incomplete$ID

animal_data3 <-
  animal_data2 %>% 
  filter(incomplete == 'FALSE')

animal_data3$sample.type <- 
  mapvalues(animal_data3$sample.type,
            from = c('Surfperch Guts',
                     'Flatfish Guts',
                     'Rockfish Guts'),
            to = c('Surfperch',
                   'Flatfish',
                   'Rockfish'))

#### Plot breakdown by Raman ID, shape, and colour ####

moddata3 <- moddata1 %>% filter(num == 1)
moddata3$raman.ID[is.na(moddata3$raman.ID)] <- 'Unknown'
moddata3$raman.ID[moddata3$raman.ID == ''] <- 'Unknown'
moddata3$raman.ID <- as.character(moddata3$raman.ID)
moddata3$raman.ID <- as.factor(moddata3$raman.ID)

moddata3$site <- as.character(moddata3$site)
moddata3$site[moddata3$sample.type == 'Blanks'] <- 'Blanks'
moddata3$site <- as.factor(moddata3$site)

moddata3$raman.ID <- mapvalues(moddata3$raman.ID,
                               from = c('Polyster', 'Salt Crystal', 
                                        'Additive', 'Plastics Additive',
                                        'Plastics Dye', 'Unknown Polymer'),
                               to = c('Polyester', 'Salt', 'Unknown Synthetic',
                                      'Unknown Synthetic', 'Unknown Synthetic',
                                      'Unknown Synthetic'))

moddata3$raman.ID <- factor(moddata3$raman.ID,
                            levels = c('Unknown',
                                       'Acrylic',
                                       'Nylon',
                                       'Kevlar',
                                       'Paint',
                                       'Polyacrylonitrile',
                                       'Polyester',
                                       'Polypropylene',
                                       'Polystyrene',
                                       'Polyurethane',
                                       'Styrene copolymer',
                                       'Unknown Synthetic',
                                       'Rayon',
                                       'Cellulose',
                                       'Wool',
                                       'Bone',
                                       'Mineral',
                                       'Salt'))

moddata3$sample.type <- factor(moddata3$sample.type,
                               levels = c('Blanks',
                                          'Plankton Jars',
                                          'Plankton Tows',
                                          'Mussels',
                                          'Clams',
                                          'Sea Cucumbers',
                                          'Crabs',
                                          'Sea Stars',
                                          'Flatfish Guts',
                                          'Flatfish Livers',
                                          'Surfperch Guts',
                                          'Surfperch Livers',
                                          'Rockfish Guts',
                                          'Rockfish: Ingested Animals',
                                          'Rockfish Livers'))


tiff(
  'Polymer Plot.tiff',
  res = 300,
  width = 16.5,
  height = 11.43,
  units = 'cm',
  pointsize = 12
)

ggplot(moddata3) +
  geom_bar(
    aes(x = sample.type,
        y = num,
        fill = raman.ID),
    colour = 'black',
    size = 0.25,
    position = 'fill',
    stat = 'identity'
  ) +
  labs(x = 'Sample Type', y = 'Proportion of Particles') +
  facet_grid(. ~ shape, scales = 'free_x', space = 'free_x') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values =
                      qualitative_hcl(n = length(unique(moddata3$raman.ID)),
                                      palette = 'Dark3')) +
  guides(fill = guide_legend(ncol = 1)) +
  theme1 +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 6),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(0.75, 'line')
  )

dev.off()

## And by colour

summary(moddata3$colour)
moddata3$colour <- mapvalues(
  moddata3$colour,
  from = c(
    'olive',
    'yellow green',
    'black',
    'blue',
    'brown',
    'clear',
    'green',
    'orange',
    'pink',
    'purple',
    'red',
    'yellow',
    'rainbow'
  ),
  to = c(
    'Olive',
    'Yellow Green',
    'Black',
    'Blue',
    'Brown',
    'Clear',
    'Green',
    'Orange',
    'Pink',
    'Purple',
    'Red',
    'Yellow',
    'Multi-colour'
  )
)
summary(moddata3$colour)

tiff(
  'Colour Plot.tiff',
  res = 300,
  width = 16.5,
  height = 11.43,
  units = 'cm',
  pointsize = 12
)

ggplot(moddata3) +
  geom_bar(
    aes(x = sample.type,
        y = num,
        fill = colour),
    colour = 'black',
    size = 0.25,
    position = 'fill',
    stat = 'identity'
  ) +
  labs(x = 'Sample Type', y = 'Proportion or Particles') +
  facet_grid(. ~ shape, scales = 'free_x', space = 'free_x') +
  scale_fill_manual(
    values = c(
      'Black',
      'Steel Blue',
      'Azure2',
      'Forest Green',
      'Orange',
      'Pink',
      'Blue Violet',
      'Red3',
      'Gold',
      'Saddle Brown',
      'Dark Olive Green',
      'Yellow Green',
      'Turquoise',
      'White'
    )
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(ncol = 1)) +
  theme1 +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 6),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(0.75, 'line')
  )

dev.off()
