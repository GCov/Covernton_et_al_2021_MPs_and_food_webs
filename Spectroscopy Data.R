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

PT$particle.type[is.na(PT$particle.type)] <-
  with(subset(PT, is.na(particle.type)), 
       ifelse(length == 'NA',
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
            to = c('Unknown',
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

#### Estimate the particle type for unkown particles ####

### THIS ISN'T WORKING AND I DON'T KNOW WHY (maybe just CV?)

## Construct a multinomial regression model for known particle types

moddata <- subset(full_spec_data,
                  particle.type != 'Unknown' &
                    particle.type != 'NA')

moddata$particle.type <- as.character(moddata$particle.type)
moddata$particle.type <- as.factor(moddata$particle.type)

## Fit intial model

particle.mod1 <-
  multinom(particle.type ~
             raman.ID,
           data = moddata)
summary(particle.mod1)

## Cross validation

## Randomly shuffle the data

cvdata <- moddata[sample(nrow(moddata)),]

cvdata$blank.match <- as.factor(cvdata$blank.match)

## Create 10 equally size folds

folds <- cut(seq(1,nrow(cvdata)),breaks=10,labels=FALSE)

## Perform 10 fold cross validation

success_rate <- numeric()


for (i in 1:10) {
  testIndexes <- which(folds == i, arr.ind = TRUE)
  testData <- cvdata[testIndexes,]
  trainData <- cvdata[-testIndexes,]
  model <- multinom(particle.type ~
                      colour + shape + sample.type,
                    data = trainData)
  test <- predict(model,
                  newdata = subset(testData,
                                   particle.type == 'Synthetic Polymer'))
  successes <- ifelse(test == subset(testData,
                                     particle.type == 
                                       'Synthetic Polymer')$particle.type,
                      1,
                      0)
  success_rate[i] <- mean(successes)
}

mean(success_rate)
range(success_rate)

pred_data <- subset(full_spec_data,
                    particle.type == 'Unknown')

predict.mlr <- predict(particle.mod1, newdata = pred_data, type = 'class')

# ## 2nd choice particles for p<0.8 for a class
# 
# predict.probs <- 
#   as.data.frame(predict(particle.mod1, 
#                         newdata = pred_data, 
#                         type = 'probs'))
# 
# predict.mlr2 <- character()
# 
# for(i in 1:length(pred_data$particle.type)) {
#   predictions <- predict.probs[i,]
#   max <- max(predictions)
#   second <- max(predictions[predictions != max])
#   name <- ifelse(max >= 0.8,
#                  names(predictions)[which(predictions == max, 
#                                           arr.ind = T)[, "col"]],
#                  names(predictions)[which(predictions == second,
#                                           arr.ind = T)[, "col"]])
#   predict.mlr2[i] <- name
# }

full_spec_data$particle.type[full_spec_data$particle.type == 'Unknown' &
                               !is.na(full_spec_data$particle.type)] <-
  predict.mlr

full_spec_data$particle.type <- as.character(full_spec_data$particle.type)
full_spec_data$particle.type <- as.factor(full_spec_data$particle.type)

## See what the predictions were for each unknown particle

pred_data$predict <- as.factor(predict.mlr)

plot(predict ~ colour, data = pred_data)

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
