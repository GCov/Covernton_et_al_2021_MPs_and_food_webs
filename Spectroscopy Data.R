##### Setup #####

## load packages
library(plyr)
library(ggplot2)
library(dplyr)
library(nnet)
library(MuMIn)
library(boot)
library(colorspace)
library(randomForest)
require(caTools)
library(ROCR)
library(extrafont)
loadfonts(device = "win")

#### Load spectroscopy data ####

## Load plankton tow data

plankton_tows <- read.csv("plankton_tows.csv", header = TRUE)
PT_field <- read.csv("PT_field_data.csv", header = TRUE)

## Clean up and combine plankton tow field data 

summary(PT_field$sample.volume)

PT <- left_join(plankton_tows, PT_field, by = 'ID')

names(PT)
PT$ID <- as.factor(PT$ID)
PT$site <- as.factor(PT$site)
PT$sample.type <- as.factor(PT$sample.type)
PT$size.fraction <- as.factor(PT$size.fraction)
PT$shape <- as.factor(PT$shape)
PT$colour <- as.factor(PT$colour)
PT$raman.ID <- as.factor(PT$raman.ID)
PT$user.id <- as.factor(PT$user.id)

## Clean up plankton tow data

names(PT)
summary(PT$size.fraction)
summary(PT$shape)
summary(PT$colour)

summary(PT$shape)

PT$colour <- mapvalues(
  PT$colour,
  from = levels(PT$colour),
  to = c("Black",
         "Blue",
         "Brown",
         "Clear",
         "Green",
         "Orange",
         "Pink",
         "Purple",
         "Red",
         "Yellow"))

PT$shape <- mapvalues(
  PT$shape,
  from = c('fibre',
           'fibre ',
           'fragment'),
  to = c('Fibre',
         'Fibre',
         'Fragment')
)
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

## Load plankton jars data

PJ <- read.csv("plankton_jars.csv", header = TRUE)  # plankton jars

PJ$ID <- as.factor(PJ$ID)
PJ$site <- as.factor(PJ$site)
PJ$sample.type <- as.factor(PJ$sample.type)
PJ$size.fraction <- as.factor(PJ$size.fraction)
PJ$shape <- as.factor(PJ$shape)
PJ$colour <- as.factor(PJ$colour)
PJ$raman.ID <- as.factor(PJ$raman.ID)
PJ$user.id <- as.factor(PJ$user.id)

## Clean up plankton jars data

names(PJ)
summary(PJ$size.fraction)
summary(PJ$shape)
summary(PJ$colour)

summary(PJ$shape)

PJ$colour <- mapvalues(PJ$colour,
                       from = levels(PJ$colour),
                       to = c("Black",
                              "Blue",
                              "Brown",
                              "Clear",
                              "Green",
                              "Orange",
                              "Pink",
                              "Purple",
                              "Multi-colour",
                              "Red"
                              ))


PJ$shape <- mapvalues(
  PJ$shape,
  from = c('fiber',
           'fibre',
           'fragment'),
  to = c('Fibre',
         'Fibre',
         'Fragment')
)
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
full_spec_data$ID <- as.factor(full_spec_data$ID)
full_spec_data$site <- as.factor(full_spec_data$site)
full_spec_data$sample.type <- as.factor(full_spec_data$sample.type)
full_spec_data$size.fraction <- as.factor(full_spec_data$size.fraction)
full_spec_data$shape <- as.factor(full_spec_data$shape)
full_spec_data$colour <- as.factor(full_spec_data$colour)
full_spec_data$raman.ID <- as.factor(full_spec_data$raman.ID)
full_spec_data$user.id <- as.factor(full_spec_data$user.id)

summary(full_spec_data$blank.match)
summary(full_spec_data$size.fraction)
summary(full_spec_data$shape)
summary(full_spec_data$colour)

full_spec_data$colour <-
  mapvalues(full_spec_data$colour,
            from = levels(full_spec_data$colour),
            to = c('Black',
                   'Black',
                   'Blue',
                   'Blue',
                   'Brown',
                   'Brown',
                   'Clear',
                   'Clear',
                   'Green',
                   'Multi-colour',
                   'Olive',
                   'Orange',
                   'Pink',
                   'Purple',
                   'Purple',
                   'Multi-colour',
                   'Red',
                   'White',
                   'Yellow',
                   'Yellow-green'))

summary(full_spec_data$colour)

summary(full_spec_data$shape)

full_spec_data$shape <- mapvalues(full_spec_data$shape,
                                 from = levels(full_spec_data$shape),
                                 to = c('Fibre',
                                        'Fibre',
                                        'Fibre',
                                        'Film',
                                        'Fragment'))
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
                     'Paint Chip',
                     'Polyster'),
            to = c(NA,
                   'Acrylic',
                   'Cellulose',
                   'Cellulose',
                   'Unknown',
                   'Nylon',
                   'Polyester',
                   'Wool',
                   'Paint',
                   'Polyester'))

full_spec_data$particle.type <- 
  mapvalues(full_spec_data$raman.ID,
            from = levels(full_spec_data$raman.ID),
            to = c('Synthetic Polymer',
                   'Synthetic Polymer',
                   'Natural',
                   'Natural Anthropogenic',
                   'Synthetic Polymer',
                   'Natural',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Synthetic Polymer',
                   'Semi-synthetic',
                   'Natural',
                   'Synthetic Polymer',
                   'Unknown',
                   'Synthetic Polymer',
                   'Natural Anthropogenic'))

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


#### Estimate the particle type for unknown particles ####

## Create random forest model for known particle types

moddata1 <- rbind(full_spec_data, PT[c(1:14, 22)], PJ)

moddata2 <- subset(moddata1,
                   particle.type != 'Unknown' &
                     particle.type != 'NA')

moddata2$particle.type <- as.character(moddata2$particle.type)
moddata2$particle.type <- as.factor(moddata2$particle.type)

## Run model on test data

rf <- randomForest(particle.type ~
                     colour + shape + user.id + sample.type + site + length,
                   data = moddata2)

rf

mtry <- tuneRF(moddata2[, c(2:3,7:9,14)], 
               moddata2$particle.type,
               stepFactor = 1.5, 
               improve = 0.01,
               trace = TRUE,
               plot = TRUE)

rf <- randomForest(particle.type ~
                     colour * shape * user.id + sample.type + site + length,
                   data = moddata2,
                   mtry = 2,
                   importance = TRUE,
                   ntree = 2000)

plot(rf$err.rate[,1], type = 'l')  ## find ideal ntree value 

rf <- randomForest(particle.type ~
                     colour * shape * user.id + sample.type + site + length,
                   data = moddata2,
                   mtry = 2,
                   importance = TRUE,
                   ntree = 250)

rf
importance(rf)

varImpPlot(rf, main = "")


pred <- predict(rf)

## Performance metrics

cm <- table(moddata2[, 15], pred)

round((sum(diag(cm))/sum(cm))*100,2)

plot_rf_confusion <- function(rf_model)
{
  conf_mat <- rf_model$confusion[, -ncol(rf_model$confusion)]
  class_size <- apply(conf_mat, 1, sum)
  
  conf_mat_perc <- (conf_mat / class_size) * 100
  
  melt_conf_mat <- reshape2::melt(conf_mat_perc, na.rm = TRUE)
  
  melt_conf_mat_raw <- reshape2::melt(conf_mat, na.rm = TRUE)
  
  melt_conf_mat$value2 <- melt_conf_mat_raw$value
  
  conf_plot <-
    ggplot(data = melt_conf_mat, aes_string('Var2', 'Var1', fill = 'value')) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = pal[3],
      high = pal[2],
      mid = pal[4],
      midpoint = 50,
      limit = c(0, 100),
      space = "Lab",
      name = "Accuracy (%)"
    ) +
    geom_text(aes(label = value2), size = 4, nudge_y = 0.2) +
    geom_text(aes(label = sprintf("%.1f", round(value, 1))), size = 4,
              nudge_x = -0.1, nudge_y = -0.2) +
    geom_text(label = "%", nudge_x = 0.3, nudge_y = -0.2, size = 4) +
    theme_minimal() +
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      size = 10,
      hjust = 1
    ),
    plot.margin = unit(c(0,0,0,0.5), "cm")) +
    theme(axis.text.y = element_text(size = 10)) +
    labs(x = 'Predicted', y = 'True') +
    coord_fixed()
  
  return(conf_plot)
}

png("Confusion Matrix.png",
    width = 14,
    height = 14, 
    units = "cm",
    res = 700)

plot_rf_confusion(rf)

dev.off()

## Now predict unknown particle types

predPT_data <- subset(PT[c(1:14, 22)], 
                      particle.type == 'Unknown')
predPT_data <- rbind(train[1, -c(16:17)] , predPT_data)
predPT_data <- predPT_data[-1,]

predPJ_data <- subset(PJ,  
                      particle.type == 'Unknown')
predPJ_data <- rbind(train[1, -c(16:17)] , predPJ_data)
predPJ_data <- predPJ_data[-1,]

predanimal_data <- subset(full_spec_data,
                    particle.type == 'Unknown')
predanimal_data <- rbind(train[1, -c(16:17)] , predanimal_data)
predanimal_data <- predanimal_data[-1,]

predictPT.rf <- predict(rf, newdata = predPT_data)

predictPJ.rf <- predict(rf, newdata = predPJ_data)

predictanimal.rf <- predict(rf, newdata = predanimal_data)


## add predictions back to original data

PT[PT$particle.type == 'Unknown' &
     !is.na(PT$particle.type), ]$particle.type <- predictPT.rf

PJ[PJ$particle.type == 'Unknown' &
     !is.na(PJ$particle.type),]$particle.type <- predictPJ.rf

full_spec_data$particle.type[full_spec_data$particle.type == 'Unknown' &
                               !is.na(full_spec_data$particle.type)] <-
  predictanimal.rf

PT$particle.type <- as.character(PT$particle.type)
PT$particle.type <- as.factor(PT$particle.type)

PJ$particle.type <- as.character(PJ$particle.type)
PJ$particle.type <- as.factor(PJ$particle.type)

full_spec_data$particle.type <- as.character(full_spec_data$particle.type)
full_spec_data$particle.type <- as.factor(full_spec_data$particle.type)

#### Make sure each sample has all levels of particle type ####

## plankton tows

PT$ID <- as.character(PT$ID)
PT$ID <- as.factor(PT$ID)

all_typesPT <- expand.grid(ID = levels(PT$ID),
                           particle.type = levels(PT$particle.type))

infoPT <- PT[c(1:3, 6, 13, 21)] %>% 
  group_by(ID, site, sample.type, size.fraction, blank.match, sample.volume) %>% 
  summarize()

all_typesPT2 <- left_join(all_typesPT, infoPT, by = c("ID"))

countsPT <- PT[ ,c(1, 4:12, 14, 22)]

PT2 <- left_join(all_typesPT2,
                 countsPT,
                 by = c('ID', 'particle.type', 'size.fraction'))

## plankton jars
PJ$ID <- as.character(PJ$ID)
PJ$ID <- as.factor(PJ$ID)

all_typesPJ <- expand.grid(ID = levels(PJ$ID),
                           particle.type = levels(PJ$particle.type))

infoPJ <- PJ[c(1:3, 6, 13)] %>% 
  group_by(ID, site, sample.type, size.fraction, blank.match) %>% 
  summarize()

all_typesPJ2 <- left_join(all_typesPJ, infoPJ, by = 'ID')

countsPJ <- PJ[ ,c(1, 4:12, 14:15)]

PJ2 <- left_join(all_typesPJ2,
                 countsPJ,
                 by = c('ID', 'particle.type', 'size.fraction'))


## animal samples

full_spec_data$ID <- as.character(full_spec_data$ID)
full_spec_data$ID <- as.factor(full_spec_data$ID)

full_spec_data$sample.type <- as.character(full_spec_data$sample.type)
full_spec_data$sample.type <- as.factor(full_spec_data$sample.type)

all_types <- expand.grid(ID = levels(full_spec_data$ID), 
                         particle.type = levels(full_spec_data$particle.type))

info <- 
  full_spec_data[ ,c(1:3, 5:6, 13)] %>% 
  group_by(ID, site, sample.type, size.fraction, blank.match) %>% 
  summarize()

all_types2 <- left_join(all_types, info, by = 'ID')

counts <- full_spec_data[ ,c(1, 4:12, 14:15)]

full_spec_data2 <- left_join(all_types2,
                             counts,
                             by = c('ID',
                                    'particle.type', 
                                    'size.fraction'))

summary(full_spec_data2$sample.type)

full_spec_data2$sample.type <-
  mapvalues(full_spec_data2$sample.type,
            from = c("Rockfish Gut",
                     "Rockfish Gut Animals"),
            to = c("Rockfish Guts",
                   "Rockfish: Ingested Animals"))

#### Blank subtract ####

## Summarize by particle type

PT2$num <- ifelse(is.na(PT2$length), 0, 1)

PT_particle_type <-
  PT2 %>%
  group_by(ID,
           site,
           sample.type,
           blank.match,
           particle.type,
           sample.volume) %>%
  summarize(count = sum(num))  # plankton tows

PJ2$num <- ifelse(is.na(PJ2$length), 0, 1)

PJ_particle_type <-
  PJ2 %>%
  group_by(ID,
           site,
           sample.type,
           blank.match,
           particle.type) %>%
  summarize(count = sum(num))  # plankton jars

full_spec_data2$num <- ifelse(is.na(full_spec_data2$length), 0, 1)

animals_particle_type <-
  full_spec_data2 %>%
  group_by(ID,
           site,
           sample.type,
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
  group_by(blank.match, particle.type) %>% 
  summarize(blank.mean = mean(count))

PJ_blanks_means <- 
  PJ_blanks %>% 
  group_by(blank.match, particle.type) %>% 
  summarize(blank.mean = mean(count))

animal_blanks_means <- 
  animal_blanks %>% 
  group_by(blank.match, particle.type) %>% 
  summarize(blank.mean = mean(count))

## Blank subtract

# plankton tows

PT_data2 <-
  left_join(PT_data, 
            PT_blanks_means, 
            by = c('blank.match', 'particle.type'))

PT_data2$blank.mean[is.na(PT_data2$blank.mean)] <- 0

PT_data2$adj.count <- floor(with(PT_data2, count - blank.mean))
PT_data2$adj.count[PT_data2$adj.count < 0] <- 0

# plankton jars

PJ_data2 <-
  left_join(PJ_data, 
            PJ_blanks_means, 
            by = c('blank.match', 'particle.type'))

PJ_data2$blank.mean[is.na(PJ_data2$blank.mean)] <- 0

PJ_data2$adj.count <- floor(with(PJ_data2, count - blank.mean))
PJ_data2$adj.count[PJ_data2$adj.count < 0] <- 0

# animal data

animal_data2 <-
  left_join(animal_data, 
            animal_blanks_means, 
            by = c('blank.match', 'particle.type'))

animal_data2$blank.mean[is.na(animal_data2$blank.mean)] <- 0

animal_data2$adj.count <- floor(with(animal_data2, count - blank.mean))
animal_data2$adj.count[animal_data2$adj.count < 0] <- 0

## Combine samples that were split

animal_data2$ID <- mapvalues(
  animal_data2$ID,
  from = c(
    'CBSS3 1/2',
    'CBSS3 2/2',
    'CBSS6-1',
    'CBSS6-2',
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
    'CBSS6',
    'CBSS6',
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

# animal_incomplete <- 
#   animal_data2 %>% 
#   group_by(ID, sample.type) %>%
#   filter(sample.type != 'Mussels' &
#            sample.type != 'Clams'&
#            sample.type != 'Surfperch Livers' &
#            sample.type != 'Flatfish Livers' &
#            sample.type != 'Rockfish Livers') %>% 
#   summarize(size.fractions = length(unique(size.fraction))) %>% 
#   filter(size.fractions < 2)
# 
# animal_data2$incomplete <- animal_data2$ID %in% animal_incomplete$ID
# 
# animal_data3 <-
#   animal_data2 %>% 
#   filter(incomplete == 'FALSE')
# 
# animal_data3$sample.type <- 
#   mapvalues(animal_data3$sample.type,
#             from = c('Surfperch Guts',
#                      'Flatfish Guts',
#                      'Rockfish Guts'),
#             to = c('Surfperch',
#                    'Flatfish',
#                    'Rockfish'))

#### Plot breakdown by Raman ID, shape, and colour ####

moddata1$num <- ifelse(is.na(moddata1$shape), 0, 1)
moddata3 <- moddata1 %>% filter(num == 1)
moddata3$raman.ID[is.na(moddata3$raman.ID)] <- 'Unknown'
moddata3$raman.ID[moddata3$raman.ID == ''] <- 'Unknown'
moddata3$raman.ID <- as.character(moddata3$raman.ID)
moddata3$raman.ID <- as.factor(moddata3$raman.ID)

moddata3$site <- as.character(moddata3$site)
moddata3$site[moddata3$sample.type == 'Blanks'] <- 'Blanks'
moddata3$site <- as.factor(moddata3$site)

moddata3$raman.ID <- mapvalues(moddata3$raman.ID,
                               from = c('Salt Crystal', 
                                        'Additive', 'Plastics Additive',
                                        'Plastics Dye', 'Unknown Polymer'),
                               to = c('Salt', 'Unknown Synthetic',
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
  res = 700,
  width = 19,
  height = 15,
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
  facet_grid(particle.type ~ shape, 
             scales = 'free_x', 
             space = 'free_x',
             labeller = label_wrap_gen(width = 10)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values =
                      sequential_hcl(n = 18, palette = "Viridis")) +
  guides(fill = guide_legend(ncol = 1)) +
  theme1 +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.75, 'line'),
    panel.spacing = unit(0.25, "cm")
  )

dev.off()

## And by colour

summary(moddata3$colour)
moddata3$colour <- mapvalues(
  moddata3$colour,
  from = c(
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
  res = 700,
  width = 19,
  height = 15,
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
  labs(x = 'Sample Type', y = 'Proportion of Particles') +
  facet_grid(particle.type ~ shape, 
             scales = 'free_x', 
             space = 'free_x',
             labeller = label_wrap_gen(width = 10)) +
  scale_fill_manual(
    values = c(
      'Grey16',
      'Steel Blue',
      'Saddle Brown',
      'Azure2',
      'Forest Green',
      'Turquoise',
      'Dark Olive Green',
      'Orange',
      'Pink',
      'Blue Violet',
      'Red3',
      'White',
      'Gold',
      'Yellow Green'
    )
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(ncol = 1)) +
  theme1 +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.75, 'line'),
    panel.spacing = unit(0.25, "cm")
  )

dev.off()


#### Combine food web data ####

## Combine lab and field data with isotopes data

foodweb1 <- 
  left_join(animal_info, 
            isotopes2,
            by = c('ID', 'site'))

foodweb1$sample.type <- as.factor(foodweb1$sample.type)

ingested_animals <- 
  foodweb1 %>% 
  filter(sample.type == 'Rockfish: Ingested Animals') %>% 
  select(1:14)

rockfish_info <- 
  foodweb1 %>% 
  filter(animal.type == 'Rockfish' & sample.type == 'Rockfish') %>% 
  select(1,2,15:30)

ingested_animals2 <-
  left_join(ingested_animals, rockfish_info, by = c('ID', 'site'))

foodweb1 <- rbind(subset(foodweb1, 
                         sample.type != 'Rockfish: Ingested Animals'),
                  ingested_animals2)

## Combine particle counts data with lab, field, and isotopes data

animal_data2$sample.type <- 
  mapvalues(animal_data2$sample.type,
            from = c("Flatfish Guts",
                     "Rockfish Guts",
                     "Surfperch Guts"),
            to = c("Flatfish",
                   "Rockfish",
                   "Surfperch"))

animal_data3 <-
  animal_data2 %>%
  group_by(ID,
           site,
           sample.type,
           blank.match,
           particle.type) %>%
  summarize(
    count = sum(count),
    blank.mean = sum(blank.mean),
    adj.count = sum(adj.count)
  )

foodweb2 <- left_join(animal_data3,
                      foodweb1,
                      by = c('ID', 'site', 'sample.type'))

foodweb2$ID <- as.factor(foodweb2$ID)
foodweb2$site <- as.factor(foodweb2$site)
foodweb2$sample.type <- as.factor(foodweb2$sample.type)
  



## Now summarize animal data

summary(subset(foodweb2, is.na(species)))
foodweb2$species[foodweb2$ID == "CBSS16"] <- "Dermasterias imbricata"

gutdata <- subset(foodweb2,
                  sample.type != 'Surfperch Livers' &
                    sample.type != 'Flatfish Livers' &
                    sample.type != 'Rockfish Livers')

gutdata$sample.type <- 
  mapvalues(gutdata$sample.type,
            from = 'Rockfish: Ingested Animals',
            to = 'Rockfish')
gutdata <-
  gutdata %>% 
  group_by(ID, site, sample.type, blank.match, particle.type, shell.l,
           shell.w, shell.h, arm.length, tissue.wet.weight,
           shell.weight, total.body.wet.weight, density.sep, species, 
           carapace.length, TL, SL, hepatopancreas.weight, sex, babies, 
           parasites, base_deltaN, sd_base_deltaN, deltaN, deltaC, 
           trophic.position) %>% 
  summarize(count = sum(count),
            blank.mean = sum(blank.mean),
            adj.count = sum(adj.count),
            tissue.dry.weight = sum(tissue.dry.weight))

liverdata <-subset(foodweb2,
                   sample.type == 'Surfperch Livers' |
                     sample.type == 'Flatfish Livers' |
                     sample.type == 'Rockfish Livers')

