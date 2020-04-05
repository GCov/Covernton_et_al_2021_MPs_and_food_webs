##### Setup #####

## load packages
library(plyr)
library(ggplot2)
library(dplyr)
library(ggridges)
library(colorspace)

## Define standard error function

se <- function(x) {
  sd(x)/sqrt(length(x))
}

#### Combine all data ####

## Combine lab and field data with isotopes data

foodweb1 <- 
  left_join(animal_info, 
            isotopes,
            by = c('ID', 'site', 'sample.type'))



## Combine particle counts data with lab, field, and isotopes data

foodweb2 <- left_join(animal_data4,
                      foodweb1,
                      by = c('ID', 'site', 'sample.type'))

foodweb2$ID <- as.factor(foodweb2$ID)
foodweb2$site <- as.factor(foodweb2$site)
foodweb2$sample.type <- as.factor(foodweb2$sample.type)
foodweb2shape <- as.factor(foodweb2$shape)
foodweb2$colour <- as.factor(foodweb2$colour)

#### Plot ####

ggplot(foodweb2) +
  geom_point(aes(x = trophic.position,
                 y = count,
                 colour = reorder(species, trophic.position, max))) +
  facet_grid(particle.type ~ site) +
  theme1
             