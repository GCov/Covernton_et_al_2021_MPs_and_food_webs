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

foodweb1$sample.type <- as.factor(foodweb1$sample.type)

ingested_animals <- 
  foodweb1 %>% 
  filter(sample.type == 'Rockfish: Ingested Animals') %>% 
  select(1:14)

rockfish_info <- 
  foodweb1 %>% 
  filter(animal.type == 'Rockfish' & sample.type == 'Rockfish') %>% 
  select(1,2,15:36)

ingested_animals2 <-
  left_join(ingested_animals, rockfish_info, by = c('ID', 'site'))

foodweb1 <- rbind(subset(foodweb1, 
                         sample.type != 'Rockfish: Ingested Animals'),
                  ingested_animals2)

## Combine particle counts data with lab, field, and isotopes data

foodweb2 <- left_join(animal_data3,
                      foodweb1,
                      by = c('ID', 'site', 'sample.type'))

foodweb2$ID <- as.factor(foodweb2$ID)
foodweb2$site <- as.factor(foodweb2$site)
foodweb2$sample.type <- as.factor(foodweb2$sample.type)
foodweb2shape <- as.factor(foodweb2$shape)
foodweb2$colour <- as.factor(foodweb2$colour)

#### Summarize ####

## Combine across size categories, colour, and shape

foodweb3 <- 
  foodweb2 %>% 
  group_by(ID, site, sample.type, particle.type, shell.l, shell.w, shell.h,
           arm.length, tissue.wet.weight, tissue.dry.weight, shell.weight,
           total.body.wet.weight, density.sep, species, carapace.length,
           TL, SL, sex, babies, parasites, deltaC, deltaN, trophic.position) %>% 
  summarize(count = sum(adj.count))

## Make sure each sample has all levels of particle type

all_types <- expand.grid(ID = levels(foodweb3$ID), 
                         particle.type = levels(foodweb3$particle.type))

info <- foodweb3[c(1:3, 5:23)]

all_types2 <- left_join(all_types, info, by = 'ID')

counts <- subset(foodweb3[c(1,4,24)], count > 0)

foodweb4 <- left_join(all_types2, counts, by = c('ID', 'particle.type'))

foodweb4$count[is.na(foodweb4$count)] <- 0

summary(foodweb4)

foodweb4 <- subset(foodweb4, !is.na(species))

gutdata <- subset(foodweb4,
                  sample.type != 'Surfperch Livers' &
                  sample.type != 'Flatfish Livers' &
                  sample.type != 'Rockfish Livers')

gutdata$sample.type <- 
  mapvalues(gutdata$sample.type,
            from = 'Rockfish: Ingested Animals',
            to = 'Rockfish')

gutdata <-
  gutdata %>% 
  group_by(ID, particle.type, site, sample.type,
           total.body.wet.weight, species, deltaC, deltaN,
           trophic.position) %>% 
  summarize(count = sum(count))

#### Plot ####

tiff('Trophic Position MP Plot.tiff',
     res = 300,
     width = 17,
     height = 14,
     units = 'cm',
     pointsize = 12)

ggplot(subset(gutdata, !is.na(trophic.position))) +
  geom_smooth(aes(x = trophic.position,
                  y = log(count + 1)),
              method = 'lm',
              colour = 'red',
              alpha = 0.5,
              size = 0.5) +
  geom_point(aes(x = trophic.position,
                 y = log(count + 1),
                 colour = reorder(species, trophic.position, mean)),
             size = 0.5) +
  facet_grid(particle.type ~ site) +
  labs(x = 'Trophic Position',
       y = 'ln(MPs + 1) (particles/ind)') +
  scale_colour_manual(values = sequential_hcl(palette = 'Viridis', n = 14)) +
  theme1

dev.off()

tiff('Animal Type Plot.tiff',
     res = 300,
     width = 16,
     height = 8,
     units = 'cm',
     pointsize = 12)

ggplot(subset(gutdata, particle.type == 'Synthetic Polymer')) +
  geom_boxplot(aes(x = sample.type,
                   y = count),
               size = 0.5) +
  facet_grid(. ~ site, scales = 'free_x') +
  labs(x = 'Type of Animal',
       y = 'MPs (particles/ind)') +
  theme1 + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
