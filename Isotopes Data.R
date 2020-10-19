##### Setup #####

## load packages
library(plyr)
library(ggplot2)
library(dplyr)
library(colorspace)

#### Process Data ####

isotopes <- read.csv('isotopes.csv', header = TRUE)

names(isotopes)

summary(isotopes)

## Add species data

isotopes <- left_join(isotopes, 
                      animal_info[, c(1, 3, 16)],
                      by = c("ID", "sample.type"))
isotopes$sample.type <- as.factor(isotopes$sample.type)
isotopes$species <- as.character(isotopes$species)

isotopes$species[isotopes$sample.type == "Flatfish Livers" &
                   isotopes$species == "Parophrys vetulus"] <- 
  "Parophrys vetulus liver"

isotopes$species[isotopes$sample.type == "Flatfish Livers" &
                   isotopes$species == "Platichthys stellatus"] <-
  "Platichthys stellatus liver"

isotopes$species[isotopes$sample.type == "Rockfish Livers" &
                   isotopes$species == "Sebastes melanops"] <-
  "Sebastes melanops liver"

isotopes$species[isotopes$sample.type == "Rockfish Livers" &
                   isotopes$species == "Sebastes caurinus"] <-
  "Sebastes caurinus liver"

isotopes$species[isotopes$sample.type == "Surfperch Livers" &
                   isotopes$species == "Cymatogaster aggregata"] <-
  "Cymatogaster aggregata liver"

isotopes$species <- as.factor(isotopes$species)

summary(isotopes$species)

isotopes$species[isotopes$ID == "CBSS16"] <- "Dermasterias imbricata"
isotopes$species[isotopes$ID == "HPRF1"] <- "Sebastes caurinus"

summary(isotopes$species)

theme1 <-
  theme_bw() +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.title = element_blank()
  )


tiff('Isotopic_Biplot.tiff',width = 19, height = 19, units = 'cm', res = 300)

ggplot(isotopes) +  # isotopic plot
  geom_point(aes(x = deltaC,
                 y = deltaN,
                 colour = reorder(species, deltaN, mean)),
             size = 0.75) +
  facet_grid(~site) +
  scale_colour_manual(values = qualitative_hcl(n = 20, 
                                               h = c(-180, 160), 
                                               c = 60, 
                                               l = 75)) +
  theme1

dev.off()

## Estimate trophic position for each individual using mussels as a baseline

baseline <-
  subset(isotopes, sample.type == 'Mussels') %>% 
  group_by(site) %>% 
  summarize(base_deltaN = mean(deltaN),
            sd_base_deltaN = sd(deltaN))

baseline

isotopes <- left_join(isotopes,
                      baseline,
                      by = 'site')

isotopes$trophic.position <-
  with(isotopes,
       ((deltaN - base_deltaN)/3.4)+1)

## Average liver and muscle values for fish

isotopes2 <-
  isotopes %>% 
  group_by(ID, site, base_deltaN, sd_base_deltaN) %>% 
  summarize(deltaN = mean(deltaN))

isotopes2$trophic.position <-
  with(isotopes2,
       ((deltaN - base_deltaN)/3.4)+1)

## Plot according to trophic position

ggplot(isotopes) +
  geom_violin(aes(x = 1,
                  y = trophic.position,
                  fill = site),
              size = 1) +
  facet_grid(. ~ reorder(sample.type, 
                              trophic.position,
                              mean)) +
  labs(x = 'Sample Type',
       y = 'Trophic Position') +
  scale_fill_manual(values = qualitative_hcl(palette = 'Dark 2', n = 3)) +
  theme1 +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


