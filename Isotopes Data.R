##### Setup #####

## Load packages
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

#### Plot Isotopes Data ####

theme1 <-
  theme_bw() +
  theme(
    axis.text = element_text(size = 7,
                             family = "sans"),
    axis.title = element_text(size = 9,
                              family = "sans"),
    strip.background = element_blank(),
    strip.text = element_text(size = 8,
                              family = "sans"),
    legend.text = element_text(size = 10,
                               family = "sans"),
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

isotopes$animal.type <- isotopes$sample.type

isotopes$animal.type <-
  mapvalues(isotopes$animal.type,
            from = c("Flatfish Livers",
                     "Surfperch Livers",
                     "Rockfish Livers"),
            to = c("Flatfish",
                   "Surfperch",
                   "Rockfish"))

isotopes$tissue.type <-
  ifelse(isotopes$sample.type == "Flatfish Livers" |
           isotopes$sample.type == "Surfperch Livers" |
           isotopes$sample.type == "Rockfish Livers",
         "Liver",
         "Muscle")

site.lab <- c("Coles Bay" = "Coles Bay", 
              "Elliot Bay" = "Elliot Beach", 
              "Victoria Harbour" = "Victoria Harbour")

isotopes$include <- isotopes$ID %in% gutdata$ID

isotopes2 <- subset(isotopes, include == "TRUE")

isotopes2$species <- as.character(isotopes2$species)
isotopes2$species <- as.factor(isotopes2$species)

isotopes2$common.name <- mapvalues(isotopes2$species,
                                   from = levels(isotopes2$species),
                                   to = c("Red Rock Crab",
                                          "Orange Sea Cucumber",
                                          "Shiner Surfperch",
                                          "Shiner Surfperch",
                                          "Leather Star",
                                          "Graceful Rock Crab",
                                          "Dungeness Crab",
                                          "Blue Mussel",
                                          "California Sea Cucumber",
                                          "English Sole",
                                          "English Sole",
                                          "Starry Flounder",
                                          "Starry Flounder",
                                          "Littleneck Clam",
                                          "Manila Clam",
                                          "Copper Rockfish",
                                          "Copper Rockfish",
                                          "Black Rockfish",
                                          "Black Rockfish"))

isotopes2$cn.short <- mapvalues(isotopes2$common.name,
                                from = levels(isotopes2$common.name),
                                to = c("RR",
                                       "OC",
                                       "SS",
                                       "LS",
                                       "GR",
                                       "DU",
                                       "BM",
                                       "CC",
                                       "ES",
                                       "SF",
                                       "LC",
                                       "MC",
                                       "CR",
                                       "BR"))

tiff('Isotopic_Biplot.tiff', 
     height = 8,
     width = 6,
     units = "in",
     res = 800,
     compression = "lzw")

ggplot(isotopes2) +  # isotopic plot
  geom_text(aes(x = deltaC,
                y = deltaN,
                label = cn.short,
                colour = tissue.type),
            size = 6/.pt,
            alpha = 0.8,
            family = "sans") +
  facet_grid(site ~ .,
             labeller = labeller(site = site.lab)) +
  labs(x = expression(paste(delta^13*"C (\211)")),
       y = expression(paste(delta^15*"N (\211)"))) +
  scale_colour_manual(values = pal[c(1,3)]) +
  scale_x_continuous(breaks = c(seq(-22, -10, 2)),
                     expand = c(0,0.3)) +
  scale_y_continuous(breaks = seq(8, 18, 2),
                     expand = c(0,0.4)) +
  theme1 +
  theme(legend.position = "bottom")

dev.off()

## Calculate mean and sd for baseline (mussels) values

baseline <-
  subset(isotopes, sample.type == 'Mussels') %>% 
  group_by(site) %>% 
  summarize(base_deltaN = mean(deltaN),
            sd_base_deltaN = sd(deltaN))

baseline

isotopes <- left_join(isotopes,
                      baseline,
                      by = 'site')

## Average liver and muscle values for fish

isotopes2 <-
  isotopes %>% 
  group_by(ID, site, base_deltaN, sd_base_deltaN) %>% 
  summarize(deltaN = mean(deltaN),
            deltaC = mean(deltaC))

isotopes2$trophic.position <-
  with(isotopes2,
       ((deltaN - base_deltaN)/3.4)+1)


