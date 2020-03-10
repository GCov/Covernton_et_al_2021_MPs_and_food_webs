library(ggplot2)

isotopes <- read.csv('isotopes.csv', header = TRUE)

names(isotopes)

summary(isotopes)

isotopes$sample.type <- factor(isotopes$sample.type,
                               levels = 
                                 c('Mussels', 'Clams', 'Sea Cucumber Muscle',
                                   'Sea Star Muscle', 'Crab Muscle',
                                   'Flatfish Muscle', 'Flatfish Liver',
                                   'Surfperch Muscle', 'Surfperch Liver',
                                   'Rockfish Muscle', 'Rockfish Liver'))

ggplot(isotopes) +
  geom_point(aes(x = deltaC,
                 y = deltaN,
                 colour = sample.type),
             size = 2) +
  facet_grid(.~site) +
  theme_classic() +
  theme(legend.title = element_blank())
