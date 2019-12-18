library(plyr)
library(dplyr)

unknown <- read.csv('unknown_particles.csv', header = TRUE)

names(unknown)

unknown <- unknown[1:10])

names(unknown)

levels(unknown$shape)

unknown$shape <- mapvalues(unknown$shape,
                           from = c('fiber', 'Fibre'),
                           to = c('fibre', 'fibre'))
levels(unknown$colour)

unknown$colour <- mapvalues(unknown$colour,
                            from = 'Blue',
                            to = 'blue')

levels(unknown$raman.ID)

unknown$raman.ID <- mapvalues(unknown$raman.ID,
                              from = c('Dye\n', 'Unknown\n'),
                              to = c('Dye', 'Unknown'))

levels(unknown$raman.ID)

unknown_summary <- 
  unknown %>% 
  group_by(shape, colour, raman.ID) %>% 
  summarize(count = length(ID))

