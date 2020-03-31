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

theme1 <-
  theme_bw() +
  theme(
    text = element_text(size = 10),
    axis.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    panel.grid = element_blank(),
    legend.title = element_blank()
  )


tiff('Isotopic_Biplot.tiff',width = 14, height = 8, units = 'cm', res = 300)

ggplot(isotopes) +  # isotopic plot
  geom_point(aes(x = deltaC,
                 y = deltaN,
                 colour = reorder(sample.type, deltaN, mean)),
             size = 0.75) +
  facet_grid(.~site) +
  scale_colour_manual(values = qualitative_hcl(palette = 'Dark 2', n = 12)) +
  theme1

dev.off()

## Estimate trophic position for each individual using mussels as a baseline

baseline <-
  subset(isotopes, sample.type == 'Mussels') %>% 
  group_by(site) %>% 
  summarize(base_deltaN = mean(deltaN))

baseline

isotopes <- left_join(isotopes,
                      baseline,
                      by = 'site')

isotopes$trophic.position <-
  with(isotopes,
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
