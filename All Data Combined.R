##### Setup #####

## load packages
library(plyr)
library(ggplot2)
library(dplyr)
library(ggridges)

## Define standard error function

se <- function(x) {
  sd(x)/sqrt(length(x))
}

#### Combine all data ####

PT_synthetic2 <- PT_synthetic
PT_synthetic2$count <- with(PT_synthetic2, count/sample.volume)

allcounts <- rbind.fill(PT_synthetic2[c(1:4, 6)],
                        PJ_synthetic,
                        MU_synthetic,
                        clams_synthetic,
                        SS_synthetic,
                        CU_synthetic,
                        CR_synthetic)

palette1 <- c('#A5A57D', '#E80D6D', '#0E6B1D', '#AC3E7E', '#91721E', '#84D6D8', 
              '#8F2C40')

png('NCAG Preliminary Counts.png', width = 20, height = 15, units = 'cm',
    pointsize = 16, res = 300)

ggplot(allcounts) +
  geom_violin(aes(x = site, 
                  y = count,
                  fill = sample.type)) +
  labs(x = '',
       y = 'MP Count per ind or L') +
  scale_fill_manual(values = palette1) +
  theme_bw() +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    text = element_text(size = 16),
    panel.spacing = unit(0.5, "lines"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()

## Ridgeline plot

ggplot(allcounts) +
  geom_density_ridges2(aes(x = count, 
                           y = sample.type,
                           fill = site),
                       alpha = 0.4,
                       panel_scaling = FALSE,
                       scale = 1.1,
                       rel_min_height = 0.01) +
  labs(x = 'MP Count per ind or L',
       y = '') +
  theme_ridges() +
  scale_fill_manual(values = palette1[5:7]) +
  scale_x_continuous(limits = c(0, 5)) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    text = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 16)
  )


## Summarize by polymer type

PT_polymer2 %>% group_by(sample.type, raman.ID) %>% summarize(count = sum(count))
PJ_polymer2 %>% group_by(sample.type, raman.ID) %>% summarize(count = sum(count))
MU_polymer2 %>% group_by(sample.type, raman.ID) %>% summarize(count = sum(count))
clams_polymer2 %>% group_by(sample.type, raman.ID) %>% summarize(count = sum(count))
SS_polymer2 %>% group_by(sample.type, raman.ID) %>% summarize(count = sum(count))
CU_polymer2 %>% group_by(sample.type, raman.ID) %>% summarize(count = sum(count))
