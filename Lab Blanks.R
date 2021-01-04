library(plyr)
library(dplyr)
library(ggplot2)

lab_blanks <- read.csv("lab_blanks.csv", header = TRUE)

summary(lab_blanks)

lab_blanks$ID <- as.factor(lab_blanks$ID)
lab_blanks$sample.type <- as.factor(lab_blanks$sample.type)
lab_blanks$shape <- as.factor(lab_blanks$shape)
lab_blanks$colour <- as.factor(lab_blanks$colour)
lab_blanks$raman.ID <- as.factor(lab_blanks$raman.ID)

summary(lab_blanks$shape)
lab_blanks$shape <- mapvalues(lab_blanks$shape,
                              from = levels(lab_blanks$shape),
                              to = c("Fibres", "Fragments"))

summary(lab_blanks$colour)
lab_blanks$colour <- mapvalues(lab_blanks$colour,
                               from = levels(lab_blanks$colour),
                               to = c("Black",
                                      "Blue",
                                      "Clear",
                                      "Green",
                                      "Pink",
                                      "Purple",
                                      "Red"))

levels(lab_blanks$raman.ID)
lab_blanks$raman.ID <- mapvalues(lab_blanks$raman.ID,
                                 from = c("", "\n"),
                                 to = c(NA, NA))

lab_blanks$num <- 1


ggplot(lab_blanks[!is.na(lab_blanks$colour),]) +
  geom_bar(
    aes(x = 1,
        y = num,
        fill = colour),
    colour = 'black',
    size = 0.25,
    position = 'fill',
    stat = 'identity'
  ) +
  labs(x = 'Sample Type', y = 'Proportion of Particles') +
  facet_grid(shape ~ sample.type, 
             scales = 'free_x', 
             space = 'free_x',
             labeller = label_wrap_gen(width = 10)) +
  scale_fill_manual(
    values = c(
      'Grey16',
      'Steel Blue',
      'Azure2',
      'Forest Green',
      'Pink',
      'Blue Violet',
      'Red3'
    )
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(ncol = 1)) +
  theme1 +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.75, 'line'),
    panel.spacing = unit(0.25, "cm")
  )
