library(ggplot2)
library(cowplot)

## Plankton tow and jar

seawaterplot <- plot_grid(PTMPplot, PJMPplot, nrow = 1, labels = c("A", "B"))

tiff("Concentrations Plot.tiff",
     width = 16,
     height = 11,
     units = "cm",
     res = 700)

plot_grid(seawaterplot, speciesplot,
          nrow = 2, labels = c("", "C"), 
          rel_heights = c(1, 1.6))

dev.off()


## Trophic level gut figure

tiff("Trophic Position MP Plot.tiff",
     width = 16,
     height = 12,
     units = "cm",
     res = 700)

plot_grid(MPTLplot, liverMPplot, ncol = 1, nrow = 2,
          labels = c("A", "B"), rel_heights = c(1, 1.3),
          align = "v")

dev.off()

## Rockfish guts plot

tiff("Rockfish Guts Plot.tiff",
     width = 14,
     height = 8,
     units = "cm",
     res = 700)

plot_grid(transferplot, emptyvsfullplot, ncol = 2, nrow = 1,
          labels = c("A", "B"), axis = "bt", align = "h",
          rel_widths = c(1, 4))

dev.off()

