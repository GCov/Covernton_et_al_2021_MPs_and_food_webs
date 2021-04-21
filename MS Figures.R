library(ggplot2)
library(cowplot)
library(extrafont)
loadfonts(device = "win")

## Plankton tow and jar

seawaterplot <- plot_grid(PTMPplot, PJMPplot, nrow = 1, labels = c("a", "b"),
                          label_size = 10)

tiff("Concentrations Plot.tiff",
     width = 16,
     height = 11,
     units = "cm",
     res = 700)

plot_grid(seawaterplot, speciesplot,
          nrow = 2, labels = c("", "c"), 
          rel_heights = c(1, 1.6),
          label_size = 10)

dev.off()


## Trophic level gut figure

tiff("Trophic Position MP Plot.tiff",
     width = 16,
     height = 12,
     units = "cm",
     res = 700)

plot_grid(MPTLplot, liverMPplot, ncol = 1, nrow = 2,
          labels = c("a", "b"), rel_heights = c(1, 1.3),
          align = "v", label_size = 10)

dev.off()

## Rockfish guts plot

tiff("Rockfish Guts Plot.tiff",
     width = 14,
     height = 6,
     units = "cm",
     res = 700)

plot_grid(transferplot, emptyvsfullplot, ncol = 2, nrow = 1,
          labels = c("a", "b"), axis = "bt", align = "h",
          rel_widths = c(1, 4), label_size = 10)

dev.off()


## BF and TMF plots

tiff("BF plot.tiff",
     width = 14,
     height = 12,
     units = "cm",
     res = 700)

plot_grid(BF.plot1, 
          BF.plot2,
          ncol = 1, 
          labels = c("a", "b", "c"),
          align = "v",
          label_size = 10,
          rel_heights = c(1, 1.6))

dev.off()
