library(ggplot2)
library(cowplot)

## Plankton tow and jar

tiff("Seawater Plot.tiff",
     width = 14,
     height = 10,
     units = "cm",
     res = 700)

plot_grid(PTMPplot, PTAPplot, PJMPplot, PJAPplot,
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))

dev.off()


## Trophic llvel gut figure

tiff("Trophic Position Gut Plot.tiff",
     width = 19,
     height = 10,
     units = "cm",
     res = 700)

plot_grid(MPTLplot, APTLplot, ncol = 2, nrow = 1,
          labels = c("A", "B"))

dev.off()

## Livers plot

tiff("Liver Plot.tiff",
     width = 23.5,
     height = 10,
     units = "cm",
     res = 700)

plot_grid(liverMPplot, liverAPplot, nrow = 1, rel_widths = c(1, 1.4),
          align = "h")

dev.off()

