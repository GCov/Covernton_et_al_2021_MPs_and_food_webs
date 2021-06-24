Thanks for checking out my data and analysis!!!

Required r packages:

plyr
ggplot2
dplyr
ggridges
colorspace
randomForest
R2jags*
coda
DHARMa
reshape2
glmmTMB
cowplot

*note: requires external download of JAGS software to run

There are six scripts. Each script starts by loading all necessary packages.

Script descriptions:

Field and Lab Data
Compiles sample data to be added to particle data later.

Isotopes Data
Loads, cleans, and plots the stable isotopes data.

Spectroscopy Data
Loads and cleans particle data, then uses a random forest model to predict unknown particles and plots breakdown by colour and type.

MP Analysis
Runs the Bayesian GLMMs.

MS Figures
Makes the figures for the manuscript.

Lab Blanks
...
