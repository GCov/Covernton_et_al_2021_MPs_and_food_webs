Thanks for checking out my data and analysis!!!

The following csv files are original data used in the analyses. For particle data,
length is particle length along longest dimension,"raman.ID" represents chemical 
identification via Raman spectroscopy, and "user.id" represents user opinion as to 
whether particles were likely natural or synthetic. The blanks ran alongside certain
samples are indicated by the same number under the "blank.match" heading.

"clams.csv"  		particle data for the clams
"crabs.csv" 		particle data for the crabs
"field_data"		collection dates and some biometrics for the animals 
			(TL = fish total length, SL = fish standard length,
			parasites = whether and how many parasitic isopods were on
			perch gills)
"flatfish.csv"		particle data from the flatfish
"isotopes.csv"		isotopic data for all samples (sample.weight is in mg)	
"lab_data.csv"		biometrics data for bivalves and seastars, as weights for all 
			samples at various stages
			(shell.l, shell.w, and shell.h = shell length, width, and height
			for bivalves, arm.length = seastar diameter across the longest
			arms, tissue.wet.weight = wet weight of sammple (g),
			tissue.dry.weight = weight of sample after drying (g),
			shell.weight = wet weight of bivalve shell, 
			total.body.wet.weight = wet weight of intact animal body,
			density.sep = whether density separation was conducted
"mussels.csv" 	 	particle data for the mussels	
"plankton_jars.csv"	particle data for the seawater jar samples
"plankton_tows.csv"	particle data for the seawater plankton tows
"PT_field_data.csv"	field data for the plankton tows
			(time = duration of tow in minutes and seconds,
			distance = tow distance (m), duration = duration of tow in
			seconds, speed = tow speed (m/s), net.opening.diameter = diameter
			of plankton net opening (m), sample.volume = calculated sample
			volume (L)
"rockfish.csv"		particle data for the rockfish
"sea_cucumbers.csv"	particle data for the sea cucumbers
"sea_stars.csv"		particle data for the sea stars
"surfperch.csv"		particle data for the surfperch
...

Required r packages:

plyr (1.8.6)
ggplot2 (3.3.5)
dplyr (1.0.7)
ggridges (0.5.3)
colorspace (2.0-2)
randomForest (4.6-14)
R2jags* (0.7-1)
coda (0.19-4)
DHARMa (0.4.4)
reshape2 (1.4.4)
cowplot (1.1.1)

*note: requires external download of JAGS software to run

There are six scripts. Each script starts by loading all necessary packages for that 
script.

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
...
