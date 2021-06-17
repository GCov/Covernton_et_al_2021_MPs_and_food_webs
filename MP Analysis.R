#### Load libraries ####

library(ggplot2)
library(R2jags)
library(coda)
library(lattice)
library(ggridges)
library(DHARMa)
library(reshape2)
library(plyr)
library(dplyr)
library(glmmTMB)


extract.post <- function(x){
  out <- data.frame(x$BUGSoutput$sims.list)
  long <- melt(out)
  long <- long[long$variable != "deviance" &
                 long$variable != "r", ]
  long$variable <- as.character(long$variable)
  long$variable <- as.factor(long$variable)
  long
}  # handy function for extracting posterior estimates of parameters

pal <- c("#0b3954","#bfd7ea","#ff6663","#e0ff4f","#fefffe")  # define palette

#### Plankton tow model ####

PTdata_synth <- subset(PT_data2, particle.type == "Synthetic")
PTdata_synth$particle.type <- as.character(PTdata_synth$particle.type)
PTdata_synth$particle.type <- as.factor(PTdata_synth$particle.type)
PTdata_synth$site <- as.character(PTdata_synth$site)
PTdata_synth$site <- as.factor(PTdata_synth$site) 

PTdata_synth$site <- mapvalues(PTdata_synth$site,
                               from = "Elliot Bay",
                               to = "Elliot Beach")

PTmod <- function() {
  # Likelihood
  for(i in 1:N) {
    y[i] ~ dpois(lambda_y[i])
    lambda_y[i] <- lambda_true[i] + lambda_blanks[i]
    true[i] ~ dpois(lambda_true[i])
    log(lambda_true[i]) <- 
      log(volume[i]) + alpha_site[site[i]]
    
    ## Fitted values
    fitted[i] ~ dpois(lambda_y[i])
  }
  
  ## Priors
  
  for(j in 1:nsite) {
    alpha_site[j] ~ dnorm(0, 0.01)
  }
}

## Generate initial values for MCMC

PTmodinit <- function()
{
  list(
    "alpha_site" = rnorm(3)
  )
}

## Keep track of parameters

PTmodparam <- c("alpha_site")

## Specify data

PTmoddata <-
  list(
    y = PTdata_synth$count,
    N = nrow(PTdata_synth),
    lambda_blanks = PTdata_synth$blank.mean,
    volume = PTdata_synth$sample.volume,
    site = as.integer(PTdata_synth$site),
    nsite = length(unique(PTdata_synth$site))
  )

## Run the model
PTmodrun1 <- jags.parallel(
  data = PTmoddata,
  inits = PTmodinit,
  parameters.to.save = PTmodparam,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 6193,
  model = PTmod
)

PTmodrun1
PTmodrun1mcmc <- as.mcmc(PTmodrun1)
xyplot(PTmodrun1mcmc, layout = c(6, ceiling(nvar(PTmodrun1mcmc)/6)))  # trace plots

#### Diagnostics ####

PTmodparam2 <- c("fitted", "true", "lambda_y", "lambda_true")

PTmodrun2 <- jags.parallel(
  data = PTmoddata,
  inits = PTmodinit,
  parameters.to.save = PTmodparam2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 6193,
  model = PTmod
)

PTmod.response <- t(PTmodrun2$BUGSoutput$sims.list$fitted)
PTmod.observed <- PTdata_synth$count
PTmod.fitted <- apply(t(PTmodrun2$BUGSoutput$sims.list$lambda_y),
                      1,
                      median)

check.PTmod <- createDHARMa(
  simulatedResponse = PTmod.response,
  observedResponse = PTmod.observed,
  fittedPredictedResponse = PTmod.fitted,
  integerResponse = T
)

plot(check.PTmod)
testDispersion(check.PTmod)
testZeroInflation(check.PTmod)
plotResiduals(check.PTmod, PTdata_synth$site)
plotResiduals(check.PTmod, PTdata_synth$sample.volume)

#### Inference ####

PTmodrun1long <- extract.post(PTmodrun1)

PTmodrun1long$variable <- mapvalues(
  PTmodrun1long$variable,
  from = levels(PTmodrun1long$variable),
  to = c("Coles Bay",
         "Elliot Beach",
         "Victoria Harbour")
)

PTmodrun1long$order <- c(nrow(PTmodrun1long):1)

tiff(
  'MP PT Model Posteriors.tiff', 
  height = 1.5,
  width = 6,
  units = "in",
  res = 800,
  compression = "lzw"
)

ggplot(PTmodrun1long) +  # plot parameter posteriors
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[3],
    colour = pal[1],
    alpha = 0.5, 
    size = 0.25
  ) +
  coord_cartesian(xlim = c(-7, -1.8)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()

#### Predictions ####

## Extract 'true' estimate

PTdata_synth$true.est <- apply(PTmodrun2$BUGSoutput$sims.list$true, 2, mean)
PTdata_synth$true.est.upper95 <- apply(PTmodrun2$BUGSoutput$sims.list$true, 2, quantile, 
                                    probs = 0.975)
PTdata_synth$true.est.lower95 <- apply(PTmodrun2$BUGSoutput$sims.list$true, 2, quantile, 
                                    probs = 0.025)

set.seed(5126)

PT_sim <- data.frame(
  site = c(1:3)
)

for(i in 1:3){
  y_true <-
    PTmodrun1$BUGSoutput$sims.list$alpha_site[, PT_sim$site[i]]
  PT_sim$mean[i] <- exp(mean(y_true))
  PT_sim$upper25[i] <- exp(quantile(y_true, 0.625))                  
  PT_sim$lower25[i] <- exp(quantile(y_true, 0.375))
  PT_sim$upper50[i] <- exp(quantile(y_true, 0.75))
  PT_sim$lower50[i] <- exp(quantile(y_true, 0.25))
  PT_sim$upper75[i] <- exp(quantile(y_true, 0.875))
  PT_sim$lower75[i] <- exp(quantile(y_true, 0.125))
  PT_sim$upper95[i] <- exp(quantile(y_true, 0.975))
  PT_sim$lower95[i] <- exp(quantile(y_true, 0.025))
}

PT_sim$site <- as.factor(PT_sim$site)

PT_sim$site <- mapvalues(PT_sim$site,
                           from = levels(PT_sim$site),
                           to = c("Coles Bay",
                                  "Elliot Beach",
                                  "Victoria Harbour"))

set.seed(123)

PTMPplot <- 
  ggplot() +
    geom_jitter(data = PTdata_synth,
               aes(x = site,
                   y = count/sample.volume),
               size = 1, shape = 1, colour = pal[1],
               height = 0,
               width = 0.1,
               alpha = 0.5) +
    geom_linerange(data = PT_sim,
                aes(x = site,
                    ymax = upper95,
                    ymin = lower95),
                size = 1,
                colour = pal[3]) +
    geom_point(data = PT_sim,
              aes(x = site,
                  y = mean),
              size = 2,
              fill = pal[3],
              shape = 21) +
    labs(x = 'Site',
         y = expression(paste('Particles '*L^-1))) +
    scale_y_continuous(limits = c(0, 0.15),
                       expand = c(0, 0.005),
                       breaks = seq(from = 0,
                                    to = 0.2,
                                    by = 0.05)) +
    theme1


#### Plankton jar model ####

PJdata_synth <- subset(PJ_data2, particle.type == "Synthetic")
PJdata_synth$particle.type <- as.character(PJdata_synth$particle.type)
PJdata_synth$particle.type <- as.factor(PJdata_synth$particle.type)
PJdata_synth$site <- as.character(PJdata_synth$site)
PJdata_synth$site <- as.factor(PJdata_synth$site) 

PJdata_synth$site <- mapvalues(PJdata_synth$site,
                               from = "Elliot Bay",
                               to = "Elliot Beach")

PJmod <- function() {
  # Likelihood
  for(i in 1:N) {
    y[i] ~ dpois(lambda_y[i])
    lambda_y[i] <- lambda_true[i] + lambda_blanks[i]
    true[i] ~ dpois(lambda_true[i])
    log(lambda_true[i]) <- 
      alpha_site[site[i]]
    
    ## Fitted values
    fitted[i] ~ dpois(lambda_y[i])
  }
  
  ## Priors
  
  for(j in 1:nsite) {
    alpha_site[j] ~ dnorm(0, 0.01)
  }
}

## Generate initial values for MCMC

PJmodinit <- function()
{
  list(
    "alpha_site" = rnorm(3)
  )
}

## Keep track of parameters

PJmodparam <- c("alpha_site")

## Specify data

PJmoddata <-
  list(
    y = PJdata_synth$count,
    N = nrow(PJdata_synth),
    lambda_blanks = PJdata_synth$blank.mean,
    site = as.integer(PJdata_synth$site),
    nsite = length(unique(PJdata_synth$site))
  )

## Run the model
PJmodrun1 <- jags.parallel(
  data = PJmoddata,
  inits = PJmodinit,
  parameters.to.save = PJmodparam,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 6193,
  model = PJmod
)

PJmodrun1
PJmodrun1mcmc <- as.mcmc(PJmodrun1)
xyplot(PJmodrun1mcmc, layout = c(6, ceiling(nvar(PJmodrun1mcmc)/6)))

#### Diagnostics ####
PJmodparam2 <- c("fitted", "true", "lambda_y")

PJmodrun2 <- jags.parallel(
  data = PJmoddata,
  inits = PJmodinit,
  parameters.to.save = PJmodparam2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 6193,
  model = PJmod
)

PJmod.response <- t(PJmodrun2$BUGSoutput$sims.list$fitted)
PJmod.observed <- PJdata_synth$count
PJmod.fitted <- apply(t(PJmodrun2$BUGSoutput$sims.list$lambda_y),
                      1,
                      median)

check.PJmod <- createDHARMa(
  simulatedResponse = PJmod.response,
  observedResponse = PJmod.observed,
  fittedPredictedResponse = PJmod.fitted,
  integerResponse = T
)

plot(check.PJmod)
testDispersion(check.PJmod)
testZeroInflation(check.PJmod)
plotResiduals(check.PJmod, PJdata_synth$site)

#### Inference ####

PJmodrun1long <- extract.post(PJmodrun1)

PJmodrun1long$variable <- mapvalues(
  PJmodrun1long$variable,
  from = levels(PJmodrun1long$variable),
  to = c("Coles Bay",
         "Elliot Beach",
         "Victoria Harbour")
)

PJmodrun1long$order <- c(nrow(PJmodrun1long):1)

tiff(
  'MP PJ Model Posteriors.tiff', 
  height = 1.5,
  width = 6,
  units = "in",
  res = 800,
  compression = "lzw"
)

ggplot(PJmodrun1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[3],
    colour = pal[1],
    alpha = 0.5, 
    size = 0.25
  ) +
  coord_cartesian(xlim = c(-25, 3)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

## Extract 'true' estimate

PJdata_synth$true.est <- apply(PJmodrun2$BUGSoutput$sims.list$true, 2, mean)
PJdata_synth$true.est.upper95 <- apply(PJmodrun2$BUGSoutput$sims.list$true, 2, quantile, 
                                       probs = 0.975)
PJdata_synth$true.est.lower95 <- apply(PJmodrun2$BUGSoutput$sims.list$true, 2, quantile, 
                                       probs = 0.025)

set.seed(5126)

PJ_sim <- data.frame(
  site = c(1:3),
  blank.mean = sample(PJdata_synth$blank.mean,
                      3,
                      replace = FALSE)
)

for(i in 1:3){
  lambda_true <-
    PJmodrun1$BUGSoutput$sims.list$alpha_site[, PJ_sim$site[i]]
  PJ_sim$mean[i] <- exp(mean(lambda_true))
  PJ_sim$upper25[i] <- exp(quantile(lambda_true, 0.625))
  PJ_sim$lower25[i] <- exp(quantile(lambda_true, 0.375))
  PJ_sim$upper50[i] <- exp(quantile(lambda_true, 0.75))
  PJ_sim$lower50[i] <- exp(quantile(lambda_true, 0.25))
  PJ_sim$upper75[i] <- exp(quantile(lambda_true, 0.875))
  PJ_sim$lower75[i] <- exp(quantile(lambda_true, 0.125))
  PJ_sim$upper95[i] <- exp(quantile(lambda_true, 0.975))
  PJ_sim$lower95[i] <- exp(quantile(lambda_true, 0.025))
}

PJ_sim$site <- as.factor(PJ_sim$site)

PJ_sim$site <- mapvalues(PJ_sim$site,
                         from = levels(PJ_sim$site),
                         to = c("Coles Bay",
                                "Elliot Beach",
                                "Victoria Harbour"))

set.seed(123)

PJMPplot <- 
  ggplot() +
    geom_jitter(data = PJdata_synth,
                aes(x = site,
                    y = count),
                size = 1, shape = 1, colour = pal[1],
                height = 0,
                alpha = 0.5) +
    geom_linerange(data = PJ_sim,
                   aes(x = site,
                       ymax = upper95,
                       ymin = lower95),
                   size = 1,
                   colour = pal[3]) +
    geom_point(data = PJ_sim,
               aes(x = site,
                   y = mean),
               size = 2,
               fill = pal[3],
               shape = 21) +
    labs(x = 'Site',
         y = expression(paste('Particles '*L^-1))) +
    scale_y_continuous(limits = c(0, 8),
                       expand = c(0, 0.5),
                       breaks = seq(from = 0,
                                    to = 8,
                                    by = 2)) +
    theme1


#### MP Model by Individual  ####  

MPgutdata <- subset(gutdata, !is.na(trophic.position) & 
                     particle.type == 'Synthetic')
MPgutdata$particle.type <- as.character(MPgutdata$particle.type)
MPgutdata$particle.type <- as.factor(MPgutdata$particle.type)
MPgutdata$site <- as.character(MPgutdata$site)
MPgutdata$site <- as.factor(MPgutdata$site)
MPgutdata$species <- as.character(MPgutdata$species)
MPgutdata$species <- as.factor(MPgutdata$species)

MPgutdata$site <- mapvalues(MPgutdata$site,
                            from = "Elliot Bay",
                            to = "Elliot Beach")

model1 <- function() {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dpois(lambda_y[i])
    
    lambda_y[i] <- lambda_true[i] + lambda_blanks[i]
    
    true[i] ~ dpois(lambda_true[i])
    
    log(lambda_true[i]) <-
      alpha_species[species[i]] +
      beta_TP[site[i]] * TP[i] +
      gamma_site[site[i]]
    
    TP[i] <-
      ((log(nit_lim - base[site[i]]) - log(nit_lim - animal[i])) / k) + 2
    animal[i] ~ dnorm(deltaN[i], 1/0.052)  # SD from standards
    
    ## Fitted values
    fitted[i] ~ dpois(lambda_y[i])
  }
  
  ## Priors
  
  for (j in 1:nspecies) {
    alpha_species[j] ~ dnorm(0, tau_species)
  }
  
  tau_species <- inverse(pow(sigma_species, 2))
  sigma_species ~ dexp(1)
  
  for (k in 1:nsite) {
    beta_TP[k] ~ dnorm(0, 1)
    gamma_site[k] ~ dnorm(0, 1)
    base[k] ~ dgamma(shape[k], rate[k])
    shape[k] <- pow(mean_base[k], 2) / pow(sd_base[k], 2)
    rate[k] <- mean_base[k] / pow(sd_base[k], 2)
  }
}

## Generate initial values for MCMC

model1init <- function()
{
  list(
    "sigma_species" = rexp(1),
    "beta_TP" = rnorm(3),
    "gamma_site" = rnorm(3)
  )
}

## Keep track of parameters

model1param <- c("alpha_species", "beta_TP", "gamma_site")

## Specify data

beta_zero <- 5.92
beta_one <- -0.27
nit_lim <- -beta_zero/beta_one
k <- -log((beta_zero - nit_lim)/(-nit_lim))

model1data <-
  list(
    y = MPgutdata$count,
    N = nrow(MPgutdata),
    lambda_blanks = MPgutdata$blank.mean,
    species = as.integer(MPgutdata$species),
    nspecies = length(unique(MPgutdata$species)),
    site = as.integer(MPgutdata$site),
    diff = MPgutdata$deltaN - MPgutdata$base_deltaN,
    nsite = length(unique(MPgutdata$site)),
    deltaN = MPgutdata$deltaN,
    mean_base = as.numeric(with(
      MPgutdata,
      tapply(base_deltaN,
             as.integer(site),
             mean)
    )),
    sd_base = as.numeric(with(
      MPgutdata, tapply(sd_base_deltaN, as.integer(site), mean)
    )),
    nit_lim = nit_lim,
    k = k
  )

## Run the model
run1 <- jags.parallel(
  data = model1data,
  inits = model1init,
  parameters.to.save = model1param,
  n.chains = 3,
  n.cluster = 3,
  n.iter = 7000,
  n.burnin = 500,
  n.thin = 2,
  jags.seed = 3234,
  model = model1
)

run1
run1mcmc <- as.mcmc(run1)
xyplot(run1mcmc, layout = c(6, ceiling(nvar(run1mcmc)/6)))

#### Diagnostics ####
model1param2 <- c("fitted", "true", "lambda_y", "TP")

run2 <- jags.parallel(
  data = model1data,
  inits = model1init,
  parameters.to.save = model1param2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 7000,
  n.burnin = 500,
  n.thin = 2,
  jags.seed = 3234,
  model = model1
)

model1.response <- t(run2$BUGSoutput$sims.list$fitted)
model1.observed <- MPgutdata$count
model1.fitted <- apply(t(run2$BUGSoutput$sims.list$lambda_y),
                         1,
                         median)

check.model1 <- createDHARMa(simulatedResponse = model1.response,
                               observedResponse = model1.observed, 
                               fittedPredictedResponse = model1.fitted,
                               integerResponse = T)

plot(check.model1)

plotResiduals(check.model1, MPgutdata$site)
plotResiduals(check.model1, MPgutdata$species)
plotResiduals(check.model1, MPgutdata$total.body.wet.weight)
plotResiduals(check.model1, apply(run2$BUGSoutput$sims.list$TP, 2, median))
testZeroInflation(check.model1)
testDispersion(check.model1)

plot(model1.observed-model1.fitted ~ log(model1.fitted))

#### Inference ####

run1long <- extract.post(run1)

run1long$variable <- mapvalues(run1long$variable,
                               from = levels(run1long$variable),
                               to = c("Red Rock Crab",
                                      "Starry Flounder",
                                      "Pacific Littleneck Clam",
                                      "Manila Clam",
                                      "Copper Rockfish",
                                      "Black Rockfish",
                                      "Orange Sea Cucumber",
                                      "Shiner Surfperch",
                                      "Leather Star",
                                      "Graceful Rock Crab",
                                      "Dungeness Crab",
                                      "Blue Mussel",
                                      "California Sea Cucumber",
                                      "English Sole",
                                      "Trophic Position:Coles Bay",
                                      "Trophic Position:Elliot Beach",
                                      "Trophic Position:Victoria Harbour",
                                      "Coles Bay",
                                      "Elliot Beach",
                                      "Victoria Harbour"
                                      ))

run1long$order <- c(nrow(run1long):1)

tiff(
  'MP Gut Model Posteriors.tiff', 
  height = 4.5,
  width = 6,
  units = "in",
  res = 800,
  compression = "lzw"
)

ggplot(run1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[3],
    colour = pal[1],
    alpha = 0.5, 
    size = 0.25
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.25,
    colour = pal[1]
  ) +
  coord_cartesian(xlim = c(-2.5, 2)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

## Extract 'true' estimate

MPgutdata$true.est <- apply(run2$BUGSoutput$sims.list$true, 2, mean)
MPgutdata$true.est.upper95 <- apply(run2$BUGSoutput$sims.list$true, 2, quantile, 
                                    probs = 0.975)
MPgutdata$true.est.lower95 <- apply(run2$BUGSoutput$sims.list$true, 2, quantile, 
                                    probs = 0.025)
MPgutdata$TP.est <- apply(run2$BUGSoutput$sims.list$TP, 2, mean)
MPgutdata$TP.est.lower95 <- apply(run2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.025)
MPgutdata$TP.est.upper95 <- apply(run2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.975)
MPgutdata$TP.est.lower95 <- apply(run2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.025)
MPgutdata$TP.est.upper75 <- apply(run2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.875)
MPgutdata$TP.est.lower75 <- apply(run2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.125)
MPgutdata$TP.est.upper50 <- apply(run2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.75)
MPgutdata$TP.est.lower50 <- apply(run2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.25)
MPgutdata$TP.est.upper25 <- apply(run2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.625)
MPgutdata$TP.est.lower25 <- apply(run2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.375)

set.seed(5126)

MPgutsim <- 
  data.frame(
    trophic.position = seq(
      from = 1,
      to = 6,
      length.out = 2000
    ),
    site = sample(c(1:3),
                  2000,
                  replace = TRUE),
    species = sample(MPgutdata$species,
                     2000,
                     replace = TRUE))

for(i in 1:nrow(MPgutsim)) {
  lambda_true <-
    run1$BUGSoutput$sims.list$beta_TP[, MPgutsim$site[i]]*
      MPgutsim$trophic.position[i] +
      run1$BUGSoutput$sims.list$gamma_site[, MPgutsim$site[i]] +
      run1$BUGSoutput$sims.list$alpha_species
  MPgutsim$mean[i] <- exp(mean(lambda_true))
  MPgutsim$upper25[i] <- exp(+quantile(lambda_true, 0.625))
  MPgutsim$lower25[i] <- exp(+quantile(lambda_true, 0.375))
  MPgutsim$upper50[i] <- exp(+quantile(lambda_true, 0.75))
  MPgutsim$lower50[i] <- exp(+quantile(lambda_true, 0.25))
  MPgutsim$upper75[i] <- exp(+quantile(lambda_true, 0.875))
  MPgutsim$lower75[i] <- exp(+quantile(lambda_true, 0.125))
  MPgutsim$upper95[i] <- exp(+quantile(lambda_true, 0.975))
  MPgutsim$lower95[i] <- exp(+quantile(lambda_true, 0.025))
}

MPgutsim$site <- as.factor(MPgutsim$site)

MPgutsim$site <- mapvalues(MPgutsim$site,
                           from = levels(MPgutsim$site),
                           to = c("Coles Bay",
                                  "Elliot Beach",
                                  "Victoria Harbour"))

#### Plot predictions ####

MPTLplot <-
  ggplot() +
  geom_ribbon(data = MPgutsim,
              aes(x = trophic.position,
                  ymax = upper95,
                  ymin = lower95),
              alpha = 0.05,
              size = 0.5,
              fill = pal[3]) +
  geom_ribbon(data = MPgutsim,
              aes(x = trophic.position,
                  ymax = upper75,
                  ymin = lower75),
              alpha = 0.25,
              size = 0.5,
              fill = pal[3]) +
  geom_ribbon(data = MPgutsim,
              aes(x = trophic.position,
                  ymax = upper50,
                  ymin = lower50),
              alpha = 0.5,
              size = 0.5,
              fill = pal[3]) +
  geom_ribbon(data = MPgutsim,
              aes(x = trophic.position,
                  ymax = upper25,
                  ymin = lower25),
              alpha = 0.75,
              size = 0.5,
              fill = pal[3]) +
  geom_line(data = MPgutsim,
            aes(x = trophic.position,
                y = mean),
            size = 0.5,
            colour = pal[1]) +
  geom_point(data = MPgutdata,
             aes(x = TP.est,
                 y = count),
             size = 2, shape = 20, alpha = 0.5, colour = pal[1]) +
  facet_wrap(~ site) +
  labs(x = "Trophic Position",
       y = expression(paste("Particles"~ind^-1))) +
  coord_cartesian(xlim = c(1, 6)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 8),
                     expand = c(0, 0.2),
                     breaks = seq(0, 30, 2)) +
  theme1

## Trophic position uncertainty by species

MPgutsim$species <- as.factor(MPgutsim$species)

MPgutsim$species <- mapvalues(
  MPgutsim$species,
  from = levels(MPgutsim$species),
  to = c(
    "Red Rock Crab",
    "Orange Sea Cucumber",
    "Shiner Surfperch",
    "Leather Star",
    "Graceful Rock Crab",
    "Dungeness Crab",
    "Blue Mussel",
    "California Sea Cucumber",
    "English Sole",
    "Starry Flounder",
    "Littleneck Clam",
    "Manila Clam",
    "Copper Rockfish",
    "Black Rockfish"
  )
)

tiff('Trophic Position Uncertainty Plot.tiff',
     height = 6,
     width = 6,
     units = "in",
     res = 800,
     compression = "lzw")

ggplot(MPgutdata) +
  geom_pointrange(aes(x = reorder(species, TP.est, mean),
                      y = TP.est,
                      ymin = TP.est.lower95,
                      ymax = TP.est.upper95),
                  position = position_jitter(height = 0,
                                             width = 0.5),
                  size = 0.5,
                  fatten = 1.5,
                  shape = 21,
                  alpha = 0.8,
                  colour = pal[1],
                  fill = pal[3]) +
  facet_wrap(~ site, ncol = 1) +
  labs(x = 'Site',
       y = "Trophic Position") +
  theme1 +
  theme(axis.text.x = element_text(angle = 25,
                                   hjust = 1),
        panel.grid.major.x = element_line(colour = pal[1],
                                          size = 0.2,
                                          linetype = 'dashed'))

dev.off()

## MP concentration by species

set.seed(6614)

MPgutsim2 <-
  as.data.frame(MPgutdata %>% group_by(species) %>% summarize(mean.TP = mean(TP.est)))

for(i in 1:nrow(MPgutsim2)) {
  x2 <- as.numeric()
  for (j in 1:3) {
    x1 <-
      run1$BUGSoutput$sims.list$beta_TP[, j] *
        MPgutsim2$mean.TP[i] +
        run1$BUGSoutput$sims.list$gamma_site[, j] +
        run1$BUGSoutput$sims.list$alpha_species[, i]
    lambda_true <- c(x1, x2)
    x2 <- x1
  }
  MPgutsim2$mean[i] <- exp(mean(lambda_true))
  MPgutsim2$upper25[i] <- exp(quantile(lambda_true, 0.625))
  MPgutsim2$lower25[i] <- exp(quantile(lambda_true, 0.375))
  MPgutsim2$upper50[i] <- exp(quantile(lambda_true, 0.750))
  MPgutsim2$lower50[i] <- exp(quantile(lambda_true, 0.250))
  MPgutsim2$upper75[i] <- exp(quantile(lambda_true, 0.875))
  MPgutsim2$lower75[i] <- exp(quantile(lambda_true, 0.125))
  MPgutsim2$upper95[i] <- exp(quantile(lambda_true, 0.975))
  MPgutsim2$lower95[i] <- exp(quantile(lambda_true, 0.025))
}

MPgutsim2 <- MPgutsim2 %>% arrange(mean)

MPgutsim2$order <- 1:nrow(MPgutsim2)

MPgutsim2 <- MPgutsim2 %>% arrange(species)

MPgutdata$order <- mapvalues(MPgutdata$species,
                             from = levels(MPgutdata$species),
                               to = MPgutsim2$order)

MPgutsim2$order <- as.character(MPgutsim2$order)
MPgutsim2$order <- as.numeric(MPgutsim2$order)

MPgutdata$order <- as.character(MPgutdata$order)
MPgutdata$order <- as.numeric(MPgutdata$order)

MPgutdata$common.names <- mapvalues(MPgutdata$species,
                                    from = levels(MPgutdata$species),
                                    to = c("Red Rock Crab",
                                           "Orange Sea Cucumber",
                                           "Shiner Surfperch",
                                           "Leather Star",
                                           "Graceful Rock Crab",
                                           "Dungeness Crab",
                                           "Blue Mussel",
                                           "California Sea Cucumber",
                                           "English Sole",
                                           "Starry Flounder",
                                           "Littleneck Clam",
                                           "Manila Clam",
                                           "Copper Rockfish",
                                           "Black Rockfish"))

MPgutsim2$common.names <- mapvalues(MPgutsim2$species,
                                    from = levels(MPgutdata$species),
                                    to = c("Red Rock Crab",
                                           "Orange Sea Cucumber",
                                           "Shiner Surfperch",
                                           "Leather Star",
                                           "Graceful Rock Crab",
                                           "Dungeness Crab",
                                           "Blue Mussel",
                                           "California Sea Cucumber",
                                           "English Sole",
                                           "Starry Flounder",
                                           "Littleneck Clam",
                                           "Manila Clam",
                                           "Copper Rockfish",
                                           "Black Rockfish"))

speciesplot <-
  ggplot() +
    geom_jitter(
      data = MPgutdata,
      aes(
        x = reorder(common.names, order),
        y = count
      ),
      width = 0.25,
      height = 0,
      colour = pal[1],
      size = 1,
      shape = 1,
      alpha = 0.5
    ) +
    geom_linerange(
      data = MPgutsim2,
      aes(x = reorder(common.names, order),
          ymax = upper95,
          ymin = lower95),
      size = 1,
      colour = pal[3]
    ) +
    geom_point(
      data = MPgutsim2,
      aes(x = reorder(common.names, order),
          y = mean),
      size = 2,
      colour = pal[1],
      fill = pal[3],
      shape = 21
    ) +
    labs(x = "",
         y = expression(paste('Particles '*ind^-1))) +
    scale_y_continuous(
      expand = c(0, 0.5),
      limits = c(0, 8),
      breaks = seq(0, 8, )
    ) +
    theme1 +
    theme(axis.text.x = element_text(angle = 45 , hjust = 1))


#### Export species concentration data ####


species.est <-
  MPgutdata %>%
  group_by(common.names, site) %>%
  summarize(
    TP = mean(TP.est),
    TP.low = min(TP.est),
    TP.high = max(TP.est),
    ind.conc = mean(true.est),
    ind.conc.low = min(true.est),
    ind.conc.high = max(true.est),
    ww.conc = mean(true.est / tissue.wet.weight),
    ww.conc.low = min(true.est / tissue.wet.weight),
    ww.conc.high = max(true.est / tissue.wet.weight),
    dw.conc = mean(true.est / tissue.dry.weight),
    dw.conc.low = min(true.est / tissue.dry.weight),
    dw.conc.high = max(true.est / tissue.dry.weight)
  )

write.csv(species.est,
          "species.est.csv")


#### Fish liver model ####

MPliverdata <- subset(liverdata, !is.na(trophic.position) & 
                      particle.type == 'Synthetic')
MPliverdata$particle.type <- as.character(MPliverdata$particle.type)
MPliverdata$particle.type <- as.factor(MPliverdata$particle.type)
MPliverdata$site <- as.character(MPliverdata$site)
MPliverdata$site <- as.factor(MPliverdata$site)
MPliverdata$species <- as.character(MPliverdata$species)
MPliverdata$species <- as.factor(MPliverdata$species)

MPliverdata$site <- mapvalues(MPliverdata$site,
                              from = "Elliot Bay",
                              to = "Elliot Beach")

ggplot(MPliverdata) +
  geom_point(aes(x = log(total.body.wet.weight),
                 y = log(tissue.dry.weight),
                 color = species)) +
  facet_grid(. ~ site)

ggplot(MPliverdata) +
  geom_point(aes(x = log(total.body.wet.weight),
                 y = log(tissue.wet.weight),
                 color = species)) +
  facet_grid(. ~ site)

ggplot(MPliverdata) +
  geom_point(aes(x = tissue.wet.weight,
                 y = tissue.dry.weight,
                 color = species)) +
  facet_grid(. ~ site)

plot(tissue.dry.weight ~ trophic.position, data = MPliverdata)

liver.mod <- function() {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dpois(lambda_y[i])
    
    lambda_y[i] <- lambda_true[i] + lambda_blanks[i]
    
    true[i] ~ dpois(lambda_true[i])
    
    log(lambda_true[i]) <-
      log(weight[i]) +
      alpha_species[species[i]] +
      beta_TP[site[i]] * TP[i] +
      gamma_site[site[i]]
    
    TP[i] <-
      ((log(nit_lim - base[site[i]]) - log(nit_lim - animal[i])) / k) + 2
    animal[i] ~ dnorm(deltaN[i], 1/0.052)  # SD from standards
    
    ## Fitted values
    fitted[i] ~ dpois(lambda_y[i])
  }
  
  ## Priors
  
  for (j in 1:nspecies) {
    alpha_species[j] ~ dnorm(0, tau_species)
  }
  
  tau_species <- inverse(pow(sigma_species, 2))
  sigma_species ~ dexp(1)
  
  for (k in 1:nsite) {
    beta_TP[k] ~ dnorm(0, 1)
    gamma_site[k] ~ dnorm(0, 1)
    base[k] ~ dgamma(shape[k], rate[k])
    shape[k] <- pow(mean_base[k], 2) / pow(sd_base[k], 2)
    rate[k] <- mean_base[k] / pow(sd_base[k], 2)
  }
}

## Generate initial values for MCMC

liver.mod.init <- function()
{
  list(
    "sigma_species" = rexp(1),
    "beta_TP" = rnorm(3),
    "gamma_site" = rnorm(3)
  )
}

## Keep track of parameters

liver.mod.params <- c("alpha_species", "beta_TP", "gamma_site")

## Specify data

liver.mod.data <-
  list(
    y = MPliverdata$count,
    N = nrow(MPliverdata),
    lambda_blanks = MPliverdata$blank.mean,
    weight = MPliverdata$tissue.wet.weight,
    site = as.integer(MPliverdata$site),
    nsite = length(unique(MPliverdata$site)),
    nspecies = length(unique(MPliverdata$species)),
    species = as.integer(MPliverdata$species),
    deltaN = MPliverdata$deltaN,
    mean_base = as.numeric(with(
      MPliverdata,
      tapply(base_deltaN,
             as.integer(site),
             mean)
    )),
    sd_base = as.numeric(with(
      MPliverdata, tapply(sd_base_deltaN, as.integer(site), mean)
    )),
    nit_lim = nit_lim,
    k = k
  )

## Run the model
liver.mod.run1 <- jags.parallel(
  data = liver.mod.data,
  inits = liver.mod.init,
  parameters.to.save = liver.mod.params,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 5000,
  n.burnin = 1,
  n.thin = 1,
  jags.seed = 3149,
  model = liver.mod
)

liver.mod.run1
liver.mod.run1mcmc <- as.mcmc(liver.mod.run1)
xyplot(liver.mod.run1mcmc, layout = c(6, ceiling(nvar(liver.mod.run1mcmc)/6)))

#### Diagnostics ####
liver.mod.params2 <- c("fitted", "true", "lambda_y", "TP")

liver.mod.run2 <- jags.parallel(
  data = liver.mod.data,
  inits = liver.mod.init,
  parameters.to.save = liver.mod.params2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 5000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3149,
  model = liver.mod
)

liver.mod.response <- t(liver.mod.run2$BUGSoutput$sims.list$fitted)
liver.mod.observed <- MPliverdata$count
liver.mod.fitted <- apply(t(liver.mod.run2$BUGSoutput$sims.list$lambda_y),
                           1,
                           median)

check.liver.mod <-
  createDHARMa(
    simulatedResponse = liver.mod.response,
    observedResponse = liver.mod.observed,
    fittedPredictedResponse = liver.mod.fitted,
    integerResponse = T
  )

plot(check.liver.mod)

plotResiduals(check.liver.mod, 
              apply(liver.mod.run2$BUGSoutput$sims.list$TP, 2, median))
plotResiduals(check.liver.mod,
              MPliverdata$tissue.dry.weight)
plotResiduals(check.liver.mod,
              MPliverdata$tissue.wet.weight)
plotResiduals(check.liver.mod,
              MPliverdata$species)
plotResiduals(check.liver.mod,
              MPliverdata$site)
plotResiduals(check.liver.mod,
              MPliverdata$total.body.wet.weight)

testDispersion(check.liver.mod)
testZeroInflation(check.liver.mod)

#### Inference ####

liver.mod.run1long <- extract.post(liver.mod.run1)

liver.mod.run1long$variable <-
  mapvalues(
    liver.mod.run1long$variable,
    from = levels(liver.mod.run1long$variable),
    to = c(
      "Shiner Surfperch",
      "English Sole",
      "Starry Flounder",
      "Copper Rockfish",
      "Black Rockfish",
      "Trophic Position:Coles Bay",
      "Trophic Position:Elliot Beach",
      "Trophic Position:Victoria Harbour",
      "Coles Bay",
      "Elliot Beach",
      "Victoria Harbour"
    )
  )

liver.mod.run1long$order <- c(nrow(liver.mod.run1long):1)

tiff(
  'MP Liver Model Posteriors.tiff', 
  height = 3.5,
  width = 6,
  units = "in",
  res = 800,
  compression = "lzw"
)

ggplot(liver.mod.run1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[3],
    colour = pal[1],
    alpha = 0.5, 
    size = 0.25
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.25,
    colour = pal[1]
  ) +
    coord_cartesian(xlim = c(-20, 10)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

## Extract 'true' estimate

MPliverdata$true.est <-
  apply(liver.mod.run2$BUGSoutput$sims.list$true, 2, mean)
MPliverdata$true.est.upper95 <-
  apply(liver.mod.run2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.975)
MPliverdata$true.est.lower95 <-
  apply(liver.mod.run2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.025)
MPliverdata$TP.est <-
  apply(liver.mod.run2$BUGSoutput$sims.list$TP, 2, mean)
MPliverdata$TP.est.lower95 <-
  apply(liver.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.025)
MPliverdata$TP.est.upper95 <-
  apply(liver.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.975)
MPliverdata$TP.est.lower95 <-
  apply(liver.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.025)
MPliverdata$TP.est.upper75 <-
  apply(liver.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.875)
MPliverdata$TP.est.lower75 <-
  apply(liver.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.125)
MPliverdata$TP.est.upper50 <-
  apply(liver.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.75)
MPliverdata$TP.est.lower50 <-
  apply(liver.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.25)
MPliverdata$TP.est.upper25 <-
  apply(liver.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.625)
MPliverdata$TP.est.lower25 <-
  apply(liver.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.375)

TP.mod.liver <- 
  lm(log(tissue.dry.weight) ~ trophic.position + species, 
     data = MPliverdata)

plot(resid(TP.mod.liver, type = "pearson") ~ fitted(TP.mod.liver))
summary(TP.mod.liver)

set.seed(5126)

MPliversim <- 
  data.frame(
  trophic.position = seq(
    from = 2,
    to = 4.5,
    length.out = 2000
  ),
  site = sample(c(1:3),
                2000,
                replace = TRUE),
  species = sample(MPliverdata$species,
                   2000,
                   replace = TRUE)
)

MPliversim$weight <-
  exp(rnorm(predict(TP.mod.liver, newdata = MPliversim), 0.6266))

MPliversim$species <- as.integer(MPliversim$species)

for(i in 1:2000){
  true.mean <-
    liver.mod.run1$BUGSoutput$sims.list$alpha_species +
      liver.mod.run1$BUGSoutput$sims.list$beta_TP[, MPliversim$site[i]]*
      MPliversim$trophic.position[i] +
      liver.mod.run1$BUGSoutput$sims.list$gamma_site[, MPliversim$site[i]]
  MPliversim$mean[i] <- exp(mean(true.mean))
  MPliversim$upper25[i] <- exp(quantile(true.mean, 0.625))
  MPliversim$lower25[i] <- exp(quantile(true.mean, 0.375))
  MPliversim$upper50[i] <- exp(quantile(true.mean, 0.750))
  MPliversim$lower50[i] <- exp(quantile(true.mean, 0.250))
  MPliversim$upper75[i] <- exp(quantile(true.mean, 0.875))
  MPliversim$lower75[i] <- exp(quantile(true.mean, 0.125))
  MPliversim$upper95[i] <- exp(quantile(true.mean, 0.975))
  MPliversim$lower95[i] <- exp(quantile(true.mean, 0.025))
}


MPliversim$site <- as.factor(MPliversim$site)

MPliversim$site <- mapvalues(
  MPliversim$site,
  from = levels(MPliversim$site),
  to = c("Coles Bay",
         "Elliot Beach",
         "Victoria Harbour")
)

MPliverdata$common.names <- mapvalues(MPliverdata$species,
                                      from = levels(MPliverdata$species),
                                      to = c("Shiner Surfperch",
                                             "English Sole",
                                             "Starry Flounder",
                                             "Copper Rockfish",
                                             "Black Rockfish"))

liverMPplot <-
  ggplot() +
  geom_ribbon(
    data = MPliversim,
    aes(x = trophic.position,
        ymax = upper95,
        ymin = lower95),
    alpha = 0.05,
    size = 0.5,
    fill = pal[3]
  ) +
  geom_ribbon(
    data = MPliversim,
    aes(x = trophic.position,
        ymax = upper75,
        ymin = lower75),
    alpha = 0.25,
    size = 0.5,
    fill = pal[3]
  ) +
  geom_ribbon(
    data = MPliversim,
    aes(x = trophic.position,
        ymax = upper50,
        ymin = lower50),
    alpha = 0.5,
    size = 0.5,
    fill = pal[3]
  ) +
  geom_ribbon(
    data = MPliversim,
    aes(x = trophic.position,
        ymax = upper25,
        ymin = lower25),
    alpha = 0.75,
    size = 0.5,
    fill = pal[3]
  ) +
  geom_line(
    data = MPliversim,
    aes(x = trophic.position,
        y = mean),
    size = 0.5,
    colour = pal[1]
  ) +
  geom_point(
    data = MPliverdata,
    aes(
      x = TP.est,
      y = count / tissue.dry.weight,
      fill = common.names
    ),
    colour = "black",
    size = 2,
    shape = 21,
    alpha = 0.5
  ) +
  scale_fill_manual(values = pal[1:5]) +
  facet_wrap( ~ site) +
  labs(x = 'Trophic Position',
       y = expression(paste('Particles g dry tissue ' * weight ^ -1))) +
  coord_cartesian(xlim = c(2, 4.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(
    trans = "log1p",
    expand = c(0, 0.2),
    breaks = c(0, 1, 10, 100, 400),
    limits = c(0, 400)
  ) +
  theme1 +
  guides(label.position = "bottom") +
  theme(legend.position = "bottom",
        panel.spacing = unit(0.5, "cm"),
        legend.spacing.x = unit(0.001, "cm"),
        legend.text = element_text(size = 9))


#### Rockfish ingested animals ####

#### Set up data ####

transferdata <- subset(foodweb2, 
                       sample.type == "Rockfish: Ingested Animals" &
                         particle.type == "Synthetic")

   

transferdata$site <- as.character(transferdata$site)
transferdata$site <- as.factor(transferdata$site)
transferdata$species <- as.character(transferdata$species)
transferdata$species <- as.factor(transferdata$species)
transferdata$sample.type <- as.character(transferdata$sample.type)
transferdata$sample.type <- as.factor(transferdata$sample.type)
transferdata$ID <- as.character(transferdata$ID)
transferdata$ID <- as.factor(transferdata$ID)

transfer.mod <- function() {
  for (i in 1:N) {
    y[i] ~ dpois(lambda_y[i])
    
    lambda_y[i] <- lambda_true[i] + lambda_blanks[i]
    
    true[i] ~ dpois(lambda_true[i])
    
    log(lambda_true[i]) <- 
      alpha_species[species[i]] +
      gamma_site[site[i]]
    
    ## Fitted values
    fitted[i] ~ dpois(lambda_y[i])
  }
  
  ## Priors
  
  for (j in 1:2) {
    alpha_species[j] ~ dnorm(0, 1)
  }
  
  for (k in 1:nsite) {
    gamma_site[k] ~ dnorm(0, 1)
  }
}

## Generate initial values for MCMC

transfer.mod.init <- function()
{
  list(
    "alpha_species" = rnorm(2),
    "gamma_site" = rnorm(3)
  )
}

## Keep track of parameters

transfer.mod.params <- 
  c("alpha_species", "gamma_site")

## Specify data

transfer.mod.data <-
  list(
    y = transferdata$count,
    N = nrow(transferdata),
    lambda_blanks = transferdata$blank.mean,
    weight = transferdata$tissue.dry.weight,
    site = as.integer(transferdata$site),
    nsite = length(unique(transferdata$site)),
    species = as.integer(transferdata$species)
  )

## Run the model
transfer.mod.run1 <- jags.parallel(
  data = transfer.mod.data,
  inits = transfer.mod.init,
  parameters.to.save = transfer.mod.params,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3242,
  model = transfer.mod
)

transfer.mod.run1
transfer.mod.run1mcmc <- as.mcmc(transfer.mod.run1)
xyplot(transfer.mod.run1mcmc, 
       layout = c(6, ceiling(nvar(transfer.mod.run1mcmc)/6)))

#### Diagnostics ####
transfer.mod.params2 <- c("fitted", "true", "lambda_y")

transfer.mod.run2 <- jags.parallel(
  data = transfer.mod.data,
  inits = transfer.mod.init,
  parameters.to.save = transfer.mod.params2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3242,
  model = transfer.mod
)

transfer.mod.response <- t(transfer.mod.run2$BUGSoutput$sims.list$fitted)
transfer.mod.observed <- transferdata$count
transfer.mod.fitted <- apply(t(transfer.mod.run2$BUGSoutput$sims.list$lambda_y),
                          1,
                          median)

check.transfer.mod <-
  createDHARMa(
    simulatedResponse = transfer.mod.response,
    observedResponse = transfer.mod.observed,
    fittedPredictedResponse = transfer.mod.fitted,
    integerResponse = T
  )

plot(check.transfer.mod)

plotResiduals(check.transfer.mod, transferdata$species)
plotResiduals(check.transfer.mod, transferdata$sample.type)
plotResiduals(check.transfer.mod, transferdata$ID)
plotResiduals(check.transfer.mod, transferdata$tissue.dry.weight)
plotResiduals(check.transfer.mod, transferdata$TL)
plotResiduals(check.transfer.mod, transferdata$site)

#### Inference ####

transfer.mod.run1long <- extract.post(transfer.mod.run1)

transfer.mod.run1long$variable <-
  mapvalues(
    transfer.mod.run1long$variable,
    from = levels(transfer.mod.run1long$variable),
    to = c(
      "Sebastes caurinus",
      "Sebastes melanops",
      "Coles Bay",
      "Elliot Beach",
      "Victoria Harbour"
    )
  )

transfer.mod.run1long$order <- c(nrow(transfer.mod.run1long):1)

tiff(
  'MP Transfer Model Posteriors.tiff', 
  height = 2,
  width = 6,
  units = "in",
  res = 800,
  compression = "lzw"
)

ggplot(transfer.mod.run1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[3],
    colour = pal[1],
    alpha = 0.5, 
    size = 0.25
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.25,
    colour = pal[1]
  ) +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(-3, 2.5)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

transferdata$true.est <-
  apply(transfer.mod.run2$BUGSoutput$sims.list$true, 2, mean)
transferdata$true.est.upper95 <-
  apply(transfer.mod.run2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.975)
transferdata$true.est.lower95 <-
  apply(transfer.mod.run2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.025)

set.seed(3256)

transfersim <- data.frame(data = 1)

for(i in 1:nrow(transfersim)) {
  x2 <- as.numeric()
  for (j in 1:2) {
    for (k in 1:3) {
      for (l in 1:6) {
        x1 <-
          transfer.mod.run1$BUGSoutput$sims.list$alpha_species[, j] +
            transfer.mod.run1$BUGSoutput$sims.list$gamma_site[, k]
        lambda_true <- c(x1, x2)
        x2 <- x1
      }
    }
  }
  transfersim$mean[i] <- exp(mean(lambda_true))
  transfersim$upper25[i] <- exp(quantile(lambda_true, 0.625))
  transfersim$lower25[i] <- exp(quantile(lambda_true, 0.375))
  transfersim$upper50[i] <- exp(quantile(lambda_true, 0.750))
  transfersim$lower50[i] <- exp(quantile(lambda_true, 0.250))
  transfersim$upper75[i] <- exp(quantile(lambda_true, 0.875))
  transfersim$lower75[i] <- exp(quantile(lambda_true, 0.125))
  transfersim$upper95[i] <- exp(quantile(lambda_true, 0.975))
  transfersim$lower95[i] <- exp(quantile(lambda_true, 0.025))
}

#### Plot predictions ####

transferplot <- 
  ggplot() +
    geom_linerange(
      data = transfersim,
      aes(x = 1,
          ymax = upper95,
          ymin = lower95),
      size = 0.5,
      colour = pal[3]
    ) +
    geom_point(
      data = transfersim,
      aes(x = 1,
          y = mean),
      size = 1.5,
      colour = pal[1],
      fill = pal[3],
      shape = 21
    ) +
    geom_jitter(
      data = transferdata,
      aes(
        x = 1,
        y = count
      ),
      width = 0.25,
      height = 0,
      colour = pal[1],
      size = 1,
      shape = 1,
      alpha = 0.5
    ) +
    labs(x = "",
         y = "Number of Particles") +
    scale_y_continuous(
      expand = c(0, 0.2),
      breaks = seq(0, 6, 2),
      limits = c(0, 6)
    ) +
    theme1 +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

## Calculate daily consumption of MPs

transferdata$foodmin <- 0.005 * transferdata$total.body.wet.weight
transferdata$foodmax <- 0.037 * transferdata$total.body.wet.weight

transferdata$MPmin <- transferdata$foodmin * transfersim$mean / transferdata$tissue.wet.weight
transferdata$MPmax <- transferdata$foodmax * transfersim$mean / transferdata$tissue.wet.weight

ggplot(transferdata) +
  geom_density(aes(x = MPmin),
               fill = pal[2],
               colour = pal[1],
               alpha = 0.5) +
  geom_density(aes(x = MPmax),
               fill = pal[3],
               colour = pal[1],
               alpha = 0.5) +
  labs(x = expression(paste("Number of microplastics consumed (particles "~day^-1*")")),
       y = "Density") +
  scale_x_continuous(limits = c(0, 1),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme1

#### Comparision of empty vs. full guts ####

rfishcompare <- subset(MPgutdata, sample.type == "Rockfish")

rfishcompare$full.stomach <- rfishcompare$ID %in% transferdata$ID

rfishcompare$ID <- as.character(rfishcompare$ID)
rfishcompare$ID <- as.factor(rfishcompare$ID)

rfishcompare$species <- as.character(rfishcompare$species)
rfishcompare$species <- as.factor(rfishcompare$species)

rfishcompare$full.stomach <- as.factor(rfishcompare$full.stomach)

rfish.mod <- function() {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dpois(lambda_y[i])
    
    lambda_y[i] <- lambda_true[i] + lambda_blanks[i]
    
    true[i] ~ dpois(lambda_true[i])
    
    log(lambda_true[i]) <-
      alpha_species[species[i]] +
      beta_TP * TP[i] +
      beta_TL * TL[i] +
      gamma_site[site[i]] +
      gamma_gut[full.stomach[i]]
    
    TP[i] <-
      ((log(nit_lim - base[site[i]]) - log(nit_lim - animal[i])) / k) + 2
    animal[i] ~ dnorm(deltaN[i], 1/0.052)  # SD from standards
    
    ## Fitted values
    fitted[i] ~ dpois(lambda_y[i])
  }
  
  ## Priors
  
  for (j in 1:nspecies) {
    alpha_species[j] ~ dnorm(0, tau_species)
  }
  
  tau_species <- inverse(pow(sigma_species, 2))
  sigma_species ~ dexp(1)
  
  beta_TP ~ dnorm(0, 1)
  
  beta_TL ~ dnorm(0, 1)
  
  for (k in 1:nsite) {
    gamma_site[k] ~ dnorm(0, 1)
    base[k] ~ dgamma(shape[k], rate[k])
    shape[k] <- pow(mean_base[k], 2) / pow(sd_base[k], 2)
    rate[k] <- mean_base[k] / pow(sd_base[k], 2)
  }
  
  for (l in 1:2) {
    gamma_gut[l] ~ dnorm(0, 1)
  }
}

## Generate initial values for MCMC

rfish.mod.init <- function()
{
  list(
    "sigma_species" = rexp(1),
    "beta_TP" = rnorm(1),
    "beta_TL" = rnorm(1),
    "gamma_site" = rnorm(3),
    "gamma_gut" = rnorm(2)
  )
}

## Keep track of parameters

rfish.mod.params <- 
  c("alpha_species", "beta_TP", "beta_TL", "gamma_site", "gamma_gut")

## Specify data

rfish.mod.data <-
  list(
    y = rfishcompare$count,
    N = nrow(rfishcompare),
    lambda_blanks = rfishcompare$blank.mean,
    species = as.integer(rfishcompare$species),
    nspecies = length(unique(rfishcompare$species)),
    site = as.integer(rfishcompare$site),
    diff = rfishcompare$deltaN - rfishcompare$base_deltaN,
    nsite = length(unique(rfishcompare$site)),
    deltaN = rfishcompare$deltaN,
    mean_base = as.numeric(with(
      rfishcompare,
      tapply(base_deltaN,
             as.integer(site),
             mean)
    )),
    sd_base = as.numeric(with(
      rfishcompare, tapply(sd_base_deltaN, as.integer(site), mean)
    )),
    nit_lim = nit_lim,
    k = k,
    TL = rfishcompare$TL,
    full.stomach = as.integer(rfishcompare$full.stomach)
  )

## Run the model
rfish.mod.run1 <- jags.parallel(
  data = rfish.mod.data,
  inits = rfish.mod.init,
  parameters.to.save = rfish.mod.params,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 25000,
  n.burnin = 500,
  n.thin = 5,
  jags.seed = 3242,
  model = rfish.mod
)

rfish.mod.run1
rfish.mod.run1mcmc <- as.mcmc(rfish.mod.run1)
xyplot(rfish.mod.run1mcmc, layout = c(6, ceiling(nvar(rfish.mod.run1mcmc)/6)))

#### Diagnostics ####
rfish.mod.params2 <- c("fitted", "true", "lambda_y")

rfish.mod.run2 <- jags.parallel(
  data = rfish.mod.data,
  inits = rfish.mod.init,
  parameters.to.save = rfish.mod.params2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 25000,
  n.burnin = 500,
  n.thin = 5,
  jags.seed = 3242,
  model = rfish.mod
)

rfish.mod.response <- t(rfish.mod.run2$BUGSoutput$sims.list$fitted)
rfish.mod.observed <- rfishcompare$count
rfish.mod.fitted <- apply(t(rfish.mod.run2$BUGSoutput$sims.list$lambda_y),
                             1,
                             median)

check.rfish.mod <-
  createDHARMa(
    simulatedResponse = rfish.mod.response,
    observedResponse = rfish.mod.observed,
    fittedPredictedResponse = rfish.mod.fitted,
    integerResponse = T
  )

plot(check.rfish.mod)
plotResiduals(check.rfish.mod, rfishcompare$full.stomach)
plotResiduals(check.rfish.mod, rfishcompare$species)
plotResiduals(check.rfish.mod, rfishcompare$site)
plotResiduals(check.rfish.mod, rfishcompare$tissue.wet.weight)
plotResiduals(check.rfish.mod, rfishcompare$TL)
plotResiduals(check.rfish.mod, rfishcompare$TP.est)

#### Inference ####

rfish.mod.run1long <- extract.post(rfish.mod.run1)

rfish.mod.run1long$variable <-
  mapvalues(
    rfish.mod.run1long$variable,
    from = levels(rfish.mod.run1long$variable),
    to = c(
      "Sebastes caurinus",
      "Sebastes melanops",
      "Total Length",
      "Trophic Position",
      "Empty Stomach",
      "Full Stomach",
      "Coles Bay",
      "Elliot Beach",
      "Victoria Harbour"
    )
  )

rfish.mod.run1long$order <- c(nrow(rfish.mod.run1long):1)

tiff(
  'MP Rockfish Gut Comparison Model Posteriors.tiff', 
  height = 3,
  width = 6,
  units = "in",
  res = 800,
  compression = "lzw"
)

ggplot(rfish.mod.run1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[3],
    colour = pal[1],
    alpha = 0.5, 
    size = 0.25
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.25,
    colour = pal[1]
  ) +
  coord_cartesian(xlim = c(-3, 3)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

rfishcompare$true.est <-
  apply(rfish.mod.run2$BUGSoutput$sims.list$true, 2, mean)
rfishcompare$true.est.upper95 <-
  apply(rfish.mod.run2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.975)
rfishcompare$true.est.lower95 <-
  apply(rfish.mod.run2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.025)

set.seed(3256)

rfishsim <- 
  expand.grid(species = c(1, 2),
             length = mean(rfishcompare$TL),
             full.stomach = c(1, 2))

for (i in 1:nrow(rfishsim)) {
  x2 <- as.numeric()
  for (j in 1:3) {
    for (k in 1:1000) {
      x1 <-
        rfish.mod.run1$BUGSoutput$sims.list$alpha_species[, rfishsim$species[i]] +
        rfish.mod.run1$BUGSoutput$sims.list$beta_TP *
        seq(min(rfishcompare$TP.est),
            max(rfishcompare$TP.est),
            length.out = 1000)[k] +
        rfish.mod.run1$BUGSoutput$sims.list$beta_TL *
        seq(min(rfishcompare$TL),
            max(rfishcompare$TL),
            length.out = 1000)[k] +
        rfish.mod.run1$BUGSoutput$sims.list$gamma_site[, j] +
        rfish.mod.run1$BUGSoutput$sims.list$gamma_gut[, rfishsim$full.stomach[i]]
      lambda_true <- c(x1, x2)
      x2 <- x1
    }
  }
  rfishsim$mean[i] <- exp(mean(lambda_true))
  rfishsim$upper25[i] <- exp(quantile(lambda_true, 0.625))
  rfishsim$lower25[i] <- exp(quantile(lambda_true, 0.375))
  rfishsim$upper50[i] <- exp(quantile(lambda_true, 0.75))
  rfishsim$lower50[i] <- exp(quantile(lambda_true, 0.25))
  rfishsim$upper75[i] <- exp(quantile(lambda_true, 0.875))
  rfishsim$lower75[i] <- exp(quantile(lambda_true, 0.125))
  rfishsim$upper95[i] <- exp(quantile(lambda_true, 0.975))
  rfishsim$lower95[i] <- exp(quantile(lambda_true, 0.025))
}


rfishsim$full.stomach <- as.factor(rfishsim$full.stomach)

rfishsim$full.stomach <- mapvalues(
  rfishsim$full.stomach,
  from = levels(rfishsim$full.stomach),
  to = c("Empty Stomach",
         "Full Stomach")
)

rfishcompare$full.stomach <- mapvalues(
  rfishcompare$full.stomach,
  from = levels(rfishcompare$full.stomach),
  to = c("Empty Stomach",
         "Full Stomach")
)

rfishsim$species <- as.factor(rfishsim$species)

rfishsim$species <- mapvalues(
  rfishsim$species,
  from = levels(rfishsim$species),
  to = c("Sebastes caurinus",
         "Sebastes melanops")
)

rfishsim$common.names <- 
  mapvalues(rfishsim$species,
            from = levels(rfishsim$species),
            to = c("Copper Rockfish",
                   "Black Rockfish"))

#### Plot predictions ####

emptyvsfullplot <-
  ggplot() +
    geom_jitter(
      data = rfishcompare,
      aes(
        x = full.stomach,
        y = count
      ),
      width = 0.25,
      height = 0,
      colour = pal[1],
      size = 1,
      shape = 1,
      alpha = 0.5
    ) +
    geom_linerange(
      data = rfishsim,
      aes(x = full.stomach,
          ymax = upper95,
          ymin = lower95),
      size = 0.5,
      colour = pal[3]
    ) +
    geom_point(
      data = rfishsim,
      aes(x = full.stomach,
          y = mean),
      size = 1.5,
      colour = pal[1],
      fill = pal[3],
      shape = 21
    ) +
    facet_wrap(~ common.names) +
    labs(x = "",
         y = "") +
    scale_y_continuous(
      expand = c(0, 0.2)
    ) +
    theme1


#### Export animal size data ####

write.csv(
  MPgutdata %>%
    group_by(species, site) %>%
    summarize(
      min.shell.l = min(shell.l),
      max.shell.l = max(shell.l),
      mean.shell.l = mean(shell.l),
      min.arm.length = min(arm.length),
      max.arm.length = max(arm.length),
      mean.arm.length = mean(arm.length),
      min.carapace.length = min(carapace.length),
      max.carapace.length = max(carapace.length),
      mean.carapace.length = mean(carapace.length),
      min.total.length = min(TL),
      max.total.length = max(TL),
      mean.total.length = mean(TL),
      min.body.weight = min(total.body.wet.weight),
      max.body.weight = max(total.body.wet.weight),
      mean.body.weight = mean(total.body.wet.weight),
      sample.size = length(count)
    ),
  "animalsizes.csv"
)


#### Calculate bioaccumulation factor ####

mean.water <-
  PJdata_synth %>%
  group_by(site) %>%
  summarize(water.conc = mean(true.est))

MPgutdata2 <- left_join(MPgutdata, mean.water, "site")

MPgutdata2$BF <- with(MPgutdata2,
                      ((true.est / total.body.wet.weight) * 1000) / water.conc)

MPgutdata2$feeding.strategy <-
  mapvalues(MPgutdata2$species,
            from = levels(MPgutdata2$species),
            to = c("Predator",
                   "Suspension feeder",
                   "Predator",
                   "Predator",
                   "Predator",
                   "Predator",
                   "Suspension feeder",
                   "Deposit feeder",
                   "Predator",
                   "Predator",
                   "Suspension feeder",
                   "Suspension feeder",
                   "Predator",
                   "Predator"))

mean.TP <- 
  MPgutdata2 %>%
  group_by(species) %>% 
  summarize(mean.TP = mean(TP.est))

#### Bioaccumulation plot ####

BF.plot1 <-
  ggplot(MPgutdata2) +
  geom_point(
    aes(x = TP.est,
        y = BF,
        fill = feeding.strategy),
    size = 1,
    shape = 21,
    colour = pal[1],
  ) +
  labs(x = "Trophic Position",
       y = "Bioaccumulation Factor") +
  scale_fill_manual(values = pal[2:4]) +
  scale_y_continuous(
    trans = "log1p",
    limits = c(0, 8000),
    breaks = c(0, 1, 10, 100, 1000, 5000),
    expand = c(0, 0)
  ) +
  theme1

BF.plot2 <-
  ggplot(MPgutdata2) +
  geom_boxplot(
    aes(x = reorder(common.names, TP.est, mean),
        y = BF,
        fill = feeding.strategy),
    size = 0.5,
    outlier.size = 0.5,
    colour = pal[1]
  ) +
  scale_fill_manual(values = pal[c(2,3,4)]) +
  labs(x = "Species",
       y = "Bioaccumulation Factor") +
  scale_y_continuous(
    trans = "log1p",
    limits = c(0, 8000),
    breaks = c(0, 1, 10, 100, 1000, 5000),
    expand = c(0, 0)
  ) +
  theme1 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

rfishgutdata <- subset(MPgutdata2, 
                       sample.type == "Rockfish")

BFmod1 <- glmmTMB(log(BF) ~ TP.est + (1 | species), 
                  data = rfishgutdata)

BFmod1.res <- simulateResiduals(BFmod1)

plot(BFmod1.res)

summary(BFmod1)

BFmod1.fitted <- predict(BFmod1, se.fit = TRUE, re.form = NULL)

rfishgutdata$fit <- exp(BFmod1.fitted$fit)
rfishgutdata$lower <- exp(BFmod1.fitted$fit - (1.96 * BFmod1.fitted$se.fit))
rfishgutdata$upper <- exp(BFmod1.fitted$fit + (1.96 * BFmod1.fitted$se.fit))

ggplot(rfishgutdata) +
  geom_ribbon(aes(x = TP.est,
                  ymin = lower,
                  ymax = upper,
                  fill = species),
              alpha = 0.5) +
  geom_line(aes(x = TP.est,
                y = fit,
                colour = species)) +
  geom_point(aes(x = TP.est,
                 y = BF),
             size = 1,
             colour = species) +
  labs(x = "Trophic Position",
       y = "Biaccumulation Factor") +
  scale_colour_manual(values = pal[c(2,3)]) +
  scale_fill_manual(values = pal[c(2,3)]) +
  theme1
  

#### Calculate trophic magnification factor for fish livers ####

mean(exp(liver.mod.run1$BUGSoutput$sims.list$beta_TP))

quantile(exp(liver.mod.run1$BUGSoutput$sims.list$beta_TP), probs = 0.025)
quantile(exp(liver.mod.run1$BUGSoutput$sims.list$beta_TP), probs = 0.975)
  
