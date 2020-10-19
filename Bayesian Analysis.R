#### Load libraries ####

library(ggplot2)
library(R2jags)
library(coda)
library(lattice)
library(ggridges)
library(DHARMa)
library(reshape2)
library(plyr)

pal <- c("#FFC857",  # Maximum yellow red
         "#E9724C",  # Burnt sienna
         "#C5283D",  # Cardinal
         "#481D24",  # Dark sienna
         "#255F85")  # Blue sapphire

#### Plankton tow model ####

PTdata_synth <- subset(PT_data3, particle.type == "Synthetic Polymer")
PTdata_synth$particle.type <- as.character(PTdata_synth$particle.type)
PTdata_synth$particle.type <- as.factor(PTdata_synth$particle.type)
PTdata_synth$site <- as.character(PTdata_synth$site)
PTdata_synth$site <- as.factor(PTdata_synth$site) 

PTmod <- function() {
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
    alpha_site[j] ~ dnorm(0, 1)
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
    y = PTdata_synth$orig.count,
    N = nrow(PTdata_synth),
    lambda_blanks = PTdata_synth$blank.mean,
    site = as.integer(PTdata_synth$site),
    nsite = length(unique(PTdata_synth$site))
  )

## Run the model
PTmodrun1 <- jags.parallel(
  data = PTmoddata,
  inits = PTmodinit,
  parameters.to.save = PTmodparam,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 6193,
  model = PTmod
)

PTmodrun1
PTmodrun1mcmc <- as.mcmc(PTmodrun1)
xyplot(PTmodrun1mcmc, layout = c(6, ceiling(nvar(PTmodrun1mcmc)/6)))

#### Diagnostics ####
PTmodparam2 <- c("fitted", "true", "lambda_y")

PTmodrun2 <- jags.parallel(
  data = PTmoddata,
  inits = PTmodinit,
  parameters.to.save = PTmodparam2,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 6193,
  model = PTmod
)

PTmod.response <- t(PTmodrun2$BUGSoutput$sims.list$fitted)
PTmod.observed <- PTdata_synth$orig.count
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

#### Inference ####

PTmodrun1long <- extract.post(PTmodrun1)

PTmodrun1long$variable <- mapvalues(
  PTmodrun1long$variable,
  from = levels(PTmodrun1long$variable),
  to = c("Coles Bay",
         "Elliott Bay",
         "Victoria Harbour")
)

PTmodrun1long$order <- c(nrow(PTmodrun1long):1)

png(
  'MP PT Model Posteriors.png',
  width = 9,
  height = 7,
  units = 'cm',
  res = 500
)

ggplot(PTmodrun1long) +
  geom_density_ridges(
    aes(x = exp(value),
        y = reorder(variable, order, mean)),
    fill = pal[5],
    colour = pal[5],
    alpha = 0.5, 
    size = 0.25
  ) +
  coord_cartesian(xlim = c(0, 8)) +
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
  site = c(1:3),
  blank.mean = sample(PTdata_synth$blank.mean,
                      3,
                      replace = FALSE)
)

for(i in 1:3){
  lambda_true <-
    exp(
      PTmodrun1$BUGSoutput$sims.list$alpha_site[, PT_sim$site[i]]
    )
  lambda_blanks = PT_sim$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  true <- as.numeric(rpois(lambda_true, lambda_true))
  y <- as.numeric(rpois(lambda_y, lambda_y))
  PT_sim$mean[i] <- mean(true)
  PT_sim$upper25[i] <- quantile(true, 0.625)
  PT_sim$lower25[i] <- quantile(true, 0.375)
  PT_sim$upper50[i] <- quantile(true, 0.75)
  PT_sim$lower50[i] <- quantile(true, 0.25)
  PT_sim$upper75[i] <- quantile(true, 0.875)
  PT_sim$lower75[i] <- quantile(true, 0.125)
  PT_sim$upper95[i] <- quantile(true, 0.975)
  PT_sim$lower95[i] <- quantile(true, 0.025)
  PT_sim$yupper95[i] <- quantile(y, 0.975)
  PT_sim$ylower95[i] <- quantile(y, 0.025)
}

PT_sim$site <- as.factor(PT_sim$site)

PT_sim$site <- mapvalues(PT_sim$site,
                           from = levels(PT_sim$site),
                           to = c("Coles Bay",
                                  "Elliot Bay",
                                  "Victoria Harbour"))

tiff('Plankton Tows MP Bayesian Plot.tiff',
     res = 500,
     width = 9,
     height = 8,
     units = 'cm',
     pointsize = 12)

set.seed(123)

ggplot() +
  geom_linerange(data = PT_sim,
              aes(x = site,
                  ymax = upper95,
                  ymin = lower95),
              alpha = 0.05,
              size = 0.5,
              colour = pal[3]) +
  geom_linerange(data = PT_sim,
              aes(x = site,
                  ymax = upper75,
                  ymin = lower75),
              alpha = 0.25,
              size = 0.5,
              colour = pal[3]) +
  geom_linerange(data = PT_sim,
              aes(x = site,
                  ymax = upper50,
                  ymin = lower50),
              alpha = 0.5,
              size = 0.5,
              colour = pal[3]) +
  geom_linerange(data = PT_sim,
              aes(x = site,
                  ymax = upper25,
                  ymin = lower25),
              alpha = 0.75,
              size = 0.5,
              colour = pal[3]) +
  geom_point(data = PT_sim,
            aes(x = site,
                y = mean),
            size = 2,
            colour = pal[4],) +
  geom_jitter(data = PTdata_synth,
             aes(x = site,
                 y = orig.count),
             size = 0.75, shape = 1, alpha = 0.8, colour = pal[5],
             height = 0) +
  labs(x = 'Site',
       y = expression(paste('Particles '*L^-1))) +
  scale_y_continuous(limits = c(0, 10),
                     expand = c(0, 0.25),
                     breaks = seq(from = 0,
                                  to = 10,
                                  by = 2)) +
  theme1

dev.off()


#### Plankton jar model ####

PJdata_synth <- subset(PJ_data3, particle.type == "Synthetic Polymer")
PJdata_synth$particle.type <- as.character(PJdata_synth$particle.type)
PJdata_synth$particle.type <- as.factor(PJdata_synth$particle.type)
PJdata_synth$site <- as.character(PJdata_synth$site)
PJdata_synth$site <- as.factor(PJdata_synth$site) 

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
    alpha_site[j] ~ dnorm(0, 1)
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
    y = PJdata_synth$orig.count,
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
  n.cluster = 8,
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
  n.cluster = 8,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 6193,
  model = PJmod
)

PJmod.response <- t(PJmodrun2$BUGSoutput$sims.list$fitted)
PJmod.observed <- PJdata_synth$orig.count
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

#### Inference ####

PJmodrun1long <- extract.post(PJmodrun1)

PJmodrun1long$variable <- mapvalues(
  PJmodrun1long$variable,
  from = levels(PJmodrun1long$variable),
  to = c("Coles Bay",
         "Elliott Bay",
         "Victoria Harbour")
)

PJmodrun1long$order <- c(nrow(PJmodrun1long):1)

png(
  'MP PJ Model Posteriors.png',
  width = 9,
  height = 7,
  units = 'cm',
  res = 500
)

ggplot(PJmodrun1long) +
  geom_density_ridges(
    aes(x = exp(value),
        y = reorder(variable, order, mean)),
    fill = pal[5],
    colour = pal[5],
    alpha = 0.5, 
    size = 0.25
  ) +
  coord_cartesian(xlim = c(0, 3)) +
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
    exp(
      PJmodrun1$BUGSoutput$sims.list$alpha_site[, PJ_sim$site[i]]
    )
  lambda_blanks = PJ_sim$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  true <- as.numeric(rpois(lambda_true, lambda_true))
  y <- as.numeric(rpois(lambda_y, lambda_y))
  PJ_sim$mean[i] <- mean(true)
  PJ_sim$upper25[i] <- quantile(true, 0.625)
  PJ_sim$lower25[i] <- quantile(true, 0.375)
  PJ_sim$upper50[i] <- quantile(true, 0.75)
  PJ_sim$lower50[i] <- quantile(true, 0.25)
  PJ_sim$upper75[i] <- quantile(true, 0.875)
  PJ_sim$lower75[i] <- quantile(true, 0.125)
  PJ_sim$upper95[i] <- quantile(true, 0.975)
  PJ_sim$lower95[i] <- quantile(true, 0.025)
  PJ_sim$yupper95[i] <- quantile(y, 0.975)
  PJ_sim$ylower95[i] <- quantile(y, 0.025)
}

PJ_sim$site <- as.factor(PJ_sim$site)

PJ_sim$site <- mapvalues(PJ_sim$site,
                         from = levels(PJ_sim$site),
                         to = c("Coles Bay",
                                "Elliot Bay",
                                "Victoria Harbour"))

tiff('Plankton jars MP Bayesian Plot.tiff',
     res = 500,
     width = 9,
     height = 8,
     units = 'cm',
     pointsize = 12)

set.seed(123)

ggplot() +
  geom_linerange(data = PJ_sim,
                 aes(x = site,
                     ymax = upper95,
                     ymin = lower95),
                 alpha = 0.05,
                 size = 0.5,
                 colour = pal[3]) +
  geom_linerange(data = PJ_sim,
                 aes(x = site,
                     ymax = upper75,
                     ymin = lower75),
                 alpha = 0.25,
                 size = 0.5,
                 colour = pal[3]) +
  geom_linerange(data = PJ_sim,
                 aes(x = site,
                     ymax = upper50,
                     ymin = lower50),
                 alpha = 0.5,
                 size = 0.5,
                 colour = pal[3]) +
  geom_linerange(data = PJ_sim,
                 aes(x = site,
                     ymax = upper25,
                     ymin = lower25),
                 alpha = 0.75,
                 size = 0.5,
                 colour = pal[3]) +
  geom_point(data = PJ_sim,
             aes(x = site,
                 y = mean),
             size = 2,
             colour = pal[4],) +
  geom_jitter(data = PJdata_synth,
              aes(x = site,
                  y = orig.count),
              size = 0.75, shape = 1, alpha = 0.8, colour = pal[5],
              height = 0) +
  labs(x = 'Site',
       y = expression(paste('Particles '*L^-1))) +
  scale_y_continuous(limits = c(0, 8),
                     expand = c(0, 0.1),
                     breaks = seq(from = 0,
                                  to = 8,
                                  by = 2)) +
  theme1

dev.off()


#### MP Model by Individual  ####  

MPgutdata <- subset(gutdata, !is.na(trophic.position) & 
                     particle.type == 'Synthetic Polymer')
MPgutdata$particle.type <- as.character(MPgutdata$particle.type)
MPgutdata$particle.type <- as.factor(MPgutdata$particle.type)
MPgutdata$site <- as.character(MPgutdata$site)
MPgutdata$site <- as.factor(MPgutdata$site)
MPgutdata$species <- as.character(MPgutdata$species)
MPgutdata$species <- as.factor(MPgutdata$species)

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
      ((log(nit_lim - base[site[i]]) - log(nit_lim - deltaN[i])) / k) + 2
    
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
    "gamma_site" = rnorm(3),
    "base" = rgamma(3, 1, 1)
  )
}

## Keep track of parameters

model1param <- c("alpha_species", "beta_TP", "gamma_site", "base")

## Specify data

model1data <-
  list(
    y = MPgutdata$orig.count,
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
    beta_zero = 5.92,
    beta_one = -0.27,
    nit_lim = -beta_zero/beta_one,
    k = -log((beta_zero - nit_lim)/(-nit_lim))
  )

## Run the model
run1 <- jags.parallel(
  data = model1data,
  inits = model1init,
  parameters.to.save = model1param,
  n.chains = 3,
  n.cluster = 3,
  n.iter = 8000,
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
  n.cluster = 8,
  n.iter = 8000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3234,
  model = model1
)

model1.response <- t(run2$BUGSoutput$sims.list$fitted)
model1.observed <- MPgutdata$orig.count
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
plotResiduals(check.model1, apply(run2$BUGSoutput$sims.list$TP, 2, median))
testZeroInflation(check.model1)
testDispersion(check.model1)

plot(model1.observed-model1.fitted ~ log(model1.fitted))

#### Inference ####

extract.post <- function(x){
  out <- data.frame(x$BUGSoutput$sims.list)
  long <- melt(out)
  long <- long[long$variable != "deviance" &
                 long$variable != "r", ]
  long$variable <- as.character(long$variable)
  long$variable <- as.factor(long$variable)
  long
}

run1long <- extract.post(run1)

run1long$variable <- mapvalues(run1long$variable,
                               from = levels(run1long$variable),
                               to = c("Cancer productus",
                                      "Platichthys stellatus",
                                      "Protothaca staminea",
                                      "Ruditapes philippinarum",
                                      "Sebastes caurinus",
                                      "Sebastes melanops",
                                      "Cucumeria miniata",
                                      "Cymatogaster aggregata",
                                      "Dermasterias imbricata",
                                      "Metacarcinus gracilis",
                                      "Metacarcinus magister",
                                      "Mytilus spp.",
                                      "Parastichopus californicus",
                                      "Parophrys vetulus",
                                      "Coles Bay Base delta15N",
                                      "Elliot Base delta15N",
                                      "Victoria Harbour Base delta15N",
                                      "Trophic Position:Coles Bay",
                                      "Trophic Position:Elliot Bay",
                                      "Trophic Position:Victoria Harbour",
                                      "Coles Bay",
                                      "Elliott Bay",
                                      "Victoria Harbour"
                                      ))

run1long$order <- c(nrow(run1long):1)

png(
  'MP Gut Model Posteriors.png',
  width = 16,
  height = 12,
  units = 'cm',
  res = 500
)

ggplot(run1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[5],
    colour = pal[5],
    alpha = 0.5, 
    size = 0.25
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.25,
    colour = pal[3]
  ) +
  coord_cartesian(xlim = c(-1.5, 13)) +
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

MPgutsim <- data.frame(
  trophic.position = seq(
    from = 1,
    to = 5,
    length.out = 2000
  ),
  site = sample(c(1:3),
                2000,
                replace = TRUE),
  species = sample(c(1:14),
                   2000,
                   replace = TRUE),
  blank.mean = sample(MPgutdata$blank.mean,
                      2000,
                      replace = TRUE)
)

for(i in 1:2000){
  lambda_true <-
    exp(
      run1$BUGSoutput$sims.list$beta_TP[, MPgutsim$site[i]]*
        MPgutsim$trophic.position[i] +
        run1$BUGSoutput$sims.list$gamma_site[, MPgutsim$site[i]]
    )
  lambda_blanks = MPgutsim$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  true <- as.numeric(rpois(lambda_true, lambda_true))
  y <- as.numeric(rpois(lambda_y, lambda_y))
  MPgutsim$mean[i] <- mean(lambda_true)
  MPgutsim$upper25[i] <- quantile(lambda_true, 0.625)
  MPgutsim$lower25[i] <- quantile(lambda_true, 0.375)
  MPgutsim$upper50[i] <- quantile(lambda_true, 0.75)
  MPgutsim$lower50[i] <- quantile(lambda_true, 0.25)
  MPgutsim$upper75[i] <- quantile(lambda_true, 0.875)
  MPgutsim$lower75[i] <- quantile(lambda_true, 0.125)
  MPgutsim$upper95[i] <- quantile(lambda_true, 0.975)
  MPgutsim$lower95[i] <- quantile(lambda_true, 0.025)
  MPgutsim$yupper95[i] <- quantile(y, 0.975)
  MPgutsim$ylower95[i] <- quantile(y, 0.025)
}

MPgutsim$site <- as.factor(MPgutsim$site)

MPgutsim$site <- mapvalues(MPgutsim$site,
                           from = levels(MPgutsim$site),
                           to = c("Coles Bay",
                                  "Elliot Bay",
                                  "Victoria Harbour"))

#### Plot predictions ####

tiff('Trophic Position MP Bayesian Plot.tiff',
     res = 500,
     width = 16,
     height = 12,
     units = 'cm',
     pointsize = 12)

ggplot() +
  geom_ribbon(data = MPgutsim,
              aes(x = trophic.position,
                  ymax = upper95,
                  ymin = lower95),
              alpha = 0.05,
              size = 0.5,
              fill = pal[1]) +
  geom_ribbon(data = MPgutsim,
              aes(x = trophic.position,
                  ymax = upper75,
                  ymin = lower75),
              alpha = 0.25,
              size = 0.5,
              fill = pal[1]) +
  geom_ribbon(data = MPgutsim,
              aes(x = trophic.position,
                  ymax = upper50,
                  ymin = lower50),
              alpha = 0.5,
              size = 0.5,
              fill = pal[1]) +
  geom_ribbon(data = MPgutsim,
              aes(x = trophic.position,
                  ymax = upper25,
                  ymin = lower25),
              alpha = 0.75,
              size = 0.5,
              fill = pal[1],
              colour = pal[1]) +
  geom_line(data = MPgutsim,
            aes(x = trophic.position,
                y = mean),
            size = 0.5,
            colour = pal[4],
            alpha = 0.3) +
  geom_point(data = MPgutdata,
               aes(x = TP.est,
                 y = orig.count),
             size = 0.75, shape = 1, alpha = 0.8, colour = pal[5]) +
  geom_point(data = MPgutdata,
             aes(x = TP.est,
                 y = true.est),
             size = 1.5, shape = 1, alpha = 0.5, colour = pal[3]) +
  facet_wrap(~ site) +
  labs(x = 'Trophic Position',
       y = expression(paste('Particles '*ind^-1))) +
  coord_cartesian(xlim = c(1, 5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 12),
                     expand = c(0, 0.1)) +
  theme1

dev.off()

## Trophic position uncertainty by species

MPgutsim$species <- as.factor(MPgutsim$species)

MPgutsim$species <- mapvalues(
  MPgutsim$species,
  from = levels(MPgutsim$species),
  to = c(
    "Cancer productus",
    "Cucumeria miniata",
    "Cymatogaster aggregata",
    "Dermasterias imbricata",
    "Metacarcinus gracilis",
    "Metacarcinus magister",
    "Mytilus spp.",
    "Parastichopus californicus",
    "Parophrys vetulus",
    "Platichthys stellatus",
    "Protothaca staminea",
    "Ruditapes philippinarum",
    "Sebastes caurinus",
    "Sebastes melanops"
  )
)

tiff('Trophic Position Uncertainty Plot.tiff',
     res = 500,
     width = 9,
     height = 12,
     units = 'cm',
     pointsize = 12)

ggplot(MPgutdata) +
  geom_pointrange(aes(x = reorder(species, TP.est, mean),
                      y = TP.est,
                      ymin = TP.est.lower95,
                      ymax = TP.est.upper95),
                  position = position_jitter(height = 0,
                                             width = 0.5),
                  size = 0.5,
                  fatten = 0.25,
                  shape = 1,
                  alpha = 0.5,
                  colour = pal[5]) +
  facet_wrap(~ site, ncol = 1) +
  labs(x = 'Site',
       y = "Trophic Position") +
  theme1 +
  theme(axis.text.x = element_text(angle = 55,
                                   hjust = 1),
        panel.grid.major.x = element_line(colour = pal[4],
                                          size = 0.2,
                                          linetype = 'dashed'))

dev.off()




#### MP Model by Weight ####

weight.mod <- function() {
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
      ((log(nit_lim - base[site[i]]) - log(nit_lim - deltaN[i])) / k) + 2
    
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

weight.mod.init <- function()
{
  list(
    "sigma_species" = rexp(1),
    "beta_TP" = rnorm(3),
    "gamma_site" = rnorm(3),
    "base" = rgamma(3, 1, 1)
  )
}

## Keep track of parameters

weight.mod.param <- c("alpha_species", "beta_TP", "gamma_site", "base")

## Specify data

weight.mod.data <-
  list(
    y = MPgutdata$orig.count,
    N = nrow(MPgutdata),
    lambda_blanks = MPgutdata$blank.mean,
    weight = MPgutdata$tissue.weight,
    species = as.integer(MPgutdata$species),
    nspecies = length(unique(MPgutdata$species)),
    site = as.integer(MPgutdata$site),
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
    beta_zero = 5.92,
    beta_one = -0.27,
    nit_lim = -beta_zero/beta_one,
    k = -log((beta_zero - nit_lim)/(-nit_lim))
  )

## Run the model
weight.mod.run1 <- jags.parallel(
  data = weight.mod.data,
  inits = weight.mod.init,
  parameters.to.save = weight.mod.params,
  n.chains = 3,
  n.cluster = 3,
  n.iter = 10000,
  n.burnin = 500,
  n.thin = 4,
  jags.seed = 3234,
  model = weight.mod
)

weight.mod.run1
weight.mod.run1mcmc <- as.mcmc(weight.mod.run1)
xyplot(weight.mod.run1mcmc, layout = c(6, ceiling(nvar(weight.mod.run1mcmc)/6)))

#### Diagnostics ####
weight.mod.params2 <- c("fitted", "true", "lambda_y", "TP")

weight.mod.run2 <- jags.parallel(
  data = weight.mod.data,
  inits = weight.mod.init,
  parameters.to.save = weight.mod.params2,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 10000,
  n.burnin = 500,
  n.thin = 4,
  jags.seed = 3234,
  model = weight.mod
)

weight.mod.response <- t(weight.mod.run2$BUGSoutput$sims.list$fitted)
weight.mod.observed <- MPgutdata$orig.count
weight.mod.fitted <- apply(t(weight.mod.run2$BUGSoutput$sims.list$lambda_y),
                       1,
                       median)

check.weight.mod <-
  createDHARMa(
    simulatedResponse = weight.mod.response,
    observedResponse = weight.mod.observed,
    fittedPredictedResponse = weight.mod.fitted,
    integerResponse = T
  )

plot(check.weight.mod)

#### Inference ####

weight.mod.run1long <- extract.post(weight.mod.run1)

weight.mod.run1long$variable <- mapvalues(weight.mod.run1long$variable,
                               from = levels(weight.mod.run1long$variable),
                               to = c("Cancer productus",
                                      "Platichthys stellatus",
                                      "Protothaca staminea",
                                      "Ruditapes philippinarum",
                                      "Sebastes caurinus",
                                      "Sebastes melanops",
                                      "Cucumeria miniata",
                                      "Cymatogaster aggregata",
                                      "Dermasterias imbricata",
                                      "Metacarcinus gracilis",
                                      "Metacarcinus magister",
                                      "Mytilus spp.",
                                      "Parastichopus californicus",
                                      "Parophrys vetulus",
                                      "Coles Bay Base delta15N",
                                      "Elliot Base delta15N",
                                      "Victoria Harbour Base delta15N",
                                      "Trophic Position:Coles Bay",
                                      "Trophic Position:Elliot Bay",
                                      "Trophic Position:Victoria Harbour",
                                      "TDF",
                                      "Coles Bay",
                                      "Elliott Bay",
                                      "Victoria Harbour"
                               ))

weight.mod.run1long$order <- c(nrow(weight.mod.run1long):1)

png(
  'MP Gut Model by Weight Posteriors.png',
  width = 9,
  height = 9,
  units = 'cm',
  res = 500
)

ggplot(weight.mod.run1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[5],
    colour = pal[5],
    alpha = 0.5, 
    size = 0.25
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.25,
    colour = pal[3]
  ) +
  coord_cartesian(xlim = c(-3, 13)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

## Extract 'true' estimate

MPgutdata$true.weight.est <- 
  apply(weight.mod.run2$BUGSoutput$sims.list$true, 2, mean)

TP.mod <- lm(log(tissue.weight) ~ trophic.position + species, 
             data = MPgutdata)
plot(resid(TP.mod, type = "pearson") ~ fitted(TP.mod))
summary(TP.mod)

set.seed(5126)

MPgutweightsim <- data.frame(
  trophic.position = seq(
    from = 0,
    to = 4,
    length.out = 2000
  ),
  site = sample(c(1:3),
                2000,
                replace = TRUE),
  species = sample(MPgutdata$species,
                   2000,
                   replace = TRUE),
  blank.mean = sample(MPgutdata$blank.mean,
                      2000,
                      replace = TRUE)
)



MPgutweightsim$weight <- exp(rnorm(predict(TP.mod, newdata = MPgutweightsim),
                                   1.383))

MPgutweightsim$species <- as.integer(MPgutweightsim$species)

MPgutweightsim$stand.trophic.position <-
  (MPgutweightsim$trophic.position - mean(MPgutdata$trophic.position)) /
  sd(MPgutdata$trophic.position - mean(MPgutdata$trophic.position))

for(i in 1:2000){
  lambda_true <-
    exp(
      log(MPgutweightsim$weight[i]) +
      weight.mod.run1$BUGSoutput$sims.list$beta_TP[, MPgutweightsim$site[i]] * 
        MPgutweightsim$stand.trophic.position[i] +
        weight.mod.run1$BUGSoutput$sims.list$gamma_site[, MPgutweightsim$site[i]]
    )
  lambda_blanks = MPgutweightsim$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  true <- as.numeric(rpois(lambda_true, lambda_true))/MPgutweightsim$weight[i]
  y <- as.numeric(rpois(lambda_y, lambda_y))
  MPgutweightsim$mean[i] <- mean(true)
  MPgutweightsim$upper25[i] <- quantile(true, 0.625)
  MPgutweightsim$lower25[i] <- quantile(true, 0.375)
  MPgutweightsim$upper50[i] <- quantile(true, 0.75)
  MPgutweightsim$lower50[i] <- quantile(true, 0.25)
  MPgutweightsim$upper75[i] <- quantile(true, 0.875)
  MPgutweightsim$lower75[i] <- quantile(true, 0.125)
  MPgutweightsim$upper95[i] <- quantile(true, 0.975)
  MPgutweightsim$lower95[i] <- quantile(true, 0.025)
  MPgutweightsim$yupper95[i] <- quantile(y, 0.975)
  MPgutweightsim$ylower95[i] <- quantile(y, 0.025)
}

MPgutweightsim$site <- as.factor(MPgutweightsim$site)

MPgutweightsim$site <- mapvalues(MPgutweightsim$site,
                           from = levels(MPgutweightsim$site),
                           to = c("Coles Bay",
                                  "Elliot Bay",
                                  "Victoria Harbour"))

tiff('Trophic Position MP by Weight Bayesian Plot.tiff',
     res = 300,
     width = 16,
     height = 12,
     units = 'cm',
     pointsize = 12)

ggplot() +
  geom_ribbon(data = MPgutweightsim,
              aes(x = trophic.position,
                  ymax = upper95,
                  ymin = lower95),
              alpha = 0.05,
              size = 0.5,
              fill = pal[3]) +
  geom_ribbon(data = MPgutweightsim,
              aes(x = trophic.position,
                  ymax = upper75,
                  ymin = lower75),
              alpha = 0.25,
              size = 0.5,
              fill = pal[3]) +
  geom_ribbon(data = MPgutweightsim,
              aes(x = trophic.position,
                  ymax = upper50,
                  ymin = lower50),
              alpha = 0.5,
              size = 0.5,
              fill = pal[3]) +
  geom_ribbon(data = MPgutweightsim,
              aes(x = trophic.position,
                  ymax = upper25,
                  ymin = lower25),
              alpha = 0.75,
              size = 0.5,
              fill = pal[3]) +
  geom_line(data = MPgutweightsim,
            aes(x = trophic.position,
                y = mean),
            linetype = 'dashed', 
            size = 0.5) +
  geom_point(data = MPgutdata,
             aes(x = trophic.position,
                 y = orig.count/tissue.weight),
             size = 0.75, shape = 1, alpha = 0.8) +
  geom_point(data = MPgutdata,
             aes(x = trophic.position,
                 y = true.weight.est/tissue.weight),
             size = 1.5, shape = 1, alpha = 0.5, colour = pal[1]) +
  facet_wrap(~ site) +
  labs(x = 'Trophic Position',
       y = expression(paste('Particles '*g^-1))) +
  scale_x_continuous(limits = c(0, 4),
                     expand = c(0, 0)) +
  scale_y_continuous(trans = 'log1p',
                     expand = c(0, 0.1)) +
  theme1

dev.off()
