#### Load libraries ####

library(ggplot2)
library(R2jags)
library(coda)
library(lattice)
library(ggridges)
library(DHARMa)
library(reshape2)
library(plyr)

extract.post <- function(x){
  out <- data.frame(x$BUGSoutput$sims.list)
  long <- melt(out)
  long <- long[long$variable != "deviance" &
                 long$variable != "r", ]
  long$variable <- as.character(long$variable)
  long$variable <- as.factor(long$variable)
  long
}

pal <- c("#0f0a0a","#f5efed","#2292a4","#bdbf09","#d96c06")

#### Plankton tow model ####

PTdata_synth <- subset(PT_data2, particle.type == "Synthetic Polymer")
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
    y = PTdata_synth$count,
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
  n.cluster = 16,
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
    fill = pal[3],
    colour = pal[1],
    alpha = 0.5, 
    size = 0.25
  ) +
  coord_cartesian(xlim = c(0, 8)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "",
       y = "Parameter") +
  theme1 +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = pal[2]))
dev.off()


#### Demonstrate Process ####

PTmodrun2fitted <- 
  melt(data.frame(PTmodrun2$BUGSoutput$sims.list$fitted))

PTmodrun2fitted$blanks <- 
  rpois(PTdata_synth$blank.mean, 
                                PTdata_synth$blank.mean)

PTmodrun2fitted$true <- 
  melt(data.frame(PTmodrun2$BUGSoutput$sims.list$true))$value

ggplot(PTmodrun2fitted) +
  geom_density_ridges(
    aes(x = true,
        y = variable,
        scale = 1),
    fill = pal[3],
    colour = pal[1],
    alpha = 0.5,
    size = 0.25
  ) +
  geom_point(
    data = PTdata_synth,
    aes(x = blank.mean,
        y = c(1:15)),
    colour = pal[5],
    size = 2
  ) +
  geom_point(
    data = PTdata_synth,
    aes(x = count,
        y = c(1:15)),
    colour = pal[1],
    size = 2
  ) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 18)) +
  labs(x = "Microplastic Particle Count",
       y = "Sample") +
  theme1 +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = pal[2])
  )

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
              colour = pal[1]) +
  geom_linerange(data = PT_sim,
              aes(x = site,
                  ymax = upper75,
                  ymin = lower75),
              alpha = 0.25,
              size = 0.5,
              colour = pal[1]) +
  geom_linerange(data = PT_sim,
              aes(x = site,
                  ymax = upper50,
                  ymin = lower50),
              alpha = 0.5,
              size = 0.5,
              colour = pal[1]) +
  geom_linerange(data = PT_sim,
              aes(x = site,
                  ymax = upper25,
                  ymin = lower25),
              alpha = 0.75,
              size = 0.5,
              colour = pal[1]) +
  geom_point(data = PT_sim,
            aes(x = site,
                y = mean),
            size = 2,
            colour = pal[1],) +
  geom_jitter(data = PTdata_synth,
             aes(x = site,
                 y = count),
             size = 0.75, shape = 1, colour = pal[5],
             height = 0) +
  labs(x = 'Site',
       y = expression(paste('Particles '*L^-1))) +
  scale_y_continuous(limits = c(0, 10),
                     expand = c(0, 0.25),
                     breaks = seq(from = 0,
                                  to = 10,
                                  by = 2)) +
  theme1 +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = pal[2]))

dev.off()


#### Plankton jar model ####

PJdata_synth <- subset(PJ_data2, particle.type == "Synthetic Polymer")
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
                 colour = pal[1]) +
  geom_linerange(data = PJ_sim,
                 aes(x = site,
                     ymax = upper75,
                     ymin = lower75),
                 alpha = 0.25,
                 size = 0.5,
                 colour = pal[1]) +
  geom_linerange(data = PJ_sim,
                 aes(x = site,
                     ymax = upper50,
                     ymin = lower50),
                 alpha = 0.5,
                 size = 0.5,
                 colour = pal[1]) +
  geom_linerange(data = PJ_sim,
                 aes(x = site,
                     ymax = upper25,
                     ymin = lower25),
                 alpha = 0.75,
                 size = 0.5,
                 colour = pal[1]) +
  geom_point(data = PJ_sim,
             aes(x = site,
                 y = mean),
             size = 2,
             colour = pal[1],) +
  geom_jitter(data = PJdata_synth,
              aes(x = site,
                  y = count),
              size = 0.75, shape = 1, colour = pal[5],
              height = 0) +
  labs(x = 'Site',
       y = expression(paste('Particles '*L^-1))) +
  scale_y_continuous(limits = c(0, 10),
                     expand = c(0, 0.1),
                     breaks = seq(from = 0,
                                  to = 10,
                                  by = 2)) +
  theme1 +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = pal[2]))

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
    "gamma_site" = rnorm(3),
    "base" = rgamma(3, 1, 1)
  )
}

## Keep track of parameters

model1param <- c("alpha_species", "beta_TP", "gamma_site", "base")

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
plotResiduals(check.model1, apply(run2$BUGSoutput$sims.list$TP, 2, median))
testZeroInflation(check.model1)
testDispersion(check.model1)

plot(model1.observed-model1.fitted ~ log(model1.fitted))

#### Inference ####

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
                                      "Elliot Bay",
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
    colour = pal[1],
    alpha = 0.5, 
    size = 0.25
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.25,
    colour = pal[3]
  ) +
  coord_cartesian(xlim = c(-2, 13)) +
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
    to = 6,
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
              fill = pal[4]) +
  geom_ribbon(data = MPgutsim,
              aes(x = trophic.position,
                  ymax = upper75,
                  ymin = lower75),
              alpha = 0.25,
              size = 0.5,
              fill = pal[4]) +
  geom_ribbon(data = MPgutsim,
              aes(x = trophic.position,
                  ymax = upper50,
                  ymin = lower50),
              alpha = 0.5,
              size = 0.5,
              fill = pal[4]) +
  geom_ribbon(data = MPgutsim,
              aes(x = trophic.position,
                  ymax = upper25,
                  ymin = lower25),
              alpha = 0.75,
              size = 0.5,
              fill = pal[4]) +
  geom_line(data = MPgutsim,
            aes(x = trophic.position,
                y = mean),
            size = 0.5,
            colour = pal[1]) +
  geom_point(data = MPgutdata,
               aes(x = TP.est,
                 y = count),
             size = 1, shape = 20, alpha = 0.8, colour = pal[1]) +
  geom_point(data = MPgutdata,
             aes(x = TP.est,
                 y = true.est),
             size = 1.5, shape = 1, alpha = 0.5, colour = pal[3]) +
  facet_wrap(~ site) +
  labs(x = 'Trophic Position',
       y = expression(paste('Particles '*ind^-1))) +
  coord_cartesian(xlim = c(1, 6)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 7),
                     expand = c(0, 0.1),
                     breaks = c(seq(0, 12, 2))) +
  theme1 +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = pal[2]))

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

## MP concentration by species

set.seed(6614)

MPgutsim2 <- expand.grid(
  trophic.position = 2.9,
  site = 1,
  species = c(1:14),
  blank.mean = mean(MPgutdata$blank.mean)
)

for(i in 1:14){
  lambda_true <-
    exp(
      run1$BUGSoutput$sims.list$beta_TP[, MPgutsim2$site[i]]*
        MPgutsim2$trophic.position[i] +
        run1$BUGSoutput$sims.list$gamma_site[, MPgutsim2$site[i]]+
        run1$BUGSoutput$sims.list$alpha_species[, MPgutsim2$species[i]]
    )
  lambda_blanks = MPgutsim2$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  true <- as.numeric(rpois(lambda_true, lambda_true))
  y <- as.numeric(rpois(lambda_y, lambda_y))
  MPgutsim2$mean[i] <- mean(lambda_true)
  MPgutsim2$upper25[i] <- quantile(lambda_true, 0.625)
  MPgutsim2$lower25[i] <- quantile(lambda_true, 0.375)
  MPgutsim2$upper50[i] <- quantile(lambda_true, 0.75)
  MPgutsim2$lower50[i] <- quantile(lambda_true, 0.25)
  MPgutsim2$upper75[i] <- quantile(lambda_true, 0.875)
  MPgutsim2$lower75[i] <- quantile(lambda_true, 0.125)
  MPgutsim2$upper95[i] <- quantile(lambda_true, 0.975)
  MPgutsim2$lower95[i] <- quantile(lambda_true, 0.025)
  MPgutsim2$yupper95[i] <- quantile(y, 0.975)
  MPgutsim2$ylower95[i] <- quantile(y, 0.025)
}

MPgutsim2$species <- as.factor(MPgutsim2$species)
MPgutsim2$site <- as.factor(MPgutsim2$site)

MPgutsim2$species <- mapvalues(
  MPgutsim2$species,
  from = levels(MPgutsim2$species),
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


MPgutsim2$site <- mapvalues(
  MPgutsim2$site,
  from = levels(MPgutsim2$site),
  to = c("Coles Bay",
         "Elliot Bay",
         "Victoria Harbour")
)

tiff('Species Plot.tiff',
     res = 500,
     width = 14,
     height = 8,
     units = 'cm',
     pointsize = 12)

ggplot() +
  geom_linerange(
    data = MPgutsim2,
    aes(x = species,
        ymax = upper95,
        ymin = lower95),
    alpha = 0.05,
    size = 0.5,
    colour = pal[5]
  ) +
  geom_linerange(
    data = MPgutsim2,
    aes(x = species,
        ymax = upper75,
        ymin = lower75),
    alpha = 0.25,
    size = 0.5,
    colour = pal[5]
  ) +
  geom_linerange(
    data = MPgutsim2,
    aes(x = species,
        ymax = upper50,
        ymin = lower50),
    alpha = 0.5,
    size = 0.5,
    colour = pal[5]
  ) +
  geom_linerange(
    data = MPgutsim2,
    aes(x = species,
        ymax = upper25,
        ymin = lower25),
    alpha = 0.75,
    size = 0.5,
    colour = pal[5]
  ) +
  geom_point(
    data = MPgutsim2,
    aes(x = species,
        y = mean),
    size = 1.5,
    colour = pal[1],
    fill = pal[2],
    shape = 21
  ) +
  geom_jitter(
    data = MPgutdata,
    aes(
      x = species,
      y = count
    ),
    width = 0.25,
    height = 0,
    colour = pal[3],
    size = 1,
    shape = 1,
    alpha = 0.5
  ) +
  labs(x = "",
       y = expression(paste('Particles '*ind^-1))) +
  scale_y_continuous(
    expand = c(0, 0.1)
  ) +
  theme1 +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = pal[2]),
        legend.background = element_rect(fill = pal[2]),
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()



#### Body size model with just fish ####

fishgutdata <- subset(MPgutdata, 
                      sample.type == "Flatfish" |
                        sample.type == "Rockfish" |
                        sample.type == "Surfperch")

fishgutdata$species <- as.character(fishgutdata$species)
fishgutdata$species <- as.factor(fishgutdata$species)

names(fishgutdata)

nrow(fishgutdata)

plot(TP.est ~ TL, data = fishgutdata)
plot(TL ~ species, data = fishgutdata)
plot(count ~ TL, data = fishgutdata)
plot(count ~ log(total.body.wet.weight), data = fishgutdata)

ggplot(fishgutdata) +
  geom_point(aes(x = log(TL),
                 y = log(total.body.wet.weight),
                 colour = species)) +
  geom_smooth(aes(x = log(TL),
                  y = log(total.body.wet.weight)),
              method = 'lm') +
  facet_grid(. ~ site)


fishmodel1 <- function() {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dpois(lambda_y[i])
    
    lambda_y[i] <- lambda_true[i] + lambda_blanks[i]
    
    true[i] ~ dpois(lambda_true[i])
    
    log(lambda_true[i]) <-
      alpha_species[species[i]] +
      beta_length[site[i]] * length[i] +
      gamma_site[site[i]]
    
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
    beta_length[k] ~ dnorm(0, 1)
    gamma_site[k] ~ dnorm(0, 1)
  }
}

## Generate initial values for MCMC

fishmodel1init <- function()
{
  list(
    "beta_length" = rnorm(3),
    "gamma_site" = rnorm(3),
    "sigma_species" = rexp(1)
  )
}

## Keep track of parameters

fishmodel1param <- c("alpha_species", "beta_length", "gamma_site")

## Specify data

fishmodel1data <-
  list(
    y = fishgutdata$count,
    N = nrow(fishgutdata),
    lambda_blanks = fishgutdata$blank.mean,
    species = as.integer(fishgutdata$species),
    nspecies = length(unique(fishgutdata$species)),
    site = as.integer(fishgutdata$site),
    nsite = length(unique(fishgutdata$site)),
    length = as.numeric(scale(fishgutdata$TL), center = TRUE)
  )

## Run the model
fishrun1 <- jags.parallel(
  data = fishmodel1data,
  inits = fishmodel1init,
  parameters.to.save = fishmodel1param,
  n.chains = 3,
  n.cluster = 3,
  n.iter = 5000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3234,
  model = fishmodel1
)

fishrun1
fishrun1mcmc <- as.mcmc(fishrun1)
xyplot(fishrun1mcmc, layout = c(6, ceiling(nvar(fishrun1mcmc)/6)))

#### Diagnostics ####
fishmodel1param2 <- c("fitted", "true", "lambda_y", "TP")

fishrun2 <- jags.parallel(
  data = fishmodel1data,
  inits = fishmodel1init,
  parameters.to.save = fishmodel1param2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 5000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3234,
  model = fishmodel1
)

fishmodel1.response <- t(fishrun2$BUGSoutput$sims.list$fitted)
fishmodel1.observed <- fishgutdata$count
fishmodel1.fitted <- apply(t(fishrun2$BUGSoutput$sims.list$lambda_y),
                       1,
                       median)

check.fishmodel1 <- createDHARMa(simulatedResponse = fishmodel1.response,
                             observedResponse = fishmodel1.observed, 
                             fittedPredictedResponse = fishmodel1.fitted,
                             integerResponse = T)

plot(check.fishmodel1)

plotResiduals(check.fishmodel1, fishgutdata$site)
plotResiduals(check.fishmodel1, fishgutdata$species)
plotResiduals(check.fishmodel1, fishgutdata$TL)
testZeroInflation(check.fishmodel1)
testDispersion(check.fishmodel1)

plot(fishmodel1.observed-fishmodel1.fitted ~ log(fishmodel1.fitted))

#### Inference ####

fishrun1long <- extract.post(fishrun1)

fishrun1long$variable <- mapvalues(fishrun1long$variable,
                               from = levels(fishrun1long$variable),
                               to = c("Cymatogaster aggregata",
                                      "Parophrys vetulus",
                                      "Platichthys stellatus",
                                      "Sebastes caurinus",
                                      "Sebastes melanops",
                                      "Total length (cm):Coles Bay",
                                      "Total length (cm):Elliot Bay",
                                      "Total length (cm):Victoria Harbour",
                                      "Coles Bay",
                                      "Elliott Bay",
                                      "Victoria Harbour"
                               ))

fishrun1long$order <- c(nrow(fishrun1long):1)

png(
  'Fish Gut Body Size Model Posteriors.png',
  width = 16,
  height = 12,
  units = 'cm',
  res = 500
)

ggplot(fishrun1long) +
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
  coord_cartesian(xlim = c(-3, 2.1)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

## Extract 'true' estimate

fishgutdata$true.est <-
  apply(fishrun2$BUGSoutput$sims.list$true, 2, mean)
fishgutdata$true.est.upper95 <-
  apply(fishrun2$BUGSoutput$sims.list$true, 2,
        quantile, probs = 0.975)
fishgutdata$true.est.lower95 <-
  apply(fishrun2$BUGSoutput$sims.list$true, 2,
        quantile,
        probs = 0.025)

set.seed(5126)

fishgutsim <- data.frame(
  length = seq(
    from = 5,
    to = 35,
    length.out = 2000
  ),
  site = sample(c(1:3),
                2000,
                replace = TRUE),
  species = sample(c(1:5),
                   2000,
                   replace = TRUE),
  blank.mean = sample(fishgutdata$blank.mean,
                      2000,
                      replace = TRUE)
)

fishgutsim$length.stand <-
  (log(fishgutsim$length) - log(mean(fishgutdata$TL))) /
  sd(log(fishgutsim$length) - log(mean(fishgutdata$TL)))

for(i in 1:2000){
  lambda_true <-
    exp(
      fishrun1$BUGSoutput$sims.list$beta_length[, fishgutsim$site[i]]*
        fishgutsim$length.stand[i] +
        fishrun1$BUGSoutput$sims.list$gamma_site[, fishgutsim$site[i]]
    )
  lambda_blanks = fishgutsim$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  true <- as.numeric(rpois(lambda_true, lambda_true))
  y <- as.numeric(rpois(lambda_y, lambda_y))
  fishgutsim$mean[i] <- mean(lambda_true)
  fishgutsim$upper25[i] <- quantile(lambda_true, 0.625)
  fishgutsim$lower25[i] <- quantile(lambda_true, 0.375)
  fishgutsim$upper50[i] <- quantile(lambda_true, 0.75)
  fishgutsim$lower50[i] <- quantile(lambda_true, 0.25)
  fishgutsim$upper75[i] <- quantile(lambda_true, 0.875)
  fishgutsim$lower75[i] <- quantile(lambda_true, 0.125)
  fishgutsim$upper95[i] <- quantile(lambda_true, 0.975)
  fishgutsim$lower95[i] <- quantile(lambda_true, 0.025)
  fishgutsim$yupper95[i] <- quantile(y, 0.975)
  fishgutsim$ylower95[i] <- quantile(y, 0.025)
}

fishgutsim$site <- as.factor(fishgutsim$site)

fishgutsim$site <- mapvalues(fishgutsim$site,
                           from = levels(fishgutsim$site),
                           to = c("Coles Bay",
                                  "Elliot Bay",
                                  "Victoria Harbour"))

#### Plot predictions ####

tiff('Body Size Fish Bayesian Plot.tiff',
     res = 500,
     width = 16,
     height = 12,
     units = 'cm',
     pointsize = 12)

ggplot() +
  geom_ribbon(data = fishgutsim,
              aes(x = length,
                  ymax = upper95,
                  ymin = lower95),
              alpha = 0.05,
              size = 0.5,
              fill = pal[1]) +
  geom_ribbon(data = fishgutsim,
              aes(x = length,
                  ymax = upper75,
                  ymin = lower75),
              alpha = 0.25,
              size = 0.5,
              fill = pal[1]) +
  geom_ribbon(data = fishgutsim,
              aes(x = length,
                  ymax = upper50,
                  ymin = lower50),
              alpha = 0.5,
              size = 0.5,
              fill = pal[1]) +
  geom_ribbon(data = fishgutsim,
              aes(x = length,
                  ymax = upper25,
                  ymin = lower25),
              alpha = 0.75,
              size = 0.5,
              fill = pal[1],
              colour = pal[1]) +
  geom_line(data = fishgutsim,
            aes(x = length,
                y = mean),
            size = 0.5,
            colour = pal[4],
            alpha = 0.3) +
  geom_point(data = fishgutdata,
             aes(x = TL,
                 y = count),
             size = 0.75, shape = 1, alpha = 0.8, colour = pal[5]) +
  geom_linerange(data = fishgutdata,
                 aes(x = TL,
                     ymin = true.est.lower95,
                     ymax = true.est.upper95),
                     size = 0.5, alpha = 0.5, colour = pal[3]) +
  geom_point(data = fishgutdata,
             aes(x = TL,
                 y = true.est),
             size = 1.5, shape = 1, alpha = 0.75, colour = pal[3]) +
  facet_wrap(~ site) +
  scale_x_continuous(trans = 'log1p',
                     breaks = c(5, 10, 20, 35),
                     expand = c(0, 0)) +
  labs(x = 'Total Length (cm)',
       y = expression(paste('Particles '*ind^-1))) +
  theme1 +
  theme(panel.spacing = unit(0.5, "cm"))

dev.off()


#### Fish liver model ####

MPliverdata <- subset(liverdata, !is.na(trophic.position) & 
                      particle.type == 'Synthetic Polymer')
MPliverdata$particle.type <- as.character(MPliverdata$particle.type)
MPliverdata$particle.type <- as.factor(MPliverdata$particle.type)
MPliverdata$site <- as.character(MPliverdata$site)
MPliverdata$site <- as.factor(MPliverdata$site)
MPliverdata$species <- as.character(MPliverdata$species)
MPliverdata$species <- as.factor(MPliverdata$species)

ggplot(MPliverdata) +
  geom_point(aes(x = total.body.wet.weight,
                 y = tissue.dry.weight,
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
    "gamma_site" = rnorm(3),
    "base" = rgamma(3, 1, 1)
  )
}

## Keep track of parameters

liver.mod.params <- c("alpha_species", "beta_TP", "gamma_site", "base")

## Specify data

liver.mod.data <-
  list(
    y = MPliverdata$count,
    N = nrow(MPliverdata),
    lambda_blanks = MPliverdata$blank.mean,
    weight = MPliverdata$tissue.dry.weight,
    species = as.integer(MPliverdata$species),
    nspecies = length(unique(MPliverdata$species)),
    site = as.integer(MPliverdata$site),
    nsite = length(unique(MPliverdata$site)),
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
  n.burnin = 500,
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

testDispersion(check.liver.mod)
testZeroInflation(check.liver.mod)

#### Inference ####

liver.mod.run1long <- extract.post(liver.mod.run1)

liver.mod.run1long$variable <-
  mapvalues(
    liver.mod.run1long$variable,
    from = levels(liver.mod.run1long$variable),
    to = c(
      "Cymatogaster aggregata",
      "Parophrys vetulus",
      "Platichthys stellatus",
      "Sebastes caurinus",
      "Sebastes melanops",
      "Coles Bay Base delta15N",
      "Elliot Base delta15N",
      "Victoria Harbour Base delta15N",
      "Trophic Position:Coles Bay",
      "Trophic Position:Elliot Bay",
      "Trophic Position:Victoria Harbour",
      "Coles Bay",
      "Elliott Bay",
      "Victoria Harbour"
    )
  )

liver.mod.run1long$order <- c(nrow(liver.mod.run1long):1)

png(
  'MP Liver Model Posteriors.png',
  width = 9,
  height = 9,
  units = 'cm',
  res = 500
)

ggplot(liver.mod.run1long) +
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
  coord_cartesian(xlim = c(-6, 14)) +
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

MPliversim <- data.frame(
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
                   replace = TRUE),
  blank.mean = sample(MPliverdata$blank.mean,
                      2000,
                      replace = TRUE)
)

MPliversim$weight <-
  exp(rnorm(predict(TP.mod.liver, newdata = MPliversim), 0.6266))

MPliversim$species <- as.integer(MPliversim$species)

for(i in 1:2000){
  true.mean <-
    exp(
      liver.mod.run1$BUGSoutput$sims.list$beta_TP[, MPliversim$site[i]]*
        MPliversim$trophic.position[i] +
        liver.mod.run1$BUGSoutput$sims.list$gamma_site[, MPliversim$site[i]]
    )
  lambda_true <- true*MPliversim$weight[i]
  lambda_blanks <- MPliversim$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  true.conc <- as.numeric(rpois(lambda_true, lambda_true))/MPliversim$weight[i]
  y <- as.numeric(rpois(lambda_y, lambda_y))
  MPliversim$median[i] <- median(true.mean)
  MPliversim$upper25[i] <- quantile(true.mean, 0.625)
  MPliversim$lower25[i] <- quantile(true.mean, 0.375)
  MPliversim$upper50[i] <- quantile(true.mean, 0.750)
  MPliversim$lower50[i] <- quantile(true.mean, 0.250)
  MPliversim$upper75[i] <- quantile(true.mean, 0.875)
  MPliversim$lower75[i] <- quantile(true.mean, 0.125)
  MPliversim$upper95[i] <- quantile(true.mean, 0.975)
  MPliversim$lower95[i] <- quantile(true.mean, 0.025)
  MPliversim$yupper95[i] <- quantile(y, 0.975)
  MPliversim$ylower95[i] <- quantile(y, 0.025)
}


MPliversim$site <- as.factor(MPliversim$site)

MPliversim$site <- mapvalues(MPliversim$site,
                           from = levels(MPliversim$site),
                           to = c("Coles Bay",
                                  "Elliot Bay",
                                  "Victoria Harbour"))

#### Plot predictions ####

tiff('Trophic Position MP Liver Bayesian Plot.tiff',
     res = 500,
     width = 16,
     height = 12,
     units = 'cm',
     pointsize = 12)

ggplot() +
  geom_ribbon(
    data = MPliversim,
    aes(x = trophic.position,
        ymax = upper95,
        ymin = lower95),
    alpha = 0.05,
    size = 0.5,
    fill = pal[5]
  ) +
  geom_ribbon(
    data = MPliversim,
    aes(x = trophic.position,
        ymax = upper75,
        ymin = lower75),
    alpha = 0.25,
    size = 0.5,
    fill = pal[5]
  ) +
  geom_ribbon(
    data = MPliversim,
    aes(x = trophic.position,
        ymax = upper50,
        ymin = lower50),
    alpha = 0.5,
    size = 0.5,
    fill = pal[5]
  ) +
  geom_ribbon(
    data = MPliversim,
    aes(x = trophic.position,
        ymax = upper25,
        ymin = lower25),
    alpha = 0.75,
    size = 0.5,
    fill = pal[5]
  ) +
  geom_line(
    data = MPliversim,
    aes(x = trophic.position,
        y = median),
    size = 0.5,
    colour = pal[1],
    alpha = 0.3
  ) +
  geom_point(
    data = MPliverdata,
    aes(
      x = TP.est,
      y = true.est / tissue.wet.weight,
      fill = species
    ),
    colour = "black",
    size = 2,
    shape = 21,
    alpha = 0.5
  ) +
  scale_fill_manual(values = pal[1:5]) +
  facet_wrap(~ site) +
  labs(x = 'Trophic Position',
       y = expression(paste('Particles g dry tissue ' * weight ^ -1))) +
  coord_cartesian(xlim = c(2, 4.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(
    trans = 'log1p',
    expand = c(0, 0.01),
    breaks = c(0, 1, 10, 40)
  ) +
  theme1 +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = pal[2]),
        legend.background = element_rect(fill = pal[2]))

dev.off()


#### Rockfish ingested animals vs. gut ####

transferdata <- subset(foodweb2, 
                       sample.type == "Rockfish: Ingested Animals" &
                         particle.type == "Synthetic Polymer" |
                         sample.type == "Rockfish" &
                         particle.type == "Synthetic Polymer")

transferdata$site <- as.character(transferdata$site)
transferdata$site <- as.factor(transferdata$site)
transferdata$species <- as.character(transferdata$species)
transferdata$species <- as.factor(transferdata$species)
transferdata$sample.type <- as.character(transferdata$sample.type)
transferdata$sample.type <- as.factor(transferdata$sample.type)
transferdata$sample.type <- mapvalues(transferdata$sample.type,
                                      from = levels(transferdata$sample.type),
                                      to = c("Gut",
                                             "Gut Animals"))

transfer.mod <- function() {
  for (i in 1:N) {
    y[i] ~ dpois(lambda_y[i])
    
    lambda_y[i] <- lambda_true[i] + lambda_blanks[i]
    
    true[i] ~ dpois(lambda_true[i])
    
    log(lambda_true[i]) <- alpha_gut[sample.type[i]]
    
    ## Fitted values
    fitted[i] ~ dpois(lambda_y[i])
  }
  
  ## Priors
  for (j in 1:2) {
    alpha_gut[j] ~ dnorm(0, 1)
  }
}

## Generate initial values for MCMC

transfer.mod.init <- function()
{
  list(
    "alpha_gut" = rnorm(2)
  )
}

## Keep track of parameters

transfer.mod.params <- c("alpha_gut")

## Specify data

transfer.mod.data <-
  list(
    y = transferdata$count,
    N = nrow(transferdata),
    lambda_blanks = transferdata$blank.mean,
    sample.type = as.integer(transferdata$sample.type)
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
xyplot(transfer.mod.run1mcmc, layout = c(6, ceiling(nvar(transfer.mod.run1mcmc)/6)))

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

#### Inference ####

transfer.mod.run1long <- extract.post(transfer.mod.run1)

transfer.mod.run1long$variable <-
  mapvalues(
    transfer.mod.run1long$variable,
    from = levels(transfer.mod.run1long$variable),
    to = c(
      "Gut",
      "Gut Animals"
    )
  )

transfer.mod.run1long$order <- c(nrow(transfer.mod.run1long):1)

png(
  'MP Transfer Model Posteriors.png',
  width = 9,
  height = 9,
  units = 'cm',
  res = 500
)

ggplot(transfer.mod.run1long) +
  geom_density_ridges(
    aes(x = exp(value),
        y = reorder(variable, order, mean)),
    fill = pal[5],
    colour = pal[5],
    alpha = 0.5, 
    size = 0.25
  ) +
  coord_cartesian(xlim = c(0, 1.5)) +
  scale_x_continuous(expand = c(0,0)) +
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

transfersim <- data.frame(
  sample.type = c(1, 2),
  blank.mean = c(0.3333, 0.3333)
)

for(i in 1:2) {
  lambda_true <-
    exp(transfer.mod.run1$BUGSoutput$sims.list$alpha_gut[, transfersim$sample.type[i]])
  lambda_blanks <- transfersim$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  y <- as.numeric(rpois(lambda_y, lambda_y))
  transfersim$median[i] <- median(lambda_true)
  transfersim$upper25[i] <- quantile(lambda_true, 0.625)
  transfersim$lower25[i] <- quantile(lambda_true, 0.375)
  transfersim$upper50[i] <- quantile(lambda_true, 0.750)
  transfersim$lower50[i] <- quantile(lambda_true, 0.250)
  transfersim$upper75[i] <- quantile(lambda_true, 0.875)
  transfersim$lower75[i] <- quantile(lambda_true, 0.125)
  transfersim$upper95[i] <- quantile(lambda_true, 0.975)
  transfersim$lower95[i] <- quantile(lambda_true, 0.025)
  transfersim$yupper95[i] <- quantile(y, 0.975)
  transfersim$ylower95[i] <- quantile(y, 0.025)
}


transfersim$sample.type <- as.factor(transfersim$sample.type)

transfersim$sample.type <- mapvalues(
  transfersim$sample.type,
  from = levels(transfersim$sample.type),
  to = c("Gut",
         "Gut Animals")
)

#### Plot predictions ####

tiff('Trophic Transfer Bayesian Plot.tiff',
     res = 500,
     width = 9,
     height = 8,
     units = 'cm',
     pointsize = 12)

ggplot() +
  geom_linerange(
    data = transfersim,
    aes(x = sample.type,
        ymax = upper95,
        ymin = lower95),
    alpha = 0.05,
    size = 0.5,
    colour = pal[5]
  ) +
  geom_linerange(
    data = transfersim,
    aes(x = sample.type,
        ymax = upper75,
        ymin = lower75),
    alpha = 0.25,
    size = 0.5,
    colour = pal[5]
  ) +
  geom_linerange(
    data = transfersim,
    aes(x = sample.type,
        ymax = upper50,
        ymin = lower50),
    alpha = 0.5,
    size = 0.5,
    colour = pal[5]
  ) +
  geom_linerange(
    data = transfersim,
    aes(x = sample.type,
        ymax = upper25,
        ymin = lower25),
    alpha = 0.75,
    size = 0.5,
    colour = pal[5]
  ) +
  geom_point(
    data = transfersim,
    aes(x = sample.type,
        y = median),
    size = 1.5,
    colour = pal[1],
    fill = pal[2],
    shape = 21
  ) +
  geom_jitter(
    data = transferdata,
    aes(
      x = sample.type,
      y = count
    ),
    width = 0.25,
    height = 0,
    colour = pal[3],
    size = 1,
    shape = 1,
    alpha = 0.5
  ) +
  labs(x = "",
       y = "Number of Particles") +
  scale_y_continuous(
    expand = c(0, 0.1)
  ) +
  theme1 +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = pal[2]),
        legend.background = element_rect(fill = pal[2]))

dev.off()



