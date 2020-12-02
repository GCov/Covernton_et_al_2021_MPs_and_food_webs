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

pal <- c("#0b3954","#bfd7ea","#ff6663","#e0ff4f","#fefffe")

#### Plankton tow model ####

PTdata_AP <- subset(PT_data2, particle.type != "Natural")
PTdata_AP <- 
  PTdata_AP %>% 
  group_by(ID, site, sample.type, blank.match, sample.volume) %>% 
  summarize(count = sum(count), blank.mean = sum(blank.mean))
PTdata_AP$site <- as.character(PTdata_AP$site)
PTdata_AP$site <- as.factor(PTdata_AP$site) 

PTAPmod <- function() {
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
    alpha_site[j] ~ dnorm(0, 1)
  }
}

## Generate initial values for MCMC

PTAPmodinit <- function()
{
  list(
    "alpha_site" = rnorm(3)
  )
}

## Keep track of parameters

PTAPmodparam <- c("alpha_site")

## Specify data

PTAPmoddata <-
  list(
    y = PTdata_AP$count,
    N = nrow(PTdata_AP),
    lambda_blanks = PTdata_AP$blank.mean,
    volume = PTdata_AP$sample.volume,
    site = as.integer(PTdata_AP$site),
    nsite = length(unique(PTdata_AP$site))
  )

## Run the model
PTAPmodrun1 <- jags.parallel(
  data = PTAPmoddata,
  inits = PTAPmodinit,
  parameters.to.save = PTAPmodparam,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 6193,
  model = PTAPmod
)

PTAPmodrun1
PTAPmodrun1mcmc <- as.mcmc(PTAPmodrun1)
xyplot(PTAPmodrun1mcmc, layout = c(6, ceiling(nvar(PTAPmodrun1mcmc)/6)))

#### Diagnostics ####
PTAPmodparam2 <- c("fitted", "true", "lambda_y", "lambda_true")

PTAPmodrun2 <- jags.parallel(
  data = PTAPmoddata,
  inits = PTAPmodinit,
  parameters.to.save = PTAPmodparam2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 6193,
  model = PTAPmod
)

PTAPmod.response <- t(PTAPmodrun2$BUGSoutput$sims.list$fitted)
PTAPmod.observed <- PTdata_AP$count
PTAPmod.fitted <- apply(t(PTAPmodrun2$BUGSoutput$sims.list$lambda_y),
                      1,
                      median)

check.PTAPmod <- createDHARMa(
  simulatedResponse = PTAPmod.response,
  observedResponse = PTAPmod.observed,
  fittedPredictedResponse = PTAPmod.fitted,
  integerResponse = T
)

plot(check.PTAPmod)
testDispersion(check.PTAPmod)
testZeroInflation(check.PTAPmod)

#### Inference ####

PTAPmodrun1long <- extract.post(PTAPmodrun1)

PTAPmodrun1long$variable <- mapvalues(
  PTAPmodrun1long$variable,
  from = levels(PTAPmodrun1long$variable),
  to = c("Coles Bay",
         "Elliott Bay",
         "Victoria Harbour")
)

PTAPmodrun1long$order <- c(nrow(PTAPmodrun1long):1)

png(
  'PTAP Model Posteriors.png',
  width = 9,
  height = 7,
  units = 'cm',
  res = 500
)

ggplot(PTAPmodrun1long) +
  geom_density_ridges(
    aes(x = exp(value),
        y = reorder(variable, order, mean)),
    fill = pal[3],
    colour = pal[1],
    alpha = 0.5, 
    size = 0.25
  ) +
  coord_cartesian(xlim = c(0, 0.2)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()

#### Predictions ####

## Extract 'true' estimate

PTdata_AP$true.est <-
  apply(PTAPmodrun2$BUGSoutput$sims.list$true, 2, mean)
PTdata_AP$true.est.upper95 <-
  apply(PTAPmodrun2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.975)
PTdata_AP$true.est.lower95 <-
  apply(PTAPmodrun2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.025)

set.seed(5126)

PTAP_sim <- data.frame(
  site = c(1:3)
)

for(i in 1:3){
  y_true <-
    exp(
      PTAPmodrun1$BUGSoutput$sims.list$alpha_site[, PTAP_sim$site[i]]
    )
  PTAP_sim$mean[i] <- mean(y_true)
  PTAP_sim$upper25[i] <- quantile(y_true, 0.625)
  PTAP_sim$lower25[i] <- quantile(y_true, 0.375)
  PTAP_sim$upper50[i] <- quantile(y_true, 0.75)
  PTAP_sim$lower50[i] <- quantile(y_true, 0.25)
  PTAP_sim$upper75[i] <- quantile(y_true, 0.875)
  PTAP_sim$lower75[i] <- quantile(y_true, 0.125)
  PTAP_sim$upper95[i] <- quantile(y_true, 0.975)
  PTAP_sim$lower95[i] <- quantile(y_true, 0.025)
}

PTAP_sim$site <- as.factor(PTAP_sim$site)

PTAP_sim$site <- mapvalues(PTAP_sim$site,
                         from = levels(PTAP_sim$site),
                         to = c("Coles Bay",
                                "Elliot Bay",
                                "Victoria Harbour"))

tiff('Plankton Tows AP Plot.tiff',
     res = 500,
     width = 9,
     height = 8,
     units = 'cm',
     pointsize = 12)

set.seed(123)

ggplot() +
  geom_linerange(data = PTAP_sim,
                 aes(x = site,
                     ymax = upper95,
                     ymin = lower95),
                 alpha = 0.05,
                 size = 0.5,
                 colour = pal[1]) +
  geom_linerange(data = PTAP_sim,
                 aes(x = site,
                     ymax = upper75,
                     ymin = lower75),
                 alpha = 0.25,
                 size = 0.5,
                 colour = pal[1]) +
  geom_linerange(data = PTAP_sim,
                 aes(x = site,
                     ymax = upper50,
                     ymin = lower50),
                 alpha = 0.5,
                 size = 0.5,
                 colour = pal[1]) +
  geom_linerange(data = PTAP_sim,
                 aes(x = site,
                     ymax = upper25,
                     ymin = lower25),
                 alpha = 0.75,
                 size = 0.5,
                 colour = pal[1]) +
  geom_point(data = PTAP_sim,
             aes(x = site,
                 y = mean),
             size = 2,
             colour = pal[1],) +
  geom_jitter(data = PTdata_AP,
              aes(x = site,
                  y = count/sample.volume),
              size = 0.75, shape = 1, colour = pal[3],
              height = 0,
              width = 0.1) +
  labs(x = 'Site',
       y = expression(paste('Particles '*L^-1))) +
  scale_y_continuous(limits = c(0, 0.2),
                     expand = c(0, 0.005),
                     breaks = seq(from = 0,
                                  to = 0.2,
                                  by = 0.05)) +
  theme1

dev.off()


#### Plankton jar model ####

PJdata_AP <- subset(PJ_data2, particle.type != "Natural")
PJdata_AP <- 
  PJdata_AP %>% 
  group_by(ID, site, sample.type, blank.match) %>% 
  summarize(count = sum(count), blank.mean = sum(blank.mean))
PJdata_AP$site <- as.character(PJdata_AP$site)
PJdata_AP$site <- as.factor(PJdata_AP$site) 

PJAPmod <- function() {
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

PJAPmodinit <- function()
{
  list(
    "alpha_site" = rnorm(3)
  )
}

## Keep track of parameters

PJAPmodparam <- c("alpha_site")

## Specify data

PJAPmoddata <-
  list(
    y = PJdata_AP$count,
    N = nrow(PJdata_AP),
    lambda_blanks = PJdata_AP$blank.mean,
    site = as.integer(PJdata_AP$site),
    nsite = length(unique(PJdata_AP$site))
  )

## Run the model
PJAPmodrun1 <- jags.parallel(
  data = PJAPmoddata,
  inits = PJAPmodinit,
  parameters.to.save = PJAPmodparam,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 6193,
  model = PJAPmod
)

PJAPmodrun1
PJAPmodrun1mcmc <- as.mcmc(PJAPmodrun1)
xyplot(PJAPmodrun1mcmc, layout = c(6, ceiling(nvar(PJAPmodrun1mcmc)/6)))

#### Diagnostics ####
PJAPmodparam2 <- c("fitted", "true", "lambda_y")

PJAPmodrun2 <- jags.parallel(
  data = PJAPmoddata,
  inits = PJAPmodinit,
  parameters.to.save = PJAPmodparam2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 6193,
  model = PJAPmod
)

PJAPmod.response <- t(PJAPmodrun2$BUGSoutput$sims.list$fitted)
PJAPmod.observed <- PJdata_AP$count
PJAPmod.fitted <- apply(t(PJAPmodrun2$BUGSoutput$sims.list$lambda_y),
                      1,
                      median)

check.PJAPmod <- createDHARMa(
  simulatedResponse = PJAPmod.response,
  observedResponse = PJAPmod.observed,
  fittedPredictedResponse = PJAPmod.fitted,
  integerResponse = T
)

plot(check.PJAPmod, asFactor = TRUE)
testDispersion(check.PJAPmod)
testZeroInflation(check.PJAPmod)
plotResiduals(check.PJAPmod, PJdata_AP$site)

#### Inference ####

PJAPmodrun1long <- extract.post(PJAPmodrun1)

PJAPmodrun1long$variable <- mapvalues(
  PJAPmodrun1long$variable,
  from = levels(PJAPmodrun1long$variable),
  to = c("Coles Bay",
         "Elliott Bay",
         "Victoria Harbour")
)

PJAPmodrun1long$order <- c(nrow(PJAPmodrun1long):1)

png(
  'MP PJAP Model Posteriors.png',
  width = 9,
  height = 7,
  units = 'cm',
  res = 500
)

ggplot(PJAPmodrun1long) +
  geom_density_ridges(
    aes(x = exp(value),
        y = reorder(variable, order, mean)),
    fill = pal[3],
    colour = pal[1],
    alpha = 0.5, 
    size = 0.25
  ) +
  coord_cartesian(xlim = c(0, 10)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

## Extract 'true' estimate

PJdata_AP$true.est <- apply(PJAPmodrun2$BUGSoutput$sims.list$true, 2, mean)
PJdata_AP$true.est.upper95 <- apply(PJAPmodrun2$BUGSoutput$sims.list$true, 2, quantile, 
                                       probs = 0.975)
PJdata_AP$true.est.lower95 <- apply(PJAPmodrun2$BUGSoutput$sims.list$true, 2, quantile, 
                                       probs = 0.025)

set.seed(5126)

PJAP_sim <- data.frame(
  site = c(1:3),
  blank.mean = sample(PJdata_AP$blank.mean,
                      3,
                      replace = FALSE)
)

for(i in 1:3){
  lambda_true <-
    exp(
      PJAPmodrun1$BUGSoutput$sims.list$alpha_site[, PJAP_sim$site[i]]
    )
  lambda_blanks = PJAP_sim$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  true <- as.numeric(rpois(lambda_true, lambda_true))
  y <- as.numeric(rpois(lambda_y, lambda_y))
  PJAP_sim$mean[i] <- mean(true)
  PJAP_sim$upper25[i] <- quantile(true, 0.625)
  PJAP_sim$lower25[i] <- quantile(true, 0.375)
  PJAP_sim$upper50[i] <- quantile(true, 0.75)
  PJAP_sim$lower50[i] <- quantile(true, 0.25)
  PJAP_sim$upper75[i] <- quantile(true, 0.875)
  PJAP_sim$lower75[i] <- quantile(true, 0.125)
  PJAP_sim$upper95[i] <- quantile(true, 0.975)
  PJAP_sim$lower95[i] <- quantile(true, 0.025)
  PJAP_sim$yupper95[i] <- quantile(y, 0.975)
  PJAP_sim$ylower95[i] <- quantile(y, 0.025)
}

PJAP_sim$site <- as.factor(PJAP_sim$site)

PJAP_sim$site <- mapvalues(PJAP_sim$site,
                         from = levels(PJAP_sim$site),
                         to = c("Coles Bay",
                                "Elliot Bay",
                                "Victoria Harbour"))

tiff('Plankton jars AP Plot.tiff',
     res = 500,
     width = 9,
     height = 8,
     units = 'cm',
     pointsize = 12)

set.seed(123)

ggplot() +
  geom_linerange(data = PJAP_sim,
                 aes(x = site,
                     ymax = upper95,
                     ymin = lower95),
                 alpha = 0.05,
                 size = 0.5,
                 colour = pal[1]) +
  geom_linerange(data = PJAP_sim,
                 aes(x = site,
                     ymax = upper75,
                     ymin = lower75),
                 alpha = 0.25,
                 size = 0.5,
                 colour = pal[1]) +
  geom_linerange(data = PJAP_sim,
                 aes(x = site,
                     ymax = upper50,
                     ymin = lower50),
                 alpha = 0.5,
                 size = 0.5,
                 colour = pal[1]) +
  geom_linerange(data = PJAP_sim,
                 aes(x = site,
                     ymax = upper25,
                     ymin = lower25),
                 alpha = 0.75,
                 size = 0.5,
                 colour = pal[1]) +
  geom_point(data = PJAP_sim,
             aes(x = site,
                 y = mean),
             size = 2,
             colour = pal[1],) +
  geom_jitter(data = PJdata_AP,
              aes(x = site,
                  y = count),
              size = 0.75, shape = 1, colour = pal[3],
              height = 0,
              width = 0.1) +
  labs(x = 'Site',
       y = expression(paste('Particles '*L^-1))) +
  scale_y_continuous(limits = c(0, 30),
                     expand = c(0, 0.5),
                     breaks = seq(from = 0,
                                  to = 30,
                                  by = 5)) +
  theme1

dev.off()


#### MP Model by Individual  ####  

APgutdata <- subset(gutdata, !is.na(trophic.position) & 
                      particle.type != 'Natural')
APgutdata <-
  APgutdata %>% 
  group_by(ID, site, sample.type, blank.match, tissue.wet.weight,
           total.body.wet.weight, species, TL, base_deltaN,
           sd_base_deltaN, deltaN, deltaC, tissue.dry.weight) %>% 
  summarize(count = sum(count), blank.mean = sum(blank.mean))

APgutdata$site <- as.character(APgutdata$site)
APgutdata$site <- as.factor(APgutdata$site)
APgutdata$species <- as.character(APgutdata$species)
APgutdata$species <- as.factor(APgutdata$species)

APmodel1 <- function() {
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

APmodel1init <- function()
{
  list(
    "sigma_species" = rexp(1),
    "beta_TP" = rnorm(3),
    "gamma_site" = rnorm(3),
    "base" = rgamma(3, 1, 1)
  )
}

## Keep track of parameters

APmodel1param <- c("alpha_species", "beta_TP", "gamma_site", "base")

## Specify data

beta_zero <- 5.92
beta_one <- -0.27
nit_lim <- -beta_zero/beta_one
k <- -log((beta_zero - nit_lim)/(-nit_lim))

APmodel1data <-
  list(
    y = APgutdata$count,
    N = nrow(APgutdata),
    lambda_blanks = APgutdata$blank.mean,
    species = as.integer(APgutdata$species),
    nspecies = length(unique(APgutdata$species)),
    site = as.integer(APgutdata$site),
    diff = APgutdata$deltaN - APgutdata$base_deltaN,
    nsite = length(unique(APgutdata$site)),
    deltaN = APgutdata$deltaN,
    mean_base = as.numeric(with(
      APgutdata,
      tapply(base_deltaN,
             as.integer(site),
             mean)
    )),
    sd_base = as.numeric(with(
      APgutdata, tapply(sd_base_deltaN, as.integer(site), mean)
    )),
    nit_lim = nit_lim,
    k = k
  )

## APrun the model
APrun1 <- jags.parallel(
  data = APmodel1data,
  inits = APmodel1init,
  parameters.to.save = APmodel1param,
  n.chains = 3,
  n.cluster = 3,
  n.iter = 7000,
  n.burnin = 500,
  n.thin = 2,
  jags.seed = 3234,
  model = APmodel1
)

APrun1
APrun1mcmc <- as.mcmc(APrun1)
xyplot(APrun1mcmc, layout = c(6, ceiling(nvar(APrun1mcmc)/6)))

#### Diagnostics ####
APmodel1param2 <- c("fitted", "true", "lambda_y", "TP")

APrun2 <- jags.parallel(
  data = APmodel1data,
  inits = APmodel1init,
  parameters.to.save = APmodel1param2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 7000,
  n.burnin = 500,
  n.thin = 2,
  jags.seed = 3234,
  model = APmodel1
)

APmodel1.response <- t(APrun2$BUGSoutput$sims.list$fitted)
APmodel1.observed <- APgutdata$count
APmodel1.fitted <- apply(t(APrun2$BUGSoutput$sims.list$lambda_y),
                       1,
                       median)

check.APmodel1 <- createDHARMa(simulatedResponse = APmodel1.response,
                             observedResponse = APmodel1.observed, 
                             fittedPredictedResponse = APmodel1.fitted,
                             integerResponse = T)

plot(check.APmodel1)

plotResiduals(check.APmodel1, APgutdata$site)
plotResiduals(check.APmodel1, APgutdata$species)
plotResiduals(check.APmodel1, apply(APrun2$BUGSoutput$sims.list$TP, 2, median))
testZeroInflation(check.APmodel1)
testDispersion(check.APmodel1)

plot(APmodel1.observed-APmodel1.fitted ~ log(APmodel1.fitted))

#### Inference ####

APrun1long <- extract.post(APrun1)

APrun1long$variable <- mapvalues(APrun1long$variable,
                               from = levels(APrun1long$variable),
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

APrun1long$order <- c(nrow(APrun1long):1)

png(
  'AP Gut Model Posteriors.png',
  width = 16,
  height = 12,
  units = 'cm',
  res = 500
)

ggplot(APrun1long) +
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
    colour = pal[3]
  ) +
  coord_cartesian(xlim = c(-2, 13)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

## Extract 'true' estimate

APgutdata$true.est <- apply(APrun2$BUGSoutput$sims.list$true, 2, mean)
APgutdata$true.est.upper95 <- apply(APrun2$BUGSoutput$sims.list$true, 2, quantile, 
                                    probs = 0.975)
APgutdata$true.est.lower95 <- apply(APrun2$BUGSoutput$sims.list$true, 2, quantile, 
                                    probs = 0.025)
APgutdata$TP.est <- apply(APrun2$BUGSoutput$sims.list$TP, 2, mean)
APgutdata$TP.est.lower95 <- apply(APrun2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.025)
APgutdata$TP.est.upper95 <- apply(APrun2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.975)
APgutdata$TP.est.lower95 <- apply(APrun2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.025)
APgutdata$TP.est.upper75 <- apply(APrun2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.875)
APgutdata$TP.est.lower75 <- apply(APrun2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.125)
APgutdata$TP.est.upper50 <- apply(APrun2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.75)
APgutdata$TP.est.lower50 <- apply(APrun2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.25)
APgutdata$TP.est.upper25 <- apply(APrun2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.625)
APgutdata$TP.est.lower25 <- apply(APrun2$BUGSoutput$sims.list$TP, 2, quantile,
                                  probs = 0.375)

set.seed(5126)

APgutsim <- data.frame(
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
  blank.mean = sample(APgutdata$blank.mean,
                      2000,
                      replace = TRUE)
)

for(i in 1:2000){
  lambda_true <-
    exp(
      APrun1$BUGSoutput$sims.list$beta_TP[, APgutsim$site[i]]*
        APgutsim$trophic.position[i] +
        APrun1$BUGSoutput$sims.list$gamma_site[, APgutsim$site[i]]
    )
  lambda_blanks = APgutsim$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  true <- as.numeric(rpois(lambda_true, lambda_true))
  y <- as.numeric(rpois(lambda_y, lambda_y))
  APgutsim$mean[i] <- mean(lambda_true)
  APgutsim$upper25[i] <- quantile(lambda_true, 0.625)
  APgutsim$lower25[i] <- quantile(lambda_true, 0.375)
  APgutsim$upper50[i] <- quantile(lambda_true, 0.75)
  APgutsim$lower50[i] <- quantile(lambda_true, 0.25)
  APgutsim$upper75[i] <- quantile(lambda_true, 0.875)
  APgutsim$lower75[i] <- quantile(lambda_true, 0.125)
  APgutsim$upper95[i] <- quantile(lambda_true, 0.975)
  APgutsim$lower95[i] <- quantile(lambda_true, 0.025)
  APgutsim$yupper95[i] <- quantile(y, 0.975)
  APgutsim$ylower95[i] <- quantile(y, 0.025)
}

APgutsim$site <- as.factor(APgutsim$site)

APgutsim$site <- mapvalues(APgutsim$site,
                           from = levels(APgutsim$site),
                           to = c("Coles Bay",
                                  "Elliot Bay",
                                  "Victoria Harbour"))

#### Plot predictions ####

tiff('Trophic Position AP Plot.tiff',
     res = 500,
     width = 16,
     height = 12,
     units = 'cm',
     pointsize = 12)

ggplot() +
  geom_ribbon(data = APgutsim,
              aes(x = trophic.position,
                  ymax = upper95,
                  ymin = lower95),
              alpha = 0.05,
              size = 0.5,
              fill = pal[4]) +
  geom_ribbon(data = APgutsim,
              aes(x = trophic.position,
                  ymax = upper75,
                  ymin = lower75),
              alpha = 0.25,
              size = 0.5,
              fill = pal[4]) +
  geom_ribbon(data = APgutsim,
              aes(x = trophic.position,
                  ymax = upper50,
                  ymin = lower50),
              alpha = 0.5,
              size = 0.5,
              fill = pal[4]) +
  geom_ribbon(data = APgutsim,
              aes(x = trophic.position,
                  ymax = upper25,
                  ymin = lower25),
              alpha = 0.75,
              size = 0.5,
              fill = pal[4]) +
  geom_line(data = APgutsim,
            aes(x = trophic.position,
                y = mean),
            size = 0.5,
            colour = pal[1]) +
  geom_point(data = APgutdata,
             aes(x = TP.est,
                 y = count),
             size = 1, shape = 20, alpha = 0.8, colour = pal[1]) +
  facet_wrap(~ site) +
  labs(x = 'Trophic Position',
       y = expression(paste('Particles '*ind^-1))) +
  coord_cartesian(xlim = c(1, 6)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 12),
                     expand = c(0, 0.1),
                     breaks = c(seq(0, 12, 2))) +
  theme1

dev.off()



## AP concentration by species

set.seed(6614)

APgutsim2 <- expand.grid(
  trophic.position = 2.9,
  site = 1,
  species = c(1:14),
  blank.mean = mean(APgutdata$blank.mean)
)

for(i in 1:14){
  lambda_true <-
    exp(
      APrun1$BUGSoutput$sims.list$beta_TP[, APgutsim2$site[i]]*
        APgutsim2$trophic.position[i] +
        APrun1$BUGSoutput$sims.list$gamma_site[, APgutsim2$site[i]]+
        APrun1$BUGSoutput$sims.list$alpha_species[, APgutsim2$species[i]]
    )
  lambda_blanks = APgutsim2$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  true <- as.numeric(rpois(lambda_true, lambda_true))
  y <- as.numeric(rpois(lambda_y, lambda_y))
  APgutsim2$mean[i] <- mean(lambda_true)
  APgutsim2$upper25[i] <- quantile(lambda_true, 0.625)
  APgutsim2$lower25[i] <- quantile(lambda_true, 0.375)
  APgutsim2$upper50[i] <- quantile(lambda_true, 0.75)
  APgutsim2$lower50[i] <- quantile(lambda_true, 0.25)
  APgutsim2$upper75[i] <- quantile(lambda_true, 0.875)
  APgutsim2$lower75[i] <- quantile(lambda_true, 0.125)
  APgutsim2$upper95[i] <- quantile(lambda_true, 0.975)
  APgutsim2$lower95[i] <- quantile(lambda_true, 0.025)
  APgutsim2$yupper95[i] <- quantile(y, 0.975)
  APgutsim2$ylower95[i] <- quantile(y, 0.025)
}

APgutsim2$species <- as.factor(APgutsim2$species)
APgutsim2$site <- as.factor(APgutsim2$site)

APgutsim2$species <- mapvalues(
  APgutsim2$species,
  from = levels(APgutsim2$species),
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

tiff('AP Species Plot.tiff',
     res = 500,
     width = 14,
     height = 8,
     units = 'cm',
     pointsize = 12)

ggplot() +
  geom_linerange(
    data = APgutsim2,
    aes(x = species,
        ymax = upper95,
        ymin = lower95),
    alpha = 0.05,
    size = 0.5,
    colour = pal[1]
  ) +
  geom_linerange(
    data = APgutsim2,
    aes(x = species,
        ymax = upper75,
        ymin = lower75),
    alpha = 0.25,
    size = 0.5,
    colour = pal[1]
  ) +
  geom_linerange(
    data = APgutsim2,
    aes(x = species,
        ymax = upper50,
        ymin = lower50),
    alpha = 0.5,
    size = 0.5,
    colour = pal[1]
  ) +
  geom_linerange(
    data = APgutsim2,
    aes(x = species,
        ymax = upper25,
        ymin = lower25),
    alpha = 0.75,
    size = 0.5,
    colour = pal[1]
  ) +
  geom_point(
    data = APgutsim2,
    aes(x = species,
        y = mean),
    size = 1.5,
    colour = pal[1],
    fill = pal[5],
    shape = 21
  ) +
  geom_jitter(
    data = APgutdata,
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
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()



#### Body size model with just fish ####

fishgutAPdata <- subset(APgutdata, 
                      sample.type == "Flatfish" |
                        sample.type == "Rockfish" |
                        sample.type == "Surfperch")

fishgutAPdata$species <- as.character(fishgutAPdata$species)
fishgutAPdata$species <- as.factor(fishgutAPdata$species)

names(fishgutAPdata)

nrow(fishgutAPdata)

plot(TP.est ~ TL, data = fishgutAPdata)
plot(TL ~ species, data = fishgutAPdata)
plot(count ~ TL, data = fishgutAPdata)
plot(count ~ log(total.body.wet.weight), data = fishgutAPdata)

ggplot(fishgutAPdata) +
  geom_point(aes(x = log(TL),
                 y = log(total.body.wet.weight),
                 colour = species)) +
  geom_smooth(aes(x = log(TL),
                  y = log(total.body.wet.weight)),
              method = 'lm') +
  facet_grid(. ~ site)


fishAPmodel1 <- function() {
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

fishAPmodel1init <- function()
{
  list(
    "beta_length" = rnorm(3),
    "gamma_site" = rnorm(3),
    "sigma_species" = rexp(1)
  )
}

## Keep track of parameters

fishAPmodel1param <- c("alpha_species", "beta_length", "gamma_site")

## Specify data

fishAPmodel1data <-
  list(
    y = fishgutAPdata$count,
    N = nrow(fishgutAPdata),
    lambda_blanks = fishgutAPdata$blank.mean,
    species = as.integer(fishgutAPdata$species),
    nspecies = length(unique(fishgutAPdata$species)),
    site = as.integer(fishgutAPdata$site),
    nsite = length(unique(fishgutAPdata$site)),
    length = as.numeric(scale(fishgutAPdata$TL), center = TRUE)
  )

## Run the model
fishAPrun1 <- jags.parallel(
  data = fishAPmodel1data,
  inits = fishAPmodel1init,
  parameters.to.save = fishAPmodel1param,
  n.chains = 3,
  n.cluster = 3,
  n.iter = 5000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3234,
  model = fishAPmodel1
)

fishAPrun1
fishAPrun1mcmc <- as.mcmc(fishAPrun1)
xyplot(fishAPrun1mcmc, layout = c(6, ceiling(nvar(fishAPrun1mcmc)/6)))

#### Diagnostics ####
fishAPmodel1param2 <- c("fitted", "true", "lambda_y", "TP")

fishAPrun2 <- jags.parallel(
  data = fishAPmodel1data,
  inits = fishAPmodel1init,
  parameters.to.save = fishAPmodel1param2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 5000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3234,
  model = fishAPmodel1
)

fishAPmodel1.response <- t(fishAPrun2$BUGSoutput$sims.list$fitted)
fishAPmodel1.observed <- fishgutAPdata$count
fishAPmodel1.fitted <- apply(t(fishAPrun2$BUGSoutput$sims.list$lambda_y),
                           1,
                           median)

check.fishAPmodel1 <- createDHARMa(simulatedResponse = fishAPmodel1.response,
                                 observedResponse = fishAPmodel1.observed, 
                                 fittedPredictedResponse = fishAPmodel1.fitted,
                                 integerResponse = T)

plot(check.fishAPmodel1)

plotResiduals(check.fishAPmodel1, fishgutAPdata$site)
plotResiduals(check.fishAPmodel1, fishgutAPdata$species)
plotResiduals(check.fishAPmodel1, fishgutAPdata$TL)
plotResiduals(check.fishAPmodel1, log(fishgutAPdata$total.body.wet.weight))
testZeroInflation(check.fishAPmodel1)
testDispersion(check.fishAPmodel1)

plot(fishAPmodel1.observed-fishAPmodel1.fitted ~ log(fishAPmodel1.fitted))

#### Inference ####

fishAPrun1long <- extract.post(fishAPrun1)

fishAPrun1long$variable <- mapvalues(fishAPrun1long$variable,
                                   from = levels(fishAPrun1long$variable),
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

fishAPrun1long$order <- c(nrow(fishAPrun1long):1)

png(
  'Fish Gut Body Size Model Posteriors.png',
  width = 16,
  height = 12,
  units = 'cm',
  res = 500
)

ggplot(fishAPrun1long) +
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
    colour = pal[3]
  ) +
  coord_cartesian(xlim = c(-3, 2.1)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

## Extract 'true' estimate

fishgutAPdata$true.est <-
  apply(fishAPrun2$BUGSoutput$sims.list$true, 2, mean)
fishgutAPdata$true.est.upper95 <-
  apply(fishAPrun2$BUGSoutput$sims.list$true, 2,
        quantile, probs = 0.975)
fishgutAPdata$true.est.lower95 <-
  apply(fishAPrun2$BUGSoutput$sims.list$true, 2,
        quantile,
        probs = 0.025)

set.seed(5126)

fishgutAPsim <- data.frame(
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
  blank.mean = sample(fishgutAPdata$blank.mean,
                      2000,
                      replace = TRUE)
)

fishgutAPsim$length.stand <-
  (log(fishgutAPsim$length) - log(mean(fishgutAPdata$TL))) /
  sd(log(fishgutAPsim$length) - log(mean(fishgutAPdata$TL)))

for(i in 1:2000){
  lambda_true <-
    exp(
      fishAPrun1$BUGSoutput$sims.list$beta_length[, fishgutAPsim$site[i]]*
        fishgutAPsim$length.stand[i] +
        fishAPrun1$BUGSoutput$sims.list$gamma_site[, fishgutAPsim$site[i]]
    )
  lambda_blanks = fishgutAPsim$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  true <- as.numeric(rpois(lambda_true, lambda_true))
  y <- as.numeric(rpois(lambda_y, lambda_y))
  fishgutAPsim$mean[i] <- mean(lambda_true)
  fishgutAPsim$upper25[i] <- quantile(lambda_true, 0.625)
  fishgutAPsim$lower25[i] <- quantile(lambda_true, 0.375)
  fishgutAPsim$upper50[i] <- quantile(lambda_true, 0.75)
  fishgutAPsim$lower50[i] <- quantile(lambda_true, 0.25)
  fishgutAPsim$upper75[i] <- quantile(lambda_true, 0.875)
  fishgutAPsim$lower75[i] <- quantile(lambda_true, 0.125)
  fishgutAPsim$upper95[i] <- quantile(lambda_true, 0.975)
  fishgutAPsim$lower95[i] <- quantile(lambda_true, 0.025)
  fishgutAPsim$yupper95[i] <- quantile(y, 0.975)
  fishgutAPsim$ylower95[i] <- quantile(y, 0.025)
}

fishgutAPsim$site <- as.factor(fishgutAPsim$site)

fishgutAPsim$site <- mapvalues(fishgutAPsim$site,
                             from = levels(fishgutAPsim$site),
                             to = c("Coles Bay",
                                    "Elliot Bay",
                                    "Victoria Harbour"))

#### Plot predictions ####

tiff('Body Size Fish AP Plot.tiff',
     res = 500,
     width = 16,
     height = 12,
     units = 'cm',
     pointsize = 12)

ggplot() +
  geom_ribbon(data = fishgutAPsim,
              aes(x = length,
                  ymax = upper95,
                  ymin = lower95),
              alpha = 0.05,
              size = 0.5,
              fill = pal[2]) +
  geom_ribbon(data = fishgutAPsim,
              aes(x = length,
                  ymax = upper75,
                  ymin = lower75),
              alpha = 0.25,
              size = 0.5,
              fill = pal[2]) +
  geom_ribbon(data = fishgutAPsim,
              aes(x = length,
                  ymax = upper50,
                  ymin = lower50),
              alpha = 0.5,
              size = 0.5,
              fill = pal[2]) +
  geom_ribbon(data = fishgutAPsim,
              aes(x = length,
                  ymax = upper25,
                  ymin = lower25),
              alpha = 0.75,
              size = 0.5,
              fill = pal[2],
              colour = pal[2]) +
  geom_line(data = fishgutAPsim,
            aes(x = length,
                y = mean),
            size = 0.5,
            colour = pal[1],
            alpha = 0.3) +
  geom_point(data = fishgutAPdata,
             aes(x = TL,
                 y = count),
             size = 0.75, shape = 1, alpha = 0.8, colour = pal[1]) +
  geom_linerange(data = fishgutAPdata,
                 aes(x = TL,
                     ymin = true.est.lower95,
                     ymax = true.est.upper95),
                 size = 0.5, alpha = 0.5, colour = pal[3]) +
  geom_point(data = fishgutAPdata,
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

APliverdata <- subset(liverdata, !is.na(trophic.position) & 
                        particle.type != 'Natural')
APliverdata <- 
  APliverdata %>% 
  group_by(ID, site, sample.type, blank.match, tissue.wet.weight,
           tissue.dry.weight, total.body.wet.weight, animal.type,
           species, TL, base_deltaN, sd_base_deltaN, deltaN, deltaC) %>% 
  summarize(count = sum(count), blank.mean = sum(blank.mean))

APliverdata$site <- as.character(APliverdata$site)
APliverdata$site <- as.factor(APliverdata$site)
APliverdata$species <- as.character(APliverdata$species)
APliverdata$species <- as.factor(APliverdata$species)

ggplot(APliverdata) +
  geom_point(aes(x = total.body.wet.weight,
                 y = tissue.dry.weight,
                 color = species)) +
  facet_grid(. ~ site)

ggplot(APliverdata) +
  geom_point(aes(x = tissue.wet.weight,
                 y = tissue.dry.weight,
                 color = species)) +
  facet_grid(. ~ site)

plot(tissue.dry.weight ~ trophic.position, data = APliverdata)

liverAP.mod <- function() {
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

liverAP.mod.init <- function()
{
  list(
    "sigma_species" = rexp(1),
    "beta_TP" = rnorm(3),
    "gamma_site" = rnorm(3),
    "base" = rgamma(3, 1, 1)
  )
}

## Keep track of parameters

liverAP.mod.params <- c("alpha_species", "beta_TP", "gamma_site", "base")

## Specify data

liverAP.mod.data <-
  list(
    y = APliverdata$count,
    N = nrow(APliverdata),
    lambda_blanks = APliverdata$blank.mean,
    weight = APliverdata$tissue.dry.weight,
    species = as.integer(APliverdata$species),
    nspecies = length(unique(APliverdata$species)),
    site = as.integer(APliverdata$site),
    nsite = length(unique(APliverdata$site)),
    deltaN = APliverdata$deltaN,
    mean_base = as.numeric(with(
      APliverdata,
      tapply(base_deltaN,
             as.integer(site),
             mean)
    )),
    sd_base = as.numeric(with(
      APliverdata, tapply(sd_base_deltaN, as.integer(site), mean)
    )),
    nit_lim = nit_lim,
    k = k
  )

## Run the model
liverAP.mod.run1 <- jags.parallel(
  data = liverAP.mod.data,
  inits = liverAP.mod.init,
  parameters.to.save = liverAP.mod.params,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 5000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3149,
  model = liverAP.mod
)

liverAP.mod.run1
liverAP.mod.run1mcmc <- as.mcmc(liverAP.mod.run1)
xyplot(liverAP.mod.run1mcmc, layout = c(6, ceiling(nvar(liverAP.mod.run1mcmc)/6)))

#### Diagnostics ####
liverAP.mod.params2 <- c("fitted", "true", "lambda_y", "TP")

liverAP.mod.run2 <- jags.parallel(
  data = liverAP.mod.data,
  inits = liverAP.mod.init,
  parameters.to.save = liverAP.mod.params2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 5000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3149,
  model = liverAP.mod
)

liverAP.mod.response <- t(liverAP.mod.run2$BUGSoutput$sims.list$fitted)
liverAP.mod.observed <- APliverdata$count
liverAP.mod.fitted <- apply(t(liverAP.mod.run2$BUGSoutput$sims.list$lambda_y),
                          1,
                          median)

check.liverAP.mod <-
  createDHARMa(
    simulatedResponse = liverAP.mod.response,
    observedResponse = liverAP.mod.observed,
    fittedPredictedResponse = liverAP.mod.fitted,
    integerResponse = T
  )

plot(check.liverAP.mod)

plotResiduals(check.liverAP.mod, 
              apply(liverAP.mod.run2$BUGSoutput$sims.list$TP, 2, median))
plotResiduals(check.liverAP.mod,
              APliverdata$tissue.dry.weight)
plotResiduals(check.liverAP.mod,
              APliverdata$tissue.wet.weight)

testDispersion(check.liverAP.mod)
testZeroInflation(check.liverAP.mod)

#### Inference ####

liverAP.mod.run1long <- extract.post(liverAP.mod.run1)

liverAP.mod.run1long$variable <-
  mapvalues(
    liverAP.mod.run1long$variable,
    from = levels(liverAP.mod.run1long$variable),
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

liverAP.mod.run1long$order <- c(nrow(liverAP.mod.run1long):1)

png(
  'AP Liver Model Posteriors.png',
  width = 9,
  height = 9,
  units = 'cm',
  res = 500
)

ggplot(liverAP.mod.run1long) +
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
    colour = pal[3]
  ) +
  coord_cartesian(xlim = c(-6, 14)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

## Extract 'true' estimate

APliverdata$true.est <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$true, 2, mean)
APliverdata$true.est.upper95 <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.975)
APliverdata$true.est.lower95 <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.025)
APliverdata$TP.est <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$TP, 2, mean)
APliverdata$TP.est.lower95 <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.025)
APliverdata$TP.est.upper95 <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.975)
APliverdata$TP.est.lower95 <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.025)
APliverdata$TP.est.upper75 <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.875)
APliverdata$TP.est.lower75 <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.125)
APliverdata$TP.est.upper50 <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.75)
APliverdata$TP.est.lower50 <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.25)
APliverdata$TP.est.upper25 <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.625)
APliverdata$TP.est.lower25 <-
  apply(liverAP.mod.run2$BUGSoutput$sims.list$TP, 2, quantile,
        probs = 0.375)

TP.mod.liver <- 
  lm(log(tissue.dry.weight) ~ trophic.position + species, 
     data = APliverdata)

plot(resid(TP.mod.liver, type = "pearson") ~ fitted(TP.mod.liver))
summary(TP.mod.liver)

set.seed(5126)

APliversim <- data.frame(
  trophic.position = seq(
    from = 2,
    to = 4.5,
    length.out = 2000
  ),
  site = sample(c(1:3),
                2000,
                replace = TRUE),
  species = sample(APliverdata$species,
                   2000,
                   replace = TRUE),
  blank.mean = sample(APliverdata$blank.mean,
                      2000,
                      replace = TRUE)
)

APliversim$weight <-
  exp(rnorm(predict(TP.mod.liver, newdata = APliversim), 0.6266))

APliversim$species <- as.integer(APliversim$species)

for(i in 1:2000){
  true.mean <-
    exp(
      liverAP.mod.run1$BUGSoutput$sims.list$beta_TP[, APliversim$site[i]]*
        APliversim$trophic.position[i] +
        liverAP.mod.run1$BUGSoutput$sims.list$gamma_site[, APliversim$site[i]]
    )
  lambda_true <- true*APliversim$weight[i]
  lambda_blanks <- APliversim$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  true.conc <- as.numeric(rpois(lambda_true, lambda_true))/APliversim$weight[i]
  y <- as.numeric(rpois(lambda_y, lambda_y))
  APliversim$median[i] <- median(true.mean)
  APliversim$upper25[i] <- quantile(true.mean, 0.625)
  APliversim$lower25[i] <- quantile(true.mean, 0.375)
  APliversim$upper50[i] <- quantile(true.mean, 0.750)
  APliversim$lower50[i] <- quantile(true.mean, 0.250)
  APliversim$upper75[i] <- quantile(true.mean, 0.875)
  APliversim$lower75[i] <- quantile(true.mean, 0.125)
  APliversim$upper95[i] <- quantile(true.mean, 0.975)
  APliversim$lower95[i] <- quantile(true.mean, 0.025)
  APliversim$yupper95[i] <- quantile(y, 0.975)
  APliversim$ylower95[i] <- quantile(y, 0.025)
}


APliversim$site <- as.factor(APliversim$site)

APliversim$site <- mapvalues(APliversim$site,
                             from = levels(APliversim$site),
                             to = c("Coles Bay",
                                    "Elliot Bay",
                                    "Victoria Harbour"))

#### Plot predictions ####

tiff('Trophic Position AP Liver Bayesian Plot.tiff',
     res = 500,
     width = 16,
     height = 12,
     units = 'cm',
     pointsize = 12)

ggplot() +
  geom_ribbon(
    data = APliversim,
    aes(x = trophic.position,
        ymax = upper95,
        ymin = lower95),
    alpha = 0.05,
    size = 0.5,
    fill = pal[4]
  ) +
  geom_ribbon(
    data = APliversim,
    aes(x = trophic.position,
        ymax = upper75,
        ymin = lower75),
    alpha = 0.25,
    size = 0.5,
    fill = pal[4]
  ) +
  geom_ribbon(
    data = APliversim,
    aes(x = trophic.position,
        ymax = upper50,
        ymin = lower50),
    alpha = 0.5,
    size = 0.5,
    fill = pal[4]
  ) +
  geom_ribbon(
    data = APliversim,
    aes(x = trophic.position,
        ymax = upper25,
        ymin = lower25),
    alpha = 0.75,
    size = 0.5,
    fill = pal[4]
  ) +
  geom_line(
    data = APliversim,
    aes(x = trophic.position,
        y = median),
    size = 0.5,
    colour = pal[1],
    alpha = 0.3
  ) +
  geom_point(
    data = APliverdata,
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
    breaks = c(0, 1, 10, 20)
  ) +
  theme1

dev.off()


#### Rockfish ingested animals vs. gut ####

APtransferdata <- subset(foodweb2, 
                       sample.type == "Rockfish: Ingested Animals" &
                         particle.type != "Natural" |
                         sample.type == "Rockfish" &
                         particle.type != "Natural")
APtransferdata <-
  APtransferdata %>% 
  group_by(ID, site, sample.type, blank.match, tissue.wet.weight, 
           tissue.dry.weight, total.body.wet.weight, species, TL,
           base_deltaN, sd_base_deltaN, deltaN, deltaC) %>% 
  summarize(count = sum(count), blank.mean = sum(blank.mean))

APtransferdata$site <- as.character(APtransferdata$site)
APtransferdata$site <- as.factor(APtransferdata$site)
APtransferdata$species <- as.character(APtransferdata$species)
APtransferdata$species <- as.factor(APtransferdata$species)
APtransferdata$sample.type <- as.character(APtransferdata$sample.type)
APtransferdata$sample.type <- as.factor(APtransferdata$sample.type)
APtransferdata$sample.type <- mapvalues(APtransferdata$sample.type,
                                      from = levels(APtransferdata$sample.type),
                                      to = c("Gut",
                                             "Gut Animals"))

transferAP.mod <- function() {
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

transferAP.mod.init <- function()
{
  list(
    "alpha_gut" = rnorm(2)
  )
}

## Keep track of parameters

transferAP.mod.params <- c("alpha_gut")

## Specify data

transferAP.mod.data <-
  list(
    y = APtransferdata$count,
    N = nrow(APtransferdata),
    lambda_blanks = APtransferdata$blank.mean,
    sample.type = as.integer(APtransferdata$sample.type)
  )

## Run the model
transferAP.mod.run1 <- jags.parallel(
  data = transferAP.mod.data,
  inits = transferAP.mod.init,
  parameters.to.save = transferAP.mod.params,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3242,
  model = transferAP.mod
)

transferAP.mod.run1
transferAP.mod.run1mcmc <- as.mcmc(transferAP.mod.run1)
xyplot(transferAP.mod.run1mcmc, 
       layout = c(6, ceiling(nvar(transferAP.mod.run1mcmc)/6)))

#### Diagnostics ####
transferAP.mod.params2 <- c("fitted", "true", "lambda_y")

transferAP.mod.run2 <- jags.parallel(
  data = transferAP.mod.data,
  inits = transferAP.mod.init,
  parameters.to.save = transferAP.mod.params2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3242,
  model = transferAP.mod
)

transferAP.mod.response <- t(transferAP.mod.run2$BUGSoutput$sims.list$fitted)
transferAP.mod.observed <- APtransferdata$count
transferAP.mod.fitted <- apply(t(transferAP.mod.run2$BUGSoutput$sims.list$lambda_y),
                             1,
                             median)

check.transferAP.mod <-
  createDHARMa(
    simulatedResponse = transferAP.mod.response,
    observedResponse = transferAP.mod.observed,
    fittedPredictedResponse = transferAP.mod.fitted,
    integerResponse = T
  )

plot(check.transferAP.mod)

plotResiduals(check.transferAP.mod, APtransferdata$species)
plotResiduals(check.transferAP.mod, APtransferdata$sample.type)

#### Inference ####

transferAP.mod.run1long <- extract.post(transferAP.mod.run1)

transferAP.mod.run1long$variable <-
  mapvalues(
    transferAP.mod.run1long$variable,
    from = levels(transferAP.mod.run1long$variable),
    to = c(
      "Gut",
      "Gut Animals"
    )
  )

transferAP.mod.run1long$order <- c(nrow(transferAP.mod.run1long):1)

png(
  'AP Transfer Model Posteriors.png',
  width = 9,
  height = 9,
  units = 'cm',
  res = 500
)

ggplot(transferAP.mod.run1long) +
  geom_density_ridges(
    aes(x = exp(value),
        y = reorder(variable, order, mean)),
    fill = pal[3],
    colour = pal[1],
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

APtransferdata$true.est <-
  apply(transferAP.mod.run2$BUGSoutput$sims.list$true, 2, mean)
APtransferdata$true.est.upper95 <-
  apply(transferAP.mod.run2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.975)
APtransferdata$true.est.lower95 <-
  apply(transferAP.mod.run2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.025)

set.seed(3256)

transfersim <- data.frame(
  sample.type = c(1, 2),
  blank.mean = c(0.3333, 0.3333)
)

for(i in 1:2) {
  lambda_true <-
    exp(transferAP.mod.run1$BUGSoutput$sims.list$alpha_gut[, transfersim$sample.type[i]])
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

tiff('AP Transfer Plot.tiff',
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
    colour = pal[3]
  ) +
  geom_linerange(
    data = transfersim,
    aes(x = sample.type,
        ymax = upper75,
        ymin = lower75),
    alpha = 0.25,
    size = 0.5,
    colour = pal[3]
  ) +
  geom_linerange(
    data = transfersim,
    aes(x = sample.type,
        ymax = upper50,
        ymin = lower50),
    alpha = 0.5,
    size = 0.5,
    colour = pal[3]
  ) +
  geom_linerange(
    data = transfersim,
    aes(x = sample.type,
        ymax = upper25,
        ymin = lower25),
    alpha = 0.75,
    size = 0.5,
    colour = pal[3]
  ) +
  geom_point(
    data = transfersim,
    aes(x = sample.type,
        y = median),
    size = 1.5,
    colour = pal[1],
    fill = pal[5],
    shape = 21
  ) +
  geom_jitter(
    data = APtransferdata,
    aes(
      x = sample.type,
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
    expand = c(0, 0.1)
  ) +
  theme1

dev.off()


#### Comparision of empty vs. full guts ####

APtransferdata2 <-
  APtransferdata %>% 
  group_by(ID) %>% 
  summarize(num = length(count)) %>% 
  filter(num > 1)

rfishAPcompare <- subset(APgutdata, sample.type == "Rockfish")

rfishAPcompare$full.stomach <- rfishAPcompare$ID %in% APtransferdata2$ID

rfishAPcompare$ID <- as.character(rfishAPcompare$ID)
rfishAPcompare$ID <- as.factor(rfishAPcompare$ID)

rfishAPcompare$species <- as.character(rfishAPcompare$species)
rfishAPcompare$species <- as.factor(rfishAPcompare$species)

rfishAPcompare$full.stomach <- as.factor(rfishAPcompare$full.stomach)


rfishAP.mod <- function() {
  for (i in 1:N) {
    y[i] ~ dpois(lambda_y[i])
    
    lambda_y[i] <- lambda_true[i] + lambda_blanks[i]
    
    true[i] ~ dpois(lambda_true[i])
    
    log(lambda_true[i]) <- alpha_gut[full.stomach[i]]
    
    ## Fitted values
    fitted[i] ~ dpois(lambda_y[i])
  }
  
  ## Priors
  for (j in 1:2) {
    alpha_gut[j] ~ dnorm(0, 1)
  }
}

## Generate initial values for MCMC

rfishAP.mod.init <- function()
{
  list(
    "alpha_gut" = rnorm(2)
  )
}

## Keep track of parameters

rfishAP.mod.params <- c("alpha_gut")

## Specify data

rfishAP.mod.data <-
  list(
    y = rfishAPcompare$count,
    N = nrow(rfishAPcompare),
    lambda_blanks = rfishAPcompare$blank.mean,
    full.stomach = as.integer(rfishAPcompare$full.stomach)
  )

## Run the model
rfishAP.mod.run1 <- jags.parallel(
  data = rfishAP.mod.data,
  inits = rfishAP.mod.init,
  parameters.to.save = rfishAP.mod.params,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3242,
  model = rfishAP.mod
)

rfishAP.mod.run1
rfishAP.mod.run1mcmc <- as.mcmc(rfishAP.mod.run1)
xyplot(rfishAP.mod.run1mcmc, layout = c(6, ceiling(nvar(rfishAP.mod.run1mcmc)/6)))

#### Diagnostics ####
rfishAP.mod.params2 <- c("fitted", "true", "lambda_y")

rfishAP.mod.run2 <- jags.parallel(
  data = rfishAP.mod.data,
  inits = rfishAP.mod.init,
  parameters.to.save = rfishAP.mod.params2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3242,
  model = rfishAP.mod
)

rfishAP.mod.response <- t(rfishAP.mod.run2$BUGSoutput$sims.list$fitted)
rfishAP.mod.observed <- rfishAPcompare$count
rfishAP.mod.fitted <- apply(t(rfishAP.mod.run2$BUGSoutput$sims.list$lambda_y),
                          1,
                          median)

check.rfishAP.mod <-
  createDHARMa(
    simulatedResponse = rfishAP.mod.response,
    observedResponse = rfishAP.mod.observed,
    fittedPredictedResponse = rfishAP.mod.fitted,
    integerResponse = T
  )

plot(check.rfishAP.mod)
plotResiduals(check.rfishAP.mod, rfishAPcompare$full.stomach)

#### Inference ####

rfishAP.mod.run1long <- extract.post(rfishAP.mod.run1)

rfishAP.mod.run1long$variable <-
  mapvalues(
    rfishAP.mod.run1long$variable,
    from = levels(rfishAP.mod.run1long$variable),
    to = c(
      "Empty Stomach",
      "Full Stomach"
    )
  )

rfishAP.mod.run1long$order <- c(nrow(rfishAP.mod.run1long):1)

png(
  'Rockfish Gut AP Comparison Model Posteriors.png',
  width = 9,
  height = 9,
  units = 'cm',
  res = 500
)

ggplot(rfishAP.mod.run1long) +
  geom_density_ridges(
    aes(x = exp(value),
        y = reorder(variable, order, mean)),
    fill = pal[3],
    colour = pal[1],
    alpha = 0.5, 
    size = 0.25
  ) +
  coord_cartesian(xlim = c(0, 3)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

rfishAPcompare$true.est <-
  apply(rfishAP.mod.run2$BUGSoutput$sims.list$true, 2, mean)
rfishAPcompare$true.est.upper95 <-
  apply(rfishAP.mod.run2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.975)
rfishAPcompare$true.est.lower95 <-
  apply(rfishAP.mod.run2$BUGSoutput$sims.list$true, 2, quantile,
        probs = 0.025)

set.seed(3256)

rfishAPsim <- data.frame(
  full.stomach = c(1, 2),
  blank.mean = c(1, 1)
)

for(i in 1:2) {
  lambda_true <-
    exp(rfishAP.mod.run1$BUGSoutput$sims.list$alpha_gut[, rfishAPsim$full.stomach[i]])
  lambda_blanks <- rfishAPsim$blank.mean[i]
  lambda_y <- lambda_true + lambda_blanks
  y <- as.numeric(rpois(lambda_y, lambda_y))
  rfishAPsim$median[i] <- median(lambda_true)
  rfishAPsim$upper25[i] <- quantile(lambda_true, 0.625)
  rfishAPsim$lower25[i] <- quantile(lambda_true, 0.375)
  rfishAPsim$upper50[i] <- quantile(lambda_true, 0.750)
  rfishAPsim$lower50[i] <- quantile(lambda_true, 0.250)
  rfishAPsim$upper75[i] <- quantile(lambda_true, 0.875)
  rfishAPsim$lower75[i] <- quantile(lambda_true, 0.125)
  rfishAPsim$upper95[i] <- quantile(lambda_true, 0.975)
  rfishAPsim$lower95[i] <- quantile(lambda_true, 0.025)
  rfishAPsim$yupper95[i] <- quantile(y, 0.975)
  rfishAPsim$ylower95[i] <- quantile(y, 0.025)
}


rfishAPsim$full.stomach <- as.factor(rfishAPsim$full.stomach)

rfishAPsim$full.stomach <- mapvalues(
  rfishAPsim$full.stomach,
  from = levels(rfishAPsim$full.stomach),
  to = c("Empty Stomach",
         "Full Stomach")
)

rfishAPcompare$full.stomach <- mapvalues(
  rfishAPcompare$full.stomach,
  from = levels(rfishAPcompare$full.stomach),
  to = c("Empty Stomach",
         "Full Stomach")
)

#### Plot predictions ####

tiff('Rockfish Gut AP Comparison Plot.tiff',
     res = 500,
     width = 9,
     height = 8,
     units = 'cm',
     pointsize = 12)

ggplot() +
  geom_linerange(
    data = rfishAPsim,
    aes(x = full.stomach,
        ymax = upper95,
        ymin = lower95),
    alpha = 0.05,
    size = 0.5,
    colour = pal[3]
  ) +
  geom_linerange(
    data = rfishAPsim,
    aes(x = full.stomach,
        ymax = upper75,
        ymin = lower75),
    alpha = 0.25,
    size = 0.5,
    colour = pal[3]
  ) +
  geom_linerange(
    data = rfishAPsim,
    aes(x = full.stomach,
        ymax = upper50,
        ymin = lower50),
    alpha = 0.5,
    size = 0.5,
    colour = pal[3]
  ) +
  geom_linerange(
    data = rfishAPsim,
    aes(x = full.stomach,
        ymax = upper25,
        ymin = lower25),
    alpha = 0.75,
    size = 0.5,
    colour = pal[3]
  ) +
  geom_point(
    data = rfishAPsim,
    aes(x = full.stomach,
        y = median),
    size = 1.5,
    colour = pal[1],
    fill = pal[3],
    shape = 21
  ) +
  geom_jitter(
    data = rfishAPcompare,
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
  labs(x = "",
       y = "Number of Particles") +
  scale_y_continuous(
    expand = c(0, 0.1)
  ) +
  theme1

dev.off()