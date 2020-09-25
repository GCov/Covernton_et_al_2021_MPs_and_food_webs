library(ggplot2)
library(R2jags)
library(coda)
library(lattice)
library(ggridges)
library(DHARMa)
library(reshape2)

pal <- c("#c52f01",  #rust
         "#b88100",  #dark goldenrod 
         "#0076a9",  #celadon blue
         "#5a8400",  #avocado
         "#011e44")  #oxford blue

MPgutdata <- subset(gutdata, !is.na(trophic.position) & 
                     particle.type == 'Synthetic Polymer')
MPgutdata$particle.type <- as.character(MPgutdata$particle.type)
MPgutdata$particle.type <- as.factor(MPgutdata$particle.type)
MPgutdata$site <- as.character(MPgutdata$site)
MPgutdata$site <- as.factor(MPgutdata$site)
MPgutdata$species <- as.character(MPgutdata$species)
MPgutdata$species <- as.factor(MPgutdata$species)

#### Fit model  ####

model1 <- function() {
  # Likelihood
  for(i in 1:N) {
    y[i] ~ dnegbin(p[i], r)
    p[i] <- r/(r + mu[i])
    mu[i] <- true[i] + contam[i]
    contam[i] ~ dpois(blank.mean[i])
    log(true[i]) <- 
     alpha_species[species[i]] + 
      beta_TP[site[i]]*TP[i] + 
      gamma_site[site[i]]
    
    ## Fitted values
    fitted[i] ~ dnegbin(p[i], r)
  }
  
  ## Priors
  r ~ dexp(1)
  
  for(j in 1:nspecies) {
    alpha_species[j] ~ dnorm(0, tau_species)
  }
  tau_species <- inverse(pow(sigma_species, 2))
  sigma_species ~ dexp(1)
  
  for(k in 1:nsite) {
    beta_TP[k] ~ dnorm(0, 1)
    gamma_site[k] ~ dnorm(0, 1)
  }
}

## Generate initial values for MCMC

model1init <- function()
{
  list(
    "r" = rexp(0.1),
    "sigma_species" = rexp(1),
    "beta_TP" = rnorm(3),
    "gamma_site" = rnorm(3)
  )
}

## Keep track of parameters

model1param <- c("r", "alpha_species", "beta_TP", "gamma_site")

## Specify data

model1data <-
  list(
    y = MPgutdata$orig.count,
    N = nrow(MPgutdata),
    blank.mean = MPgutdata$blank.mean,
    species = as.integer(MPgutdata$species),
    nspecies = length(unique(MPgutdata$species)),
    site = as.integer(MPgutdata$site),
    TP = as.numeric(scale(MPgutdata$trophic.position, center = TRUE)),
    nsite = length(unique(MPgutdata$site))
  )

## Run the model
run1 <- jags.parallel(
  data = model1data,
  inits = model1init,
  parameters.to.save = model1param,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 5000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3234,
  model = model1
)

run1
run1mcmc <- as.mcmc(run1)
xyplot(run1mcmc, layout = c(6, ceiling(nvar(run1mcmc)/6)))

#### Diagnostics ####
model1param2 <- c("fitted", "contam")

run2 <- jags.parallel(
  data = model1data,
  inits = model1init,
  parameters.to.save = model1param2,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 5000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 3234,
  model = model1
)

model1.response <- t(run2$BUGSoutput$sims.list$fitted)
model1.observed <- MPgutdata$orig.count
model1.fitted <- apply(t(run2$BUGSoutput$sims.list$fitted),
                         1,
                         median)

check.model1 <- createDHARMa(simulatedResponse = model1.response,
                               observedResponse = model1.observed, 
                               fittedPredictedResponse = model1.fitted,
                               integerResponse = T)

plot(check.model1)

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
  width = 9,
  height = 9,
  units = 'cm',
  res = 500
)

ggplot(run1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[1],
    colour = pal[5],
    alpha = 0.75
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.5,
    colour = pal[3]
  ) +
  coord_cartesian(xlim = c(-2.5, 3)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()


#### Predictions ####

## Extract 'true' estimate

MPgutdata$true.est <- 
  apply(run2$BUGSoutput$sims.list$fitted - 
          run2$BUGSoutput$sims.list$contam, 2, mean)

set.seed(5126)

MPgutsim <- data.frame(trophic.position = seq(from = 0, 
                                              to = 4, 
                                              length.out = 1000),
                       site = sample(c(1:3), 
                                     1000, 
                                     replace = TRUE),
                       species = sample(c(1:14),
                                        1000,
                                        replace = TRUE))

MPgutsim$stand.trophic.position <-
  (MPgutsim$trophic.position - mean(MPgutdata$trophic.position)) /
  sd(MPgutdata$trophic.position - mean(MPgutdata$trophic.position))

for(i in 1:1000){
  true <-
    exp(
      run1$BUGSoutput$sims.list$beta_TP[, MPgutsim$site[i]] * 
        MPgutsim$stand.trophic.position[i] +
        run1$BUGSoutput$sims.list$gamma_site[, MPgutsim$site[i]]
    )
  r <- run1$BUGSoutput$sims.list$r
  mu <- (p*r)/(1-p)
  contam <- mu - true
  MPgutsim$mean[i] <- mean(true)
  MPgutsim$upper25[i] <- quantile(true, 0.625)
  MPgutsim$lower25[i] <- quantile(true, 0.375)
  MPgutsim$upper50[i] <- quantile(true, 0.75)
  MPgutsim$lower50[i] <- quantile(true, 0.25)
  MPgutsim$upper75[i] <- quantile(true, 0.875)
  MPgutsim$lower75[i] <- quantile(true, 0.125)
  MPgutsim$upper95[i] <- quantile(true, 0.975)
  MPgutsim$lower95[i] <- quantile(true, 0.025)
  MPgutsim$blanksupper95[i] <- quantile(blanks, 0.975)
  MPgutsim$blankslower95[i] <- quantile(blanks, 0.025)
}

MPgutsim$site <- as.factor(MPgutsim$site)

MPgutsim$site <- mapvalues(MPgutsim$site,
                           from = levels(MPgutsim$site),
                           to = c("Coles Bay",
                                  "Elliot Bay",
                                  "Victoria Harbour"))

tiff('Trophic Position MP Bayesian Plot.tiff',
     res = 300,
     width = 16,
     height = 12,
     units = 'cm',
     pointsize = 12)

ggplot() +
  geom_ribbon(data = MPgutsim,
              aes(x = trophic.position,
                  ymax = upper95,
                  ymin = lower95),
              alpha = 0.3,
              size = 0.5,
              fill = pal[3]) +
  geom_line(data = MPgutsim,
            aes(x = trophic.position,
                y = mean),
            linetype = 'dashed', 
            size = 0.5) +
  geom_point(data = MPgutdata,
               aes(x = trophic.position,
                 y = orig.count),
             size = 0.5, shape = 1, alpha = 0.8) +
  geom_point(data = MPgutdata,
             aes(x = trophic.position,
                 y = adj.count),
             size = 0.5, shape = 1, alpha = 0.8, colour = pal[1]) +
  facet_wrap(~ site) +
  labs(x = 'Trophic Position',
       y = expression(paste('Particles '*ind^-1))) +
  theme1

dev.off()
