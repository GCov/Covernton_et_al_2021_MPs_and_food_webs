##### Setup #####

## load packages
library(plyr)
library(ggplot2)
library(dplyr)
library(ggridges)
library(colorspace)
library(glmmTMB)
library(MuMIn)
library(broom.mixed)
library(cowplot)
library(sjPlot)

## Define standard error function

se <- function(x) {
  sd(x)/sqrt(length(x))
}

#### Combine all data ####

## Combine lab and field data with isotopes data

foodweb1 <- 
  left_join(animal_info, 
            isotopes,
            by = c('ID', 'site', 'sample.type'))

foodweb1$sample.type <- as.factor(foodweb1$sample.type)

ingested_animals <- 
  foodweb1 %>% 
  filter(sample.type == 'Rockfish: Ingested Animals') %>% 
  select(1:14)

rockfish_info <- 
  foodweb1 %>% 
  filter(animal.type == 'Rockfish' & sample.type == 'Rockfish') %>% 
  select(1,2,15:36)

ingested_animals2 <-
  left_join(ingested_animals, rockfish_info, by = c('ID', 'site'))

foodweb1 <- rbind(subset(foodweb1, 
                         sample.type != 'Rockfish: Ingested Animals'),
                  ingested_animals2)

## Combine particle counts data with lab, field, and isotopes data

foodweb2 <- left_join(animal_data3,
                      foodweb1,
                      by = c('ID', 'site', 'sample.type'))

foodweb2$ID <- as.factor(foodweb2$ID)
foodweb2$site <- as.factor(foodweb2$site)
foodweb2$sample.type <- as.factor(foodweb2$sample.type)
foodweb2shape <- as.factor(foodweb2$shape)
foodweb2$colour <- as.factor(foodweb2$colour)

#### Summarize ####

## First summarize plankton tow data across size categories, colour and shape

PT_data3 <-
  PT_data2 %>% 
  group_by(ID, site, sample.type, particle.type, sample.volume) %>% 
  summarize(count = sum(adj.count))

## Make sure each sample has all levels of particle type

PT_data3$ID <- as.character(PT_data3$ID)
PT_data3$ID <- as.factor(PT_data3$ID)

all_typesPT <- expand.grid(ID = levels(PT_data3$ID),
                           particle.type = levels(PT_data3$particle.type))

infoPT <- PT_data3[c(1:3, 5)] %>% 
  group_by(ID, site, sample.type) %>% 
  summarize(sample.volume = unique(sample.volume))

all_typesPT2 <- left_join(all_typesPT, infoPT, by = 'ID')

countsPT <- subset(PT_data3[c(1,4,6)], count > 0)

PT_data4 <- left_join(all_typesPT2, countsPT, by = c('ID', 'particle.type'))

PT_data4$count[is.na(PT_data4$count)] <- 0

summary(PT_data4)

## Do the same for PJ

## First summarize plankton jar data across size categories, colour and shape

PJ_data3 <-
  PJ_data2 %>% 
  group_by(ID, site, sample.type, particle.type) %>% 
  summarize(count = sum(adj.count))

## Make sure each sample has all levels of particle type

PJ_data3$ID <- as.character(PJ_data3$ID)
PJ_data3$ID <- as.factor(PJ_data3$ID)

all_typesPJ <- expand.grid(ID = levels(PJ_data3$ID),
                           particle.type = levels(PJ_data3$particle.type))

infoPJ <- PJ_data3[c(1:3)] %>% 
  group_by(ID, site, sample.type) %>% 
  summarize()

all_typesPJ2 <- left_join(all_typesPJ, infoPJ, by = 'ID')

countsPJ <- subset(PJ_data3[c(1,4:5)], count > 0)

PJ_data4 <- left_join(all_typesPJ2, countsPJ, by = c('ID', 'particle.type'))

PJ_data4$count[is.na(PJ_data4$count)] <- 0

summary(PJ_data4)


## Now summarize animal data

## Combine across size categories, colour, and shape

foodweb3 <- 
  foodweb2 %>% 
  group_by(ID, site, sample.type, particle.type, shell.l, shell.w, shell.h,
           arm.length, tissue.wet.weight, tissue.dry.weight, shell.weight,
           total.body.wet.weight, density.sep, species, carapace.length,
           TL, SL, sex, babies, parasites, deltaC, deltaN, trophic.position) %>% 
  summarize(count = sum(adj.count))

## Make sure each sample has all levels of particle type

all_types <- expand.grid(ID = levels(foodweb3$ID), 
                         particle.type = levels(foodweb3$particle.type))

info <- 
  foodweb3[c(1:3, 5:23)] %>% 
  group_by(ID, site, sample.type, shell.l, shell.w, shell.h, arm.length,
           tissue.wet.weight, tissue.dry.weight, shell.weight,
           total.body.wet.weight, density.sep, species, carapace.length,
           TL, SL, sex, babies, parasites, deltaC, deltaN, trophic.position) %>% 
  summarize()

all_types2 <- left_join(all_types, info, by = 'ID')

counts <- subset(foodweb3[c(1,4,24)], count > 0)

foodweb4 <- left_join(all_types2, counts, by = c('ID', 'particle.type'))

foodweb4$count[is.na(foodweb4$count)] <- 0

summary(foodweb4)

foodweb4 <- subset(foodweb4, !is.na(species))

gutdata <- subset(foodweb4,
                  sample.type != 'Surfperch Livers' &
                  sample.type != 'Flatfish Livers' &
                  sample.type != 'Rockfish Livers')

gutdata$sample.type <- 
  mapvalues(gutdata$sample.type,
            from = 'Rockfish: Ingested Animals',
            to = 'Rockfish')

gutdata <-
  gutdata %>% 
  group_by(ID, particle.type, site, sample.type,
           total.body.wet.weight, species, deltaC, deltaN,
           trophic.position) %>% 
  summarize(count = sum(count),
            tissue.weight = sum(tissue.dry.weight))

#### Model ####

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}  # specify overdispersion assessment function

## Model gut data according to indivual

gutdata2 <- subset(gutdata, !is.na(trophic.position) & 
                     particle.type != 'Natural')

Imod1 <- glmmTMB(count ~ trophic.position*site*particle.type + 
                          (1 | sample.type),
                        family = poisson(), data = gutdata2, REML = TRUE)

overdisp_fun(Imod1)  # model is overdispersed

Imod2 <- update(Imod1, ziformula = ~ 1)  # refit with ZI

AICc(Imod1, Imod2)  # better fit with ZI according to AICc

## Diagnostics
plot(resid(Imod2) ~ fitted(Imod2))
plot(resid(Imod2) ~ gutdata2$trophic.position)
plot(resid(Imod2) ~ gutdata2$site)
plot(resid(Imod2) ~ gutdata2$particle.type)
plot(resid(Imod2) ~ gutdata2$sample.type)  # these don't look too bad

Imod3 <- update(Imod2, family = nbinom1(), ziformula = ~ 0)  # try with NB distribution

AICc(Imod2, Imod3)  # better fit as NB model

## Diagnostics
plot(resid(Imod3) ~ fitted(Imod3))
plot(resid(Imod3) ~ gutdata2$trophic.position)
plot(resid(Imod3) ~ gutdata2$site)
plot(resid(Imod3) ~ gutdata2$particle.type)
plot(resid(Imod3) ~ gutdata2$sample.type)  # these don't look too bad

## Try dropping terms

drop1(Imod3)  # can drop 3-way interaction

Imod4 <- update(Imod3, . ~ . -trophic.position:site:particle.type)

AICc(Imod3, Imod4)  # confirms better fit

drop1(Imod4)  # can drop trophic.position:site

Imod5 <- update(Imod4, . ~ . -trophic.position:site)

AICc(Imod4, Imod5)  # Imod5 is better fit

drop1(Imod5)  # can drop trophic.position:particle.type

Imod6 <- update(Imod5, . ~ . -trophic.position:particle.type)

AICc(Imod5, Imod6)  # Imod6 is better fit

drop1(Imod6)  # can drop trophic.position

Imod7 <- update(Imod6, . ~ . -trophic.position)

AICc(Imod6, Imod7)  # Imod7 is better fit

drop1(Imod7) # stop here

Imod7$call  # view final model formula

Imod8 <- update(Imod7, REML = TRUE)  # refit with REML


## Diagnostics
plot(resid(Imod8) ~ fitted(Imod8))
plot(resid(Imod8) ~ gutdata2$trophic.position)
plot(resid(Imod8) ~ gutdata2$site)
plot(resid(Imod8) ~ gutdata2$particle.type)
plot(resid(Imod8) ~ gutdata2$sample.type)  # these don't look too bad (?)


## Inference

tidy(Imod8)


## Predict

Imod.pred <- gutdata2

preds1 <- predict(Imod8, type = 'conditional', re.form = NA, se.fit = TRUE)

Imod.pred$fitted <- preds1$fit

Imod.pred$lower <- with(Imod.pred, fitted - 1.96*preds1$se.fit)
Imod.pred$upper <- with(Imod.pred, fitted + 1.96*preds1$se.fit)  
# 95% CI for preds
# Note that predictions are for averaged REs and ZI



## Now model according to gut dry weight

Cmod1 <- glmmTMB(
  count ~ trophic.position * site * particle.type +
    (1 | sample.type) + offset(log(tissue.weight)),
  family = poisson(link = log),
  data = gutdata2,
  REML = FALSE
)

plot(resid(Cmod1, type = 'pearson') ~ fitted(Cmod1))

overdisp_fun(Cmod1)  # overdispersed

Cmod2 <- update(Cmod1, 
                ziformula = ~. + offset(log(tissue.weight)))  # Add ZI structure

plot(resid(Cmod2) ~ fitted(Cmod2))  # variance still looks a bit weird

AICc(Cmod1, Cmod2)  # better fit with ZI

Cmod3 <- update(Cmod2, family = nbinom1())  # try NB error structure

AICc(Cmod2, Cmod3)  # Cmod3 is a better fit

plot(resid(Cmod3) ~ fitted(Cmod3))  # seem to be some outliers here

## test for outliers

library(DHARMa)

res <- simulateResiduals(Cmod3)
plot(res, asFactor = T)  # can keep the outliers in

## Try dropping terms

drop1(Cmod3)  # can drop trophic.position:site:particle.type

Cmod4 <- update(Cmod3, . ~ . -trophic.position:site:particle.type)

AICc(Cmod3, Cmod4)  # Cmod4 is better fit

drop1(Cmod4)  # can drop site:particle.type

Cmod5 <- update(Cmod4, . ~ . -site:particle.type)

AICc(Cmod4, Cmod5)  # Cmod5 is better fit

drop1(Cmod5)  # can drop trophic.position:particle.type

Cmod6 <- update(Cmod5, . ~ . -trophic.position:particle.type)

AICc(Cmod5, Cmod6)  # Cmod6 is a better fit

drop1(Cmod6)  # stop here


## Diagnostics

plot(resid(Cmod6) ~ fitted(Cmod6))
plot(resid(Cmod6) ~ gutdata2$trophic.position)
plot(resid(Cmod6) ~ gutdata2$site)
plot(resid(Cmod6) ~ gutdata2$particle.type)
plot(resid(Cmod6) ~ gutdata2$sample.type)

res <- simulateResiduals(Cmod6)
plot(res, asFactor = T)  # significant outlier

gutdata2$residuals <- resid(Cmod6)
plot(gutdata2$residuals)

subset(gutdata2, residuals < -15)  # remove CBMU1

gutdata3 <- subset(gutdata2, ID != CBMU1)

# Now go through  the modeling process again without the outlier

Cmod7 <-
  glmmTMB(
    count ~ trophic.position * site * particle.type +
      (1 | sample.type) + offset(log(tissue.weight)),
    family = poisson(link = log),
    data = gutdata3,
    REML = FALSE
  )

overdisp_fun(Cmod7)  # still overdispersed

Cmod8 <- update(Cmod7, ziformula = ~ . + offset(log(tissue.weight)))

AICc(Cmod7, Cmod8)  # Better fit with ZI

Cmod9 <- update(Cmod8, family = nbinom1())  # NB model won't converge

drop1(Cmod8)  # can drop the 3-way interaction

Cmod9 <- update(Cmod8, . ~ . -trophic.position:site:particle.type)

drop1(Cmod9)  # can drop trophic.position:particle.type

Cmod10 <- update(Cmod9, . ~ . -trophic.position:particle.type)

drop1(Cmod10)  # stop here

Cmod11 <- update(Cmod10, REML = TRUE)

Cmod11$call

plot(resid(Cmod11) ~ fitted(Cmod11))
plot(resid(Cmod11) ~ gutdata3$trophic.position)
plot(resid(Cmod11) ~ gutdata3$site)
plot(resid(Cmod11) ~ gutdata3$particle.type)
plot(resid(Cmod11) ~ gutdata3$sample.type)

res <- simulateResiduals(Cmod11)
plot(res, asFactor = T)  # good


## Inference

tidy(Cmod11)
summary(Cmod11)

plot_model(Cmod11)

Cmod.test1 <- update(Cmod11, . ~ . -site:particle.type)
Cmod.test2 <- update(Cmod11, . ~ . -trophic.position:site)

anova(Cmod11, Cmod.test1)  # site:particle.type significant, p<0.001
anova(Cmod11, Cmod.test2)  # trophic.position:site signficant, p=0.005

## Predictions

Cmod.pred <- gutdata3

preds2 <- predict(Cmod10, type = 'conditional', re.form = NA, se.fit = TRUE)

Cmod.pred$fitted <- preds2$fit

Cmod.pred$lower <- with(Cmod.pred, fitted - 1.96*preds2$se.fit)
Cmod.pred$upper <- with(Cmod.pred, fitted + 1.96*preds2$se.fit)  
# 95% CI for preds
# Note that predictions are for averaged REs and ZI

Cmod.pred$lower[Cmod.pred$lower < 0] <- 0


## Model liver data

liverdata <- 
  foodweb4 %>% 
  filter(!is.na(trophic.position) & particle.type != 'Natural' &
           particle.type != 'Semi-synthetic') %>%
  filter(sample.type == 'Flatfish Livers' |
           sample.type == 'Surfperch Livers' |
           sample.type == 'Rockfish Livers')

liverdata$sample.type <- as.character(liverdata$sample.type)
liverdata$sample.type <- as.factor(liverdata$sample.type)

liverdata$particle.type <- as.character(liverdata$particle.type)
liverdata$particle.type <- as.factor(liverdata$particle.type)

liverdata$species <- as.character(liverdata$species)
liverdata$species <- as.factor(liverdata$species)

Lmod1 <- glmmTMB(
  count ~ trophic.position * site * particle.type +
    (1 | sample.type),
  family = poisson(link = log),
  data = liverdata,
  REML = FALSE
)

plot(resid(Lmod1, type = 'pearson') ~ fitted(Lmod1))
overdisp_fun(Lmod1)  # not overdispersed

res <- simulateResiduals(Lmod1)
plot(res, asFactor = T)

summary(Lmod1)

Lmod2 <- update(Lmod1, family = nbinom1())  # convergence issues, NB not good

## removing terms

drop1(Lmod1)  # can remove 3-way interaction

Lmod2 <- update(Lmod1, . ~ . -trophic.position:site:particle.type)

AICc(Lmod1, Lmod2)  # Lmod2 better fit

drop1(Lmod2)  # try removing site:particle.type

Lmod3 <- update(Lmod2, . ~ . -site:particle.type)

AICc(Lmod2, Lmod3)  # Lmod 3 is better fit

drop1(Lmod3)  # try removing trophic.position:site

Lmod4 <- update(Lmod3, . ~ . -trophic.position:site)

AICc(Lmod3, Lmod4)  # Lmod4 is better fit

drop1(Lmod4)  # can remove site

Lmod5 <- update(Lmod4, . ~ . -site)

AICc(Lmod4, Lmod5)  # Lmod5 is better fit

drop1(Lmod5)  # stop here

Lmod6 <- update(Lmod5, REML = TRUE)  # refit with REML

summary(Lmod6)


## Diagnostics

plot(resid(Lmod6) ~ fitted(Lmod6))
plot(resid(Lmod6) ~ liverdata$trophic.position)
plot(resid(Lmod6) ~ liverdata$site)
plot(resid(Lmod6) ~ liverdata$particle.type)
plot(resid(Lmod6) ~ liverdata$sample.type)


## Inference

tidy(Lmod6)  # trophic level not significant, p = 0.09


## Predictions

Lmod.pred <- liverdata

preds3 <- predict(Lmod6, type = 'response', re.form = NA, se.fit = TRUE)

Lmod.pred$fitted <- preds3$fit

Lmod.pred$lower <- with(Lmod.pred, fitted - 1.96*preds3$se.fit)
Lmod.pred$upper <- with(Lmod.pred, fitted + 1.96*preds3$se.fit)  
# 95% CI for preds
# Note that predictions are for averaged REs

Lmod.pred$lower[Lmod.pred$lower < 0] <- 0  # Make sure 95% CI doesn't go under 0



## Liver concentration in terms of dry tissue weight

LCmod1 <- glmmTMB(
  count ~ trophic.position * site * particle.type +
    (1 | sample.type) + offset(log(tissue.dry.weight)),
  family = poisson(link = log),
  data = liverdata,
  REML = FALSE
)

plot(resid(LCmod1, type = 'pearson') ~ fitted(LCmod1))
overdisp_fun(LCmod1)  # overdispersed

LCmod2 <- update(LCmod1, 
                 ziformula = ~ . + offset(log(tissue.dry.weight)))
# model won't converge

res <- simulateResiduals(LCmod1)
plot(res, asFactor = T)  # significant skew in residuals

plot(resid(LCmod1, type = 'pearson') ~ liverdata$species)
plot(resid(LCmod1, type = 'pearson') ~ liverdata$total.body.wet.weight)
plot(resid(LCmod1, type = 'pearson') ~ liverdata$TL)  
# body size and species may help explain skew

LCmod2 <- update(LCmod1, . ~ . + total.body.wet.weight + species)  
# add total length and species to model

plot(resid(LCmod2, type = 'pearson') ~ fitted(LCmod2))
overdisp_fun(LCmod2)  # still overdispersed
AICc(LCmod1, LCmod2)  # better fit with body length and species

LCmod3 <- update(LCmod2, family = nbinom2())  # try NB

AICc(LCmod2, LCmod3)  # Poisson model is better fit

summary(LCmod2)

drop1(LCmod2)  # drop trophic.position:site:particle.type

LCmod3 <- update(LCmod2, . ~ . -trophic.position:site:particle.type)

AICc(LCmod2, LCmod3)  # LCmod3 is a better fit

drop1(LCmod3)  # drop trophic.position:site

LCmod4 <- update(LCmod3, . ~ . -trophic.position:site)

AICc(LCmod3, LCmod4) # LCmod4 is a better fit

drop1(LCmod4)  # can drop site:particle.type

LCmod5 <- update(LCmod4, . ~ . -site:particle.type)

AICc(LCmod4, LCmod5)  # LCmod5 is a better fit

drop1(LCmod5)  # stop here


## Diagnostics

plot(resid(LCmod5) ~ fitted(LCmod5))
plot(resid(LCmod5) ~ liverdata$trophic.position)
plot(resid(LCmod5) ~ liverdata$site)
plot(resid(LCmod5) ~ liverdata$particle.type)
plot(resid(LCmod5) ~ liverdata$sample.type)  # residual plots looks good


## Inference

LCmod5$call
summary(LCmod5)
plot_model(LCmod5)


## Prediction

LCmod.pred <- liverdata

preds4 <- predict(LCmod3, type = 'link', se.fit = TRUE)

LCmod.pred$fitted <- preds4$fit

LCmod.pred$lower <- with(LCmod.pred, fitted - 1.96*preds3$se.fit)
LCmod.pred$upper <- with(LCmod.pred, fitted + 1.96*preds3$se.fit)  
# 95% CI for preds
# Note that predictions are for averaged REs




#### Plot ####

## Gut concentration in terms of individual against trophic position

tiff('Trophic Position MP Plot.tiff',
     res = 300,
     width = 16,
     height = 12,
     units = 'cm',
     pointsize = 12)

ggplot(Imod.pred) +
  geom_ribbon(aes(x = trophic.position,
                  ymax = upper,
                  ymin = lower),
              fill = 'turquoise',
              alpha = 0.3,
              size = 0.5) +
  geom_line(aes(x = trophic.position,
                y = fitted),
            linetype = 'dashed', 
            size = 0.5) +
  geom_point(aes(x = trophic.position,
                 y = count,
                 colour = reorder(sample.type, trophic.position, mean)),
             size = 1, shape = 1, alpha = 0.8) +
  facet_grid(particle.type ~ site) +
  labs(x = 'Trophic Position',
       y = expression(paste('Particles '*ind^-1))) +
  scale_colour_manual(values = qualitative_hcl(palette = 'Dark3', n = 8)) +
  theme1

dev.off()



## Gut concentration in terms of g dry weight of gut

tiff('Trophic Position MP Weight Plot.tiff',
     res = 300,
     width = 16.7,
     height = 11,
     units = 'cm',
     pointsize = 12)

ggplot(Cmod.pred) +
  geom_ribbon(aes(x = trophic.position,
                  ymax = upper/tissue.weight,
                  ymin = lower/tissue.weight),
              fill = 'purple',
              alpha = 0.3,
              size = 0.5) +
  geom_line(aes(x = trophic.position,
                y = fitted/tissue.weight),
            linetype = 'dashed', 
            size = 0.5) +
  geom_point(aes(x = trophic.position,
                 y = count/tissue.weight,
                 colour = reorder(sample.type, trophic.position, mean)),
             size = 1, shape = 1, alpha = 1) +
  scale_y_continuous(trans = 'log1p', 
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  facet_grid(particle.type ~ site,
             labeller = label_wrap_gen(width = 9)) +
  labs(x = 'Trophic Position',
       y = expression(paste('Particles '*g^-1))) +
  scale_colour_manual(values = qualitative_hcl(palette = 'Dark3', 
                                               n = 8, rev = TRUE)) +
  theme1

dev.off()



## Liver concentration in terms of particles/individual

tiff('Trophic Position MP Liver Plot.tiff',
     res = 300,
     width = 16,
     height = 8,
     units = 'cm',
     pointsize = 12)

ggplot(Lmod.pred) +
  geom_ribbon(aes(x = trophic.position,
                  ymax = upper,
                  ymin = lower),
              fill = 'yellow',
              alpha = 0.3,
              size = 0.5) +
  geom_line(aes(x = trophic.position,
                y = fitted),
            linetype = 'dashed', 
            size = 0.5) +
  geom_point(aes(x = trophic.position,
                 y = count,
                 colour = reorder(sample.type, trophic.position, mean)),
             size = 1, shape = 1, alpha = 1) +
  facet_grid(particle.type ~ site,
             labeller = label_wrap_gen(width = 9)) +
  labs(x = 'Trophic Position',
       y = expression(paste('Particles '*ind^-1))) +
  scale_colour_manual(values = qualitative_hcl(palette = 'Dark3', n = 3)) +
  theme1

dev.off()


## Liver concentration in terms of dry tissue weight

tiff('Trophic Position MP Liver Weight Plot.tiff',
     res = 300,
     width = 16,
     height = 8,
     units = 'cm',
     pointsize = 12)

ggplot(LCmod.pred) +
  geom_ribbon(aes(x = trophic.position,
                  ymax = exp(upper)/tissue.dry.weight,
                  ymin = exp(lower)/tissue.dry.weight,
                  fill = species),
              alpha = 0.3,
              size = 0.5) +
  geom_line(aes(x = trophic.position,
                y = exp(fitted)/tissue.dry.weight,
                colour = species),
            linetype = 'dashed', 
            size = 0.5) +
  geom_point(aes(x = trophic.position,
                 y = count/tissue.dry.weight,
                 colour = reorder(species, trophic.position, mean)),
             size = 1, shape = 1, alpha = 1) +
  scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 10, 100, 1000)) +
  facet_grid(particle.type ~ site,
             labeller = label_wrap_gen(width = 9)) +
  labs(x = 'Trophic Position',
       y = expression(paste('Particles '*ind^-1))) +
  scale_colour_manual(values = qualitative_hcl(palette = 'Dark3', n = 5)) +
  scale_fill_manual(values = qualitative_hcl(palette = 'Dark3', n = 5)) +
  theme1

dev.off()



## Plot Plankton Jars


A <-
  ggplot(subset(PJ_data4, particle.type != 'Natural')) +
  geom_boxplot(aes(x = site,
                   y = count),
               size = 0.5,
               fill = 'purple',
               alpha = 0.5,
               outlier.size = 0.5) +
  facet_grid(particle.type ~ ., scales = 'free',
             space = 'free_x',
             labeller = label_wrap_gen(width = 10)) +
  labs(x = '',
       y = expression(paste('Particles ' ~ L ^ -1))) +
  theme1 +
  theme(strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1))


## Plot plankton tows

B <-
  ggplot(subset(PT_data4, particle.type != 'Natural')) +
  geom_boxplot(aes(x = site,
                   y = count/sample.volume),
               size = 0.5,
               fill = 'orange',
               alpha = 0.5,
               outlier.size = 0.5) +
  facet_grid(particle.type ~ ., scales = 'free',
             space = 'free_x',
             labeller = label_wrap_gen(width = 10)) +
  labs(x = '',
       y = expression(paste('Particles '~L^-1))) +
  theme1 +
  theme(strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1))


## Plot animal gut data

gutdata$sample.type <- factor(gutdata$sample.type, 
                              levels = c('Mussels',
                                         'Clams',
                                         'Sea Cucumbers',
                                         'Crabs',
                                         'Sea Stars',
                                         'Flatfish',
                                         'Surfperch',
                                         'Rockfish'))

C <-
  ggplot(subset(gutdata, particle.type != 'Natural')) +
  geom_boxplot(
    aes(x = sample.type,
        y = count),
    size = 0.5,
    fill = 'blue',
    alpha = 0.5,
    outlier.size = 0.5
  ) +
  facet_grid(particle.type ~ site, scales = 'free',
             space = 'free_x',
             labeller = label_wrap_gen(width = 10)) +
  labs(x = 'Type of Animal',
       y = expression(paste('Particles '~ind^-1))) +
  theme1 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 7))

water_plots <- plot_grid(A, B, labels = c('A', 'B'), nrow = 1,
                         align = 'v', rel_heights = c(1,1))

tiff('All Samples Plot.tiff',
     res = 300,
     width = 16.7,
     height = 10,
     units = 'cm',
     pointsize = 12)

plot_grid(water_plots, C, labels = c('', 'C'), nrow = 1,
          rel_widths = c(1, 1.4))

dev.off()

