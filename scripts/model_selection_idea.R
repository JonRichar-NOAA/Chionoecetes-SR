## some ideas for model fitting

library(tidyverse)
library(mgcv)
library(MuMin)

## set up a data frame with all the data
dat <- data.frame(year = ,
                  recr = ,
                  spawn = ,
                  cod = ,
                  fhs = ,
                  temp = ,
                  wind = ) # or whatever parsimonious set of variables you want to use! 

## examine distribution for recruits

hist(dat$recr)

## if skewed, may need a family other than Gaussian, e.g. Gamma

## now fit the full model

mod1 <- gam(recr ~ s(spawn, k = 4) + s(spawn, k = 4) + s(cod, k = 4) + s(fhs, k = 4) + s(temp, k = 4) + s(wind, k = 4),
            data = dat)

summary(mod1)

## drop the least important covariate (the one with highest p-value)

## e.g., 
mod2 <- gam(recr ~ s(spawn, k = 4) + s(spawn, k = 4) + s(cod, k = 4) + s(fhs, k = 4) + s(temp, k = 4),
            data = dat)

AICc(mod1, mod2)

## if mod2 is better, drop the least important covariate from mod2 and run as mod3. Keep doing this until the AICc score stops improving!

## if you want to account for autocorrelation in the data, could use GAMM with a random year term, or parameterize 
## the best model in GLS with autocorrelated residuals, setting the order of polynomial terms roughly equal to edf from the GAM!