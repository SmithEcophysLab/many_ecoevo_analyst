# euc_data_analysis.R
# analysis of Eucalyptus restoration data
# for the Many EcoEvo analysts project
# big question: How does grass cover influence Eucalyptus seedling recruitment?
# data summary: field observation data from multiple seasons of plant community composition
# data also includes a number of ancillary data
# using Zuur (2009) and https://stats.idre.ucla.edu/r/dae/zip/

## load libraries
library(tidyverse)
library(car)
# library(pscl)
library(glmmTMB) # can incorporate random effects terms
# library(lmtest)

## load data
data = read.csv('../data/Euc_data.csv')

## explore data
head(data)
colnames(data)
levels(data$Date) # 25 dates October 2006 to September 2007
levels(data$Season) # 3 seasons
levels(data$Property) # 18 properties

### data we need to know to answer the quuestion
hist(data$Euc_canopy_cover) # canopy cover of eucs, lots of zeroes!
hist(data$euc_sdlgs0_50cm) # number of small seedlings
hist(data$euc_sdlgs50cm.2m) # number of medium seedlings
hist(data$euc_sdlgs.2m) # big seedlings
data$euc_total = data$euc_sdlgs0_50cm + data$euc_sdlgs50cm.2m + data$euc_sdlgs.2m # number of seedlings total
hist(data$euc_total) # number of seedlings total
data$grass_cover = data$ExoticAnnualGrass_cover + data$ExoticPerennialGrass_cover + # total grass cover
  data$NativePerennialGrass_cover
hist(data$grass_cover) # total grass cover

### data likely to incluence Euc cover that we might need to account for
data$total_cover = rowSums(data[,10:24]) # cover of everthing
hist(data$total_cover) # cover of everthing, slightly skewed distribution with center around 100
hist(data$Distance_to_Eucalypt_canopy.m.) # distance to nearest euc canopy (some REALLY close, < 5 m)
hist(data$annual_precipitation) # annual precip >> big range
hist(data$PET) ## potential ET, not as big of a range
data$AI = data$annual_precipitation / data$PET # aridity index
hist(data$AI) # very wide range (30-90%)
data$grass_frac = data$grass_cover / data$total_cover # fraction of grass cover
hist(data$grass_frac)

## need a variable to account for resampling
data$repeated = paste(data$Property, data$Quadrat.no, sep ='')

## some important points
# what is "recruitment"? more the number of seedlings than total cover
# what is a grass? a grass like graminoid =/= a grass
# likely an interplay between climate and grass cover that needs to be teased apart
# count data Zuur pp 216
# lots of zeroes!! Zuur pp 278
# resampling!

## analyze the data (focus on counts)

### simple glm
plot(data$euc_total ~ data$grass_frac)
glm1 = glm(euc_total ~ grass_frac, family = poisson, data = data)
plot(resid(glm1) ~ fitted(glm1)) # not great!
summary(glm1)
Anova(glm1)

### zero inflated glm
#### simple model
zip1 = glmmTMB(euc_total ~ grass_frac + (1|repeated), 
               zi = ~ grass_frac,
               family = poisson,
               data = data)
plot(resid(zip1) ~ fitted(zip1)) # okay
summary(zip1) # use poisson, strong predictor of grass fraction (negative)
Anova(zip1)
AIC(zip1)

#### is it different from null model? from https://stats.idre.ucla.edu/r/dae/zip/, df = number of predictors
zip_null = glmmTMB(euc_total ~ 1 + (1|repeated), 
                   zi = ~ 1,
                   family = poisson,
                   data = data)
AIC(zip_null)
pchisq(2 * (logLik(zip1) - logLik(zip_null)), df = 1, lower.tail = FALSE) # yes

#### more complex model accounting for site aridity
zip2 = glmmTMB(euc_total ~ grass_frac + AI + (1|repeated), 
               zi = ~ grass_frac + AI,
               family = poisson,
               data = data)
plot(resid(zip2) ~ fitted(zip2)) # okay
summary(zip2) # use poisson, grass fraction still importatnt, but not AI
Anova(zip2)
AIC(zip2) # best so far

#### is it different from grass-only model? from https://stats.idre.ucla.edu/r/dae/zip/, df = number of predictors
pchisq(2 * (logLik(zip2) - logLik(zip1)), df = 2, lower.tail = FALSE) # yes

#### more complex model accounting for distance to euc forest
zip3 = glmmTMB(euc_total ~ grass_frac + Distance_to_Eucalypt_canopy.m. + (1|repeated), # won't converge! >> distance is same as random variable
               zi = ~ grass_frac + Distance_to_Eucalypt_canopy.m.,
               family = poisson,
               data = data)


## summarise across sampling dates
data_group_by_quadrat = group_by(data, repeated)
data_quadrat_mean = summarise(data_group_by_quadrat, 
                              euc_total_mean = mean(euc_total),
                              grass_frac_mean = mean(grass_frac),
                              AI_mean = mean(AI),
                              Distance_to_Eucalypt_canopy.m._mean = mean(Distance_to_Eucalypt_canopy.m.),
                              Euc_canopy_cover_mean = mean(Euc_canopy_cover))
hist(data_quadrat_mean$euc_total_mean)
hist(data_quadrat_mean$Euc_canopy_cover_mean)

## analyze
#### simple model
zip4 = glmmTMB(euc_total_mean ~ grass_frac_mean, 
               zi = ~ grass_frac_mean,
               data = data_quadrat_mean)
plot(resid(zip4) ~ fitted(zip4)) # okay
summary(zip4) # use poisson, strong predictor of grass fraction (negative)
Anova(zip4)
AIC(zip4)

#### is it different from null model? from https://stats.idre.ucla.edu/r/dae/zip/, df = number of predictors
zip_null2 = glmmTMB(euc_total_mean ~ 1, 
                   zi = ~ 1,
                   data = data_quadrat_mean)
AIC(zip_null)
pchisq(2 * (logLik(zip4) - logLik(zip_null2)), df = 1, lower.tail = FALSE) # no!

#### simple model + AI
zip5 = glmmTMB(euc_total_mean ~ grass_frac_mean + AI_mean, 
               zi = ~ grass_frac_mean,
               data = data_quadrat_mean)
plot(resid(zip5) ~ fitted(zip5)) # okay
summary(zip5) # use poisson, strong predictor of grass fraction (negative)
Anova(zip5)
AIC(zip5) # better than zip4

#### simple model + AI + Distance_to_Eucalypt_canopy.m._mean
zip6 = glmmTMB(euc_total_mean ~ grass_frac_mean + AI_mean + Distance_to_Eucalypt_canopy.m._mean, 
               zi = ~ grass_frac_mean,
               data = data_quadrat_mean)
plot(resid(zip6) ~ fitted(zip6)) # okay
summary(zip6) # use poisson, strong predictor of grass fraction (negative)
Anova(zip6)
AIC(zip6) # veryslightly better than zip5

## same analyses with Euc_canopy_cover_mean
zip7 = glmmTMB(Euc_canopy_cover_mean ~ grass_frac_mean, 
               zi = ~ grass_frac_mean,
               data = data_quadrat_mean)
plot(resid(zip7) ~ fitted(zip7)) # okay
summary(zip7) # use poisson, strong predictor of grass fraction (negative)
Anova(zip7)
AIC(zip7)

#### is it different from null model? from https://stats.idre.ucla.edu/r/dae/zip/, df = number of predictors
zip_null3 = glmmTMB(Euc_canopy_cover_mean ~ 1, 
                    zi = ~ 1,
                    data = data_quadrat_mean)
AIC(zip_null)
pchisq(2 * (logLik(zip7) - logLik(zip_null3)), df = 1, lower.tail = FALSE) # yes

#### simple model + AI
zip8 = glmmTMB(Euc_canopy_cover_mean ~ grass_frac_mean + AI_mean, 
               zi = ~ grass_frac_mean + AI_mean,
               data = data_quadrat_mean)
plot(resid(zip8) ~ fitted(zip8)) # okay
summary(zip8) # use poisson, strong predictor of grass fraction (negative)
Anova(zip8)
AIC(zip8) # better than zip7

#### simple model + AI + Distance_to_Eucalypt_canopy.m._mean
zip9 = glmmTMB(Euc_canopy_cover_mean ~ grass_frac_mean + AI_mean + Distance_to_Eucalypt_canopy.m._mean, 
               zi = ~ grass_frac_mean + AI_mean + Distance_to_Eucalypt_canopy.m._mean,
               data = data_quadrat_mean)
plot(resid(zip9) ~ fitted(zip9)) # okay
summary(zip9) # use poisson, strong predictor of grass fraction (negative)
Anova(zip9)
AIC(zip9) # veryslightly better than zip8


## simple plots
plot(Euc_canopy_cover_mean ~ grass_frac_mean, data = data_quadrat_mean)
plot(euc_total_mean ~ grass_frac_mean, data = data_quadrat_mean)


