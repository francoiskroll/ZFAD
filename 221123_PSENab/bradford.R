#####################################################
# ~ ZFAD: amyloid beta measurements of PSEN1/2 knockouts, Bradford assay ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################



# packages ----------------------------------------------------------------

library(openxlsx)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)

library(here)



# import ------------------------------------------------------------------

bra <- read.xlsx(here('221123_PSENab', 'ABmeasurements', '230202_PS12AB_bradford.xlsx'), sheet='table')

# keep standards & sample separated
std <- bra %>%
  filter(group=='standard')

spl <- bra %>%
  filter(group!='standard')



# prepare standard curve --------------------------------------------------

# have a look how raw data looks
ggplot(std, aes(x=concentration, y=absorbance)) +
  geom_point()

# first, calculate, for each concentration of standard, the mean of the two replicate sample

# cal for calibration
cal <- std %>%
  group_by(concentration) %>% # do not group by replicate, so for each concentation the two samples get pooled
  summarise_at(vars(absorbance),
               list(
                 mean=~mean(.)
               ))

# how mean calibration curve looks like
ggplot(cal, aes(x=concentration, y=mean)) +
  geom_point()

### calculate linear regression
lmcal <- lm(mean ~ concentration, data=cal)
# absorbance = 0.0007147 * concentration + 0.8434268

# R2?
# this is simply square of correlation
cor(cal$mean, cal$concentration)^2
# 0.9061379

# on data before averages?
cor(std$absorbance, std$concentration)^2
# 0.8971334

## note, curve is more linear at lower concentrations (without the two samples with largest concentrations, 1500 & 2000)
# at 1000, largest absorbance is 1.775
max(spl$absorbance) # and max absorbance in our data is below that
# so exclude last two samples, which is closer to the range we have


### calculate linear regression v2
cal2 <- cal %>%
  filter(!concentration %in% c(1500, 2000))

lmcal2 <- lm(mean ~ concentration, data=cal2)
# absorbance = 0.001079 * concentration + 0.732858

# R2?
# this is simply square of correlation
cor(cal2$mean, cal2$concentration)^2
# 0.9647243
# so R2 indeed better without the largest concentrations

# intercept is:
int <- as.numeric(lmcal2$coefficients[1])
# slope is:
slo <- as.numeric(lmcal2$coefficients[2])


# add regression line to plot
ggplot(cal2, aes(x=concentration, y=mean)) +
  geom_point() +
  geom_smooth(method='lm') +
  coord_cartesian(ylim=c(0,3))

# or manually to check same
ggplot(cal2, aes(x=concentration, y=mean)) +
  geom_point() +
  geom_segment(x=0, y=int, xend=2000, yend=slo*2000+int) +
  coord_cartesian(ylim=c(0,3))
# yes, exactly the same


# small function to convert absorbance to concentration -------------------

abs2Conc <- function(abs,
                     slo,
                     int) {
  return( (abs-int)/slo  )
}



# calculate concentrations for samples ------------------------------------

cons <- sapply(spl$absorbance, function(ab) {
  return(abs2Conc(abs=ab, slo=slo, int=int))
})

spl <- spl %>%
  mutate(predConc=cons)


## calculate the mean of the two samples' concentrations for each group
splm <- spl %>%
  group_by(group) %>% # do not group by replicate, so for each concentation the two samples get pooled
  summarise_at(vars(predConc),
               list(
                 mean=~mean(.)
               ))
# looks like same results from Guliz



# correct concentrations --------------------------------------------------

# these are concentrations are of samples diluted 10x
# in order to be in linear range & to have enough material
# so multiply by 10 here
splm <- splm %>%
  mutate(mean=mean*10)


# deal with diluted samples -----------------------------------------------

# a solution would be to multiply their concentrations by 5 and average with concentrations of undiluted samples
# but I might *not* do this, some are relatively far off from the undiluted sample
# they are in the lower range of the standard curve, so may be imprecise


# write results -----------------------------------------------------------

colnames(splm) <- c('sample', 'concentration')

write.csv(splm, here('221123_PSENab', 'ABmeasurements', 'bradford_results.csv'), row.names=FALSE)

# units are µg/mL



# sanity check ------------------------------------------------------------

# do these numbers make sense?
# something like 3000--5000 µg/mL in ~ 20 juveniles

# rough calculations of volume of 21 dpf brain was 0.52 mm3
# dissected whole head, so will be a bit bigger here, but they are also younger than 21 dpf, so will assume OK estimate
# i.e. 0.00052 g if assuming = water
# there are ~ 20 juveniles, so 0.0104 g total
# imagine tissue is liquid, we have 0.52 mm3 * 20 juveniles = 10.4 mm3 of tissue, which is ~ 10.4 µL
# mass of tissue in the Eppendorf is 0.0104 g
# of which ~ 20% is protein, so 0.00208 g = 2,080 µg ***
# on top of which we add 100 µL lysis buffer, so ~ 110 µL total
# now in 110 µL total = 0.11 mL
# so 2080/0.11 ~ 19,000 µg/mL

# and we get ~ 3,000--5,000, so ~ 20% recovery, which seems low but reasonable

# ***
# (cf. https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=1&id=100448)
# other source (nutritional stuff): protein content of fish is ~ 20g / 100g, i.e. 20%, which confirms step above
# is it that only 2.5% of protein escapes tissue and goes in solution during lysis? this seems way too low

#####

# in fact, we also have the weight in Guliz's data
# she found ~ 0.05 g for ~ 6 brains, i.e. ~ 0.008 mg / brain
# 0.05 g = 50,000 µg
# 20% is protein, so 10,000 µg

# adult brain is ~ 3 mm3 (cf. https://elifesciences.org/articles/69988, Author Response)
# so 6 adult brains = 18 mm3, so ~ 18 µL
# add 100 µL lysis buffer, so total 118 µL in Eppendorf
# so 10,000 ug protein / 118 µL = 85 µg/µL = 85,000 µg/mL
# we found ~ 5000 µg/mL so ~ 6% recovery