##This script is to run the second analysis to determine whther the functional distance between dominant and 
#target species are affected by target species' association and environmental factors
library(tidyverse)
library(tidylog)
library(glmmTMB)
library(MuMIn)
library(car)

#Import pairwise differences between traits
trait_fdist <- read.csv("Functional trait data\\results\\trait_differences_between_2sp_traits_vary.csv", row.names = 1) |> 
  filter(trait %in% c("MaxH", "MeanLDMC"))
trait_fdist$SITE_ID <- as.factor(trait_fdist$SITE_ID)
##Lets join the results of the CHi2 tests to sla_fdist###
ass <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\Chisq_results_6Feb2024.csv", row.names = 1) |> 
  select(ID, species, association) |> 
  rename(target = species)
#remember that these associations were calculated were calculated at the plot scale. Eg in a specific plot, a species has a significant association with nurse microsites

#import siteinfo which has lat and long for each plot
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  select(ID, Lat_decimal, Long_decimal) |> 
  mutate(sin_lat = sin(Lat_decimal), 
         sin_long = sin(Long_decimal)) |> 
  select(!c(Lat_decimal, Long_decimal))

#join the associations and the coordinates to the trait differences
trait_ass_join <- trait_fdist |> 
  left_join(ass, by = c("target", "ID")) |> 
  filter(association %in% c("nurse", "bare")) |> #only work with these associations
  left_join(siteinfo, by = "ID")
trait_ass_join$association <- as.factor(trait_ass_join$association)
trait_ass_join$nurse <- as.factor(trait_ass_join$nurse)
trait_ass_join$SITE_ID <- as.factor(trait_ass_join$SITE_ID)
trait_ass_join$ID <- as.factor(trait_ass_join$ID)

####MODEL SELECTION FOR MAXH####
#MaxH model#
maxh_data <- trait_ass_join |> 
  filter(trait == "MaxH") |> 
  mutate(sqrt_trait_difference = sqrt(trait_difference), 
         log_trait_difference = log(trait_difference), 
         neginv_trait_difference = -1/(1+trait_difference))
hist(maxh_data$neginv_trait_difference)
hist(maxh_data$trait_difference)

###full formula for model selection
full_formula <- as.formula("trait_difference ~ association*graz + 
                           association*AMT + association*RASE + association*aridity + 
                           association* SAC + association*pH +
                           sin_lat + sin_long + (1|nurse_sp)")



#null model
maxh_null <- glmmTMB(trait_difference ~ 1 + (1|nurse) + (1|SITE_ID/ID), data = maxh_data) #cannot use transformations because of zeroes and negative values
#alternative model
maxh_mod <- glmmTMB(trait_difference ~ association + (1|nurse) + (1|SITE_ID/ID), data = maxh_data)

summary(maxh_mod)
Anova(maxh_mod) #significant effect
anova(maxh_null, maxh_mod) #p = 1.113e-07 ***
emmeans(maxh_mod, specs = "association")
r.squaredGLMM(maxh_mod)

#model diagnostics
simres <- simulateResiduals(maxh_mod)
plot(simres)#underispersed, HOV violated

ggplot(maxh_data, aes(x = association, y = trait_difference)) +
  geom_boxplot()

##Are the means of each association group different from 0?
#Do Whelch's t.test, which does not assume equal variances
maxh_test_nurse <- t.test(maxh_data[which(maxh_data$association == "nurse") , ]$trait_difference, 
                          mu = 0, alternative = "two.sided")
#p <0.001, true mean greater than 0

maxh_test_bare <- t.test(maxh_data[which(maxh_data$association == "bare") , ]$trait_difference, 
                         mu = 0, alternative = "two.sided")
#p <0.001, true mean greater than 0