###Compare trait values of dominant, bare associated and dominant associated species
library(tidyverse)
library(tidylog)
library(glmmTMB)
library(car)
library(DHARMa)

#import data 
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species.csv", row.names = 1) |> 
  pivot_wider(names_from = trait, values_from = value) |> 
  #calculate the C:N ratio
  mutate(C_N_ratio = percentC/percentN) |> 
  select(!c(percentC, percentN)) |> 
  pivot_longer(cols = c("MeanLL","MeanSLA","MeanLDMC","MeanLA","MaxH","MaxLS","C_N_ratio"),
               names_to = "trait", values_to = "value")

##Lets join the results of the CHi2 tests to sla_fdist###
ass <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\Chisq_results_6Feb2024.csv", row.names = 1) |> 
  select(ID, species, association) |> 
  rename(taxon= species)


FT_ass_join <- FT |> 
  left_join(ass, by = c("ID", "taxon")) |>  #ass is missing species in FT, because ass is missing the dominant species
  filter(!association %in% c("neutral", "too_rare")) |> 
  mutate(association = case_when(association == "nurse" ~ "nurse_associated", 
            association == "bare" ~ "bare_associated", 
            .default = association)) 

##Now we need to add values indicitaing which species are the dominant species 
#import the species position data
sp_positions <- read.csv("Functional trait data//Clean data//sp_positions.csv", row.names = 1) |> 
  select(!replicate) |> 
  distinct(ID, taxon, position) |> 
  filter(position == "nurse_species")

#overwrite the association column in FT_ass_join with "nurse_species" if it is a nurse in that plot
IDlist <- c(unique(sp_positions$ID))

for(i in 1:length(IDlist)) {
  
  sp_positions_plot <- sp_positions[which(sp_positions$ID == IDlist[i]) , ]
  nurse_taxa <- c(sp_positions_plot$taxon)
  
  for(t in 1:length(nurse_taxa)) {
    FT_ass_join[which(FT_ass_join$ID == IDlist[i] & FT_ass_join$taxon == nurse_taxa[t]) 
                , which(colnames(FT_ass_join) == "association")] <- "nurse_species"
  }
}
#final cleaning of FT_ass_join
FT_ass_join <- FT_ass_join |> 
  filter(!is.na(value)) 
FT_ass_join$SITE_ID <- as.factor(FT_ass_join$SITE_ID)  
FT_ass_join$association <- as.factor(FT_ass_join$association)

###Models of trait values ~ association####
#MaxH#
maxh_data <- FT_ass_join |> 
  filter(trait == "MaxH") |> 
  mutate(sqrt_value = sqrt(value), 
         log_value = log(value))

maxh_null <- glmmTMB(log_value ~ 1+ (1|SITE_ID), data = maxh_data)

maxh_mod <- glmmTMB(log_value ~ association + (1|SITE_ID), data = maxh_data)

summary(maxh_mod)
Anova(maxh_mod)
anova(maxh_null, maxh_mod) #p = 8.176e-12 ***

#significance letters
cld(glht(model = maxh_mod, mcp(association = "Tukey")))
emmeans(maxh_mod, specs = "association")

maxh_simres <- simulateResiduals(maxh_mod)
plot(maxh_simres) #with log, normality of residuals ok, HOv ok


#MaxLS#
maxls_data <- FT_ass_join |> 
  filter(trait == "MaxLS") |> 
  mutate(sqrt_value = sqrt(value), 
        log_value = log(value))

maxls_null <- glmmTMB(log_value ~ 1+ (1|SITE_ID), data = maxls_data)

maxls_mod <- glmmTMB(log_value ~ association + (1|SITE_ID), data = maxls_data)

summary(maxh_mod)
Anova(maxh_mod)
anova(maxh_null, maxh_mod) #p = 8.176e-12 ***

#significance letters
cld(glht(model = maxls_mod, mcp(association = "Tukey")))
emmeans(maxls_mod, specs = "association")

maxls_simres <- simulateResiduals(maxls_mod)
plot(maxls_simres) #with log, normality of residuals ok, HOv ok


#MeanLL#
meanll_data <- FT_ass_join |> 
  filter(trait == "MeanLL") |> 
  mutate(sqrt_value = sqrt(value), 
         log_value = log(value))

meanll_null <- glmmTMB(sqrt_value ~ 1+ (1|SITE_ID), data = meanll_data) #cannot use log transformation because LL has 0 values

meanll_mod <- glmmTMB(sqrt_value ~ association + (1|SITE_ID), data = meanll_data)

summary(meanll_mod)
Anova(meanll_mod)
anova(meanll_null, meanll_mod) #p = 0.8483


meanll_simres <- simulateResiduals(meanll_mod)
plot(meanll_simres) #a little underdispersed, HOV ok


#MeanLDMC#
meanldmc_data <- FT_ass_join |> 
  filter(trait == "MeanLDMC") |> 
  mutate(sqrt_value = sqrt(value), 
         log_value = log(value))

meanldmc_null <- glmmTMB(value ~ 1+ (1|SITE_ID), data = meanldmc_data) 

meanldmc_mod <- glmmTMB(value ~ association + (1|SITE_ID), data = meanldmc_data)

summary(meanldmc_mod)
Anova(meanldmc_mod)
anova(meanldmc_null, meanldmc_mod) #0.1325

#significance letters
cld(glht(model = meanldmc_mod, mcp(association = "Tukey")))
emmeans(meanldmc_mod, specs = "association")


meanldmc_simres <- simulateResiduals(meanldmc_mod)
plot(meanldmc_simres) #everything is good without the transformation


#MeanSLA#
meansla_data <- FT_ass_join |> 
  filter(trait == "MeanSLA") |> 
  mutate(sqrt_value = sqrt(value), 
         log_value = log(value))

meansla_null <- glmmTMB(log_value ~ 1+ (1|SITE_ID), data = meansla_data) 

meansla_mod <- glmmTMB(log_value ~ association + (1|SITE_ID), data = meansla_data)

summary(meansla_mod)
Anova(meansla_mod)
anova(meansla_null, meansla_mod) #p = 0.07813 .

meansla_simres <- simulateResiduals(meansla_mod)
plot(meansla_simres) #residuals normal, HOV good


#MeanLA#
meanla_data <- FT_ass_join |> 
  filter(trait == "MeanLA") |> 
  mutate(sqrt_value = sqrt(value), 
         log_value = log(value))

meanla_null <- glmmTMB(log_value ~ 1+ (1|SITE_ID), data = meanla_data) 

meanla_mod <- glmmTMB(log_value ~ association + (1|SITE_ID), data = meanla_data)

summary(meanla_mod)
Anova(meanla_mod)
anova(meanla_null, meanla_mod) #p = 0.4094

meanla_simres <- simulateResiduals(meanla_mod)
plot(meanla_simres) #residuals normal, HOV good


#C:N ratio#
cn_data <- FT_ass_join |> 
  filter(trait == "C_N_ratio") |> 
  mutate(sqrt_value = sqrt(value), 
         log_value = log(value))

cn_null <- glmmTMB(log_value ~ 1+ (1|SITE_ID), data = cn_data) 

cn_mod <- glmmTMB(log_value ~ association + (1|SITE_ID), data = cn_data)

summary(cn_mod)
Anova(cn_mod)
anova(cn_null, cn_mod) #p = 0.7925

cn_simres <- simulateResiduals(cn_mod)
plot(cn_simres) #residuals normal, HOV good


##Lets make some graphs####
ggplot(FT_ass_join, aes(x = association, y = value)) +
  geom_boxplot(fill = "darkslategrey", alpha = 0.6) +
  facet_wrap(~trait, scale = "free_y") +
  theme_classic()

