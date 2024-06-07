###Checking the trait data###
library(tidyverse)
library(tidylog)
library(ggplot2)

FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_graz_conserved.csv", row.names = 1)

###Histograms of all the traits
##MAxH, MAxLS, MeanLA, MeanLL
histograms1 <- FT |> filter(trait %in% c("MaxH", "MaxLS", "MeanLA", "MeanLL")) |>
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~trait, scales = "free_x") +
  theme_classic()

##MeanLDMC, MeanLSA, percentC, percentN
histograms2 <- FT |> filter(trait %in% c("MeanLDMC", "MeanSLA", "percentC", "percentN")) |>
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~trait, scales = "free_x") +
  theme_classic()

##LDMC against SLA
FT_wide <- FT |> 
  group_by(ID, taxon, trait) |> 
  #get mean trait per sp per plot
  mutate(mean_value = mean(value)) |> 
  ungroup() |> 
  select(!value) |> 
  #make wide
  pivot_wider(names_from = trait, values_from = mean_value)
  
scatter1 <- ggplot(FT_wide, aes(x = MeanSLA, y = MeanLDMC)) +
              geom_point()+
              theme_classic()

##Height against LS
scatter2 <- 
  ggplot(FT_wide , aes(x = MaxH, y = MaxLS)) +
  geom_point()+
  theme_classic()

highest_maxh <- FT_wide |> 
  filter(MaxH == max(MaxH, na.rm = TRUE)) 
highest_maxh$taxon
#highest plant is an oak tree, makes sense

highest_LS <- FT_wide |> 
  filter(MaxLS == max(MaxLS, na.rm = TRUE)) 
highest_LS$taxon
#highest is from kameeldoring, makes sense
  

#percentC against percentN
scatter3 <-  
  ggplot(FT_wide, aes(x = percentC, y = percentN)) +
  geom_point()+
  theme_classic()

##LA against SLA
scatter4 <-  
  ggplot(FT_wide, aes(x = MeanLA, y = MeanSLA)) +
  geom_point()+
  theme_classic()
##mm there is one very big leaf

biggest_leaf <- FT_wide |> 
  filter(MeanLA == max(MeanLA, na.rm = TRUE))
biggest_leaf$taxon
#It is a palm, makes sense that the leaf is big

biggest_SLA <-  FT_wide |> 
  filter(MeanSLA == max(MeanSLA, na.rm = TRUE))
biggest_SLA$taxon

##LL against LA
scatter5 <-  
  ggplot(FT_wide, aes(x = MeanLA, y = MeanLL)) +
  geom_point()+
  theme_classic()

longest_leaf <- FT_wide |> 
  filter(MeanLL > max(MeanLL, na.rm = TRUE) - 20) |> 
  select(ID, taxon, MeanLL)
unique(longest_leaf$taxon)
longest_leaf$MeanLL
#The highest LL is from ephedra foeminea
#This plant has very small leaves and photosynthetic stems https://www.conifers.org/ep/Ephedra_foeminea.php
#could this be a mistake??

EF <- FT_wide |>  
  filter(taxon == "Ephedra foeminea")
#There are no LA records for Ephedra foeminea
#Thus trait data is not complete and it is thrown out anyway in subsequent steps. 