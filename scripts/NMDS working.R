#This script is playing around with NMDS ordination and permanova

library(vegan)
library(tidyverse)
library(tidylog)
library(ggplot)

set.seed(100)

#load data
#From the filled trait data for plotspecific species
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species.csv", row.names = 1)
##!remember that these species were not necessarily filled from the same graz level

##We need to standardise each trait value to mean 0 and variance of 1 with Z = (x-mean)/sd
##We need only one trait value per species

#Get the mean of each trait for each sp
FT_mean <- FT |> 
  group_by(taxon,trait) |> 
  summarise(mean_value = mean(value))

#Get it into wide format
FT_wide <- FT_mean |>
  pivot_wider(names_from = trait, values_from = mean_value) |> 
  column_to_rownames(var = "taxon") |> 
  filter(!is.na(MaxH),
         !is.na(MaxLS), 
         !is.na(MeanLA),
         !is.na(MeanLDMC),
         !is.na(MeanLL),
         !is.na(MeanSLA), 
         !is.na(percentC),
         !is.na(percentN)) |> 
  mutate(C_N_ratio = percentC/percentN) |> 
  select(!c(percentC, percentN))



#standardise trait values
std_FT_wide <- FT_wide

traitlist <- c("MaxH","MaxLS","MeanLA","MeanLDMC","MeanLL","MeanSLA","C_N_ratio")

for(t in 1:length(traitlist)) {
  
  #get the grand mean and sd for each trait
  m <- FT_wide |> 
    summarise(m = mean(FT_wide[ , which(colnames(FT_wide) == traitlist[t])]))
  
  sd <- FT_wide |> 
    summarise(m = sd(FT_wide[ , which(colnames(FT_wide) == traitlist[t])]))
  
  for (i in 1:nrow(FT_wide)) {
    
    raw_value <- FT_wide[i , which(colnames(FT_wide) == traitlist[t])]
    #replace the raw value with the standardised value using m and sd
    std_FT_wide[i , which(colnames(std_FT_wide) == traitlist[t])] <- (raw_value - m)/sd
    
  }
}

FT_sp <- c(rownames(std_FT_wide))


##Add the microsite affinity of each species
mic_af <- read.csv("Functional trait data\\Clean data\\sp_positions.csv", row.names = 1) |> 
  filter(!is.na(taxon), 
         !position == "nurse_species", 
         taxon %in% FT_sp) |> 
  distinct(taxon, position) |> 
  arrange(taxon)

#if a species has more than one position, change their position to "both"
#so if a species ever, in any set of plots occurs in more than one microsite, their microsite affinity is changed to "both"
for (r in 1:length(FT_sp)) {
  a_species <- mic_af[mic_af$taxon == FT_sp[r] , ]
  
  if(nrow(a_species) > 1) {
    mic_af[mic_af$taxon == FT_sp[r] , 2] <- "both"
  }
}

mic_af <- mic_af |> 
  distinct(taxon, position)
###mmm now most species have both as position, this is a problem
#is there a way to do it plot by plot?



#do NMDS
test_nmds <- metaMDS(std_FT_wide, distance = "euclidean")
#Extract NMDS Scores 
nmds_scores <-as.data.frame(scores(test_nmds))

ggplot(nmds_scores, aes(NMDS1, NMDS2)) +
  geom_point()


#Now let's get the species that occur in each grazing level, and which of those are nurses and targets
#all the species in the trait data
sp_graz <- FT |> 
  group_by(GRAZ) |> 
  distinct(taxon) |> 
  ungroup()


##add the grazing levels to std FT wide

