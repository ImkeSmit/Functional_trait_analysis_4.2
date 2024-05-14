###Compare trait values of dominant, bare associated and dominant associated species
library(tidyverse)
library(tidylog)

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


#Now we need to get the mean trait values of nurses and bare or nurse associated targets 
summary <- FT_ass_join |> 
  group_by(association, trait) |> 
  summarise(mean_trait_value = mean(value, na.rm = T))
