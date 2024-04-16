###Do nurse traits affect interaction outcomes??###
###Nintc ~ nurse traits###
library(tidyverse)
library(tidylog)
library(vegan)
library(DescTools)

#import nint results
nint_result <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names = 1) |> 
  #transform the nint variables to range between 0 and 1
  mutate(NIntc_richness_binom = (NIntc_richness + 1)/2, 
         NIntc_cover_binom = (NIntc_cover + 1)/2, 
         NInta_richness_binom = (NInta_richness - (-1)) / (2 - (-1)), 
         NInta_cover_binom = (NInta_cover - (-1)) / (2 - (-1)))
  

##import trait data
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species.csv",
               row.names = 1) |> 
  #standardise trait values
  group_by(trait) |> 
  mutate(sd_value = sd(value), 
         mean_value = mean(value)) |> 
  ungroup() |> 
  mutate(value_std = (value - mean_value)/sd_value)

##For each replicate, we need to get the NIntc value and the traits of the nurse
#names of the traits collectd
traits_collected <- c(unique(FT$trait))
#plot Id's
IDlist <- c(unique(nint_result$ID))
#empty table with explanatory variables
modeldat <- data.frame(ID = NA, site_ID = NA, replicate_no = NA, nurse_sp = NA, NIntc_richness = NA, NIntc_cover = NA, NInta_richness = NA, 
                      NInta_cover = NA, aridity = NA, graz = NA, nurse_mean_percentN = NA, nurse_mean_percentC = NA, 
                      nurse_meanLL = NA, nurse_meanSLA = NA, nurse_meanLDMC = NA, nurse_meanLA = NA, nurse_mean_H = NA, 
                      nurse_mean_LS = NA)
#columns with trait values in modeldat
trait_mean_names <- colnames(modeldat)[which(colnames(modeldat) %like% "nurse%")]

###Loop starts here
l = 1
for (i in 1:length(IDlist)) {
  
  #isolate a plot in the facilitation and FT data
  fac_plot <- nint_result[which(nint_result$ID == IDlist[i]) , ]
  FT_plot <- FT[which(FT$ID == IDlist[i]) , ]
  
  #list of reps in this plot
  replist <- c(unique(fac_plot$Number.of.replicate))
  
  #for each rep:
  for (r in 1:length(replist)) {
    one_rep <- fac_plot[which(fac_plot$Number.of.replicate == replist[r]) , ]
    #get the nurse species
    nurse_sp <- one_rep$nurse
    
    #get the trait values of the nurse species
    for (t in 1:length(traits_collected)) {
      val <- FT_plot |> 
        filter(taxon == nurse_sp, 
               trait == traits_collected[t]) |> 
        select(value_std)
      
      #if there are no values for this trait, put NA in the matrix
      if(nrow(val) == 0) {modeldat[l, which(colnames(modeldat) == trait_mean_names[t])] <- NA} 
      #if there are 2 or more trait measurements, get the mean and put that in the matrix
      else if (nrow(val) > 1) {mean_val <- mean(val$value_std)
      modeldat[l, which(colnames(modeldat) == trait_mean_names[t])] <- mean_val} 
      #if there is only on trait value, put that in the matrix
      else {modeldat[l, which(colnames(modeldat) == trait_mean_names[t])] <- val$value_std}
    }#loop through traits end
    
    #fill the rest of the table
    modeldat[l,1] <- comm_index
    modeldat[l,2] <- nurse_sp
    modeldat[l,3] <- one_rep$ARIDITY.v3[1]
    modeldat[l,4] <- one_rep$GRAZ[1]
    modeldat[l,13] <- one_rep$SITE_ID[1]
    
    l = l+1
  }#loop through reps end
}#loop through plots end
