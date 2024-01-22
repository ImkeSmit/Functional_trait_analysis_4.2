###playing with the traitstrap package
library(traitstrap)
library(ggplot2)
library(tidyverse)

wd <- "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait data"
setwd(wd)

#import data
#This data contains species without trait values
trait_raw <- read.csv("FT_match_facilitation_plots_plotspecific_species.csv", row.names = 1) 
#make the data long format to use in traitstrap
trait_long <-  trait_raw |> 
  rename(taxon = sp_fullname) |> 
  pivot_longer(cols = c(MeanLL, MeanSLA, MeanLDMC, MeanLA, MaxH, MaxLS), names_to = "trait", values_to = "value") #make it long format
  


###Plot some of the traits
#plot SLA vs LDMC
ggplot(trait_dat, aes(x = MeanSLA, y = MeanLDMC)) +
  geom_point()

##LL vs LA
ggplot(trait_dat, aes(x = MeanLL, y = MeanLA)) +
  geom_point()

##MaxH vs MaxLS
ggplot(trait_dat, aes(x = MaxH, y = MaxLS)) +
  geom_point()

##SLA vs aridity
ggplot(trait_dat, aes(x = ARIDITY.v3, y = MeanSLA)) +
  geom_point()

##MaxH vs aridity
ggplot(trait_dat, aes(x = ARIDITY.v3, y = MaxH)) +
  geom_point()

##LDMC vs aridity
ggplot(trait_dat, aes(x = ARIDITY.v3, y = MeanLDMC)) +
  geom_point()



###get cover per species per plot (over all microsites in a plot) to use to weight trait metrics by abundance

#To get this data we have to import the raw data again
wd <- "C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation data\\Countriesv2"
###read in the facilitation data for each country
data_files <- list.files(wd)
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")
for(i in 1:length(data_files)) {                              
  assign(paste0(countrynames[i]),                                   
         read.csv2(paste0("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation data\\Countriesv2\\",
                          data_files[i])))
}

#rbind all the countries together, with only the variables we need
vars <- c("COU", "SITE_ID", "ID", "GRAZ", "ARIDITY.v3", "Microsite", "ID_Microsite", "Number.of.replicate", "Species.within.quadrat", "Cover" )
all_countries <- rbind(algeria[ , which(colnames(algeria) %in% vars)], argentina[ , which(colnames(argentina) %in% vars)], australia[ , which(colnames(australia) %in% vars)], chile[ , which(colnames(chile) %in% vars)], 
              chinachong[ , which(colnames(chinachong) %in% vars)], chinaxin[ , which(colnames(chinaxin) %in% vars)], iranabedi[ , which(colnames(iranabedi) %in% vars)], iranfarzam[ , which(colnames(iranfarzam) %in% vars)], 
              israel[ , which(colnames(israel) %in% vars)], namibiablaum[ , which(colnames(namibiablaum) %in% vars)], namibiawang[ , which(colnames(namibiawang) %in% vars)], southafrica[ , which(colnames(southafrica) %in% vars)],  
              spainmaestre[ , which(colnames(spainmaestre) %in% vars)], spainrey[ , which(colnames(spainrey) %in% vars)])


coversum <- all_countries |> 
  group_by(SITE_ID, ID, Species.within.quadrat) |> 
  summarise(sum_sp_cover = sum(Cover, na.rm = TRUE)) |> 
  ungroup()

denom <- all_countries |>
  group_by(ID) |> 
  summarise(denominator = (100*max(Number.of.replicate)*2)) 

cover_perplot <- left_join(coversum, denom, by = "ID") |> 
 mutate(percent_cover_perplot = sum_sp_cover/denominator) |> 
  rename(taxon = Species.within.quadrat, site_ID = SITE_ID)



#do the trait filling
trait_filling <- trait_fill(
  comm = cover_perplot, #community data with abundance values
  traits = trait_long, #trait data
  abundance_col = "percent_cover_perplot", #the column with the abundance data
  taxon_col = "taxon", #column with the speciesnames (must be the same in both datasets)
  value_col = "value",
  trait_col = "trait",
  
  # specifies sampling hierarchy
  scale_hierarchy = c("site_ID", "ID"),
  
  # min number of samples
  min_n_in_sample = 1 #if there is one trait value in the sample, do not search for values higher up in the hierarchy
)


#show where the traits were filled from
autoplot(trait_filling) +
  scale_fill_manual(values = c("red", "orange", "blue")) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5))

#import facilitation_species_and_positions, to see if trait filling is just throwing out the nurses.
#because the nurses are not present in cover_perplot (except if the nurse species also occurs in a microsite)
wd <- "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait data"
setwd(wd)
fac_sp <- read.csv("facilitation_species_and_positions.csv", row.names = 1)
head(fac_sp)

#look at plot for an example
filled_2 <- trait_filling[which(trait_filling$ID == 2) , ]
fac_sp_2 <- fac_sp[which(fac_sp$ID == 2) , ]
long_2 <- trait_long[which(trait_long$ID == 2) , ]
cover_2 <- cover_perplot[which(cover_perplot$ID == 2) , ]

unique(filled_2$taxon)
unique(long_2$taxon)
unique(fac_sp_2$spname)
unique(cover_2$taxon)

#look at the raw distributions
raw_dist_np <- trait_np_bootstrap(
  filled_traits = trait_filling,
  sample_size = 100,
  raw = TRUE)
 #I think if you say raw = TRUE it does the resampling for one rep, and doesn't calculate moments

raw_dist_np$GRAZ <- as.factor(raw_dist_np$GRAZ)
ggplot(raw_dist_np, aes(x = log(value), fill = GRAZ)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("red", "green", "blue", "orange")) +
  labs(x = "log(trait value)") +
  facet_wrap(facets = vars(trait), scales = "free")

#look at the bootstrapped means and higher moments
np_bootstrapped_moments <- trait_np_bootstrap(
  trait_filling, 
  sample_size = 100, #trait values are resampled 100 times
  nrep = 50) #create 50 distributions


#now we can get the confidence intervals
sum_boot_moment <- trait_summarise_boot_moments(
  np_bootstrapped_moments, 
  parametric = TRUE) #uses moment +- 1sd

#calculate BS multivariate metrics?