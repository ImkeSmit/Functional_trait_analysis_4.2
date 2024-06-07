###TRAIT FILLING###
#This is the script for trait filling and subsetting the filled FT data for the facilitation plots##
#Weight the trait filling by CoverdryFun20
library(tidyverse)
library(tidylog)
library(traitstrap)
library(glmmTMB)
library(car)

#we will import FT data for all sites that have FT data. 
#Then do trait filling
#Then subset for facilitation plots to have most complete trait data. 

#import trait data 
FT <- read.csv("Functional trait data\\Clean data\\FT_all_sites.csv", row.names = 1) |> 
  #some values of CoverBiodesert 100 exceed 100, therefore I assume it is the sum of the cover of a species over all 100 quadrats
  #therefore divide by 100
  mutate(coverBiodesert100 = coverBiodesert100/100, 
         #make sure that there are no 0 cover values.
         #values must be >0 or NA
         coverBiodesert100 = case_when(coverBiodesert100 == 0 ~ NA, .default = as.numeric(coverBiodesert100)))
FT$GRAZ <- as.factor(FT$GRAZ)
FT$SITE_ID <- as.factor(FT$SITE_ID)


###But first: see if traits are influenced by GRAZ####
#to determine wehther we need to take GRAZ into account when filling

#meanLL~graz#
ll_dat <- FT |> 
  filter(!is.na(MeanLL)) |> 
  mutate(log_MeanLL = log(MeanLL), 
         sqrt_MeanLL = sqrt(MeanLL))

ll_null <- glmmTMB(sqrt_MeanLL ~ 1 + (1|SITE_ID), data = ll_dat)

ll_mod <- glmmTMB(sqrt_MeanLL ~ GRAZ + (1|SITE_ID), data = ll_dat)
anova(ll_null, ll_mod) #p =  0.02515 *
Anova(ll_mod)


#meanSLA~graz#
sla_dat <- FT |> 
  filter(!is.na(MeanSLA)) |> 
  mutate(log_MeanSLA = log(MeanSLA), 
         sqrt_MeanSLA = sqrt(MeanSLA)) 

sla_null <- glmmTMB(sqrt_MeanSLA ~ 1 + (1|SITE_ID), data = sla_dat)

sla_mod <- glmmTMB(sqrt_MeanSLA ~ GRAZ + (1|SITE_ID), data = sla_dat)
anova(sla_null, sla_mod)#p = 0.02515 *
Anova(sla_mod)


#meanLDMC~graz#
ldmc_dat <- FT |> 
  filter(!is.na(MeanLDMC)) |> 
  mutate(log_MeanLDMC = log(MeanLDMC), 
         sqrt_MeanLDMC = sqrt(MeanLDMC)) 

ldmc_null <- glmmTMB(sqrt_MeanLDMC ~ 1 + (1|SITE_ID), data = ldmc_dat)

ldmc_mod <- glmmTMB(sqrt_MeanLDMC ~ GRAZ + (1|SITE_ID), data = ldmc_dat)
anova(ldmc_null, ldmc_mod) #p = 0.02264 *
Anova(ldmc_mod)


#meanLA~graz#
la_dat <- FT |> 
  filter(!is.na(MeanLA)) |> 
  mutate(log_MeanLA = log(MeanLA), 
         sqrt_MeanLA = sqrt(MeanLA)) 

la_null <- glmmTMB(log_MeanLA ~ 1 + (1|SITE_ID), data = la_dat)

la_mod <- glmmTMB(log_MeanLA ~ GRAZ + (1|SITE_ID), data = la_dat)
anova(la_null, la_mod) #p = 0.4067
Anova(la_mod)


#maxH~graz#
h_dat <- FT |> 
  filter(!is.na(MaxH)) |> 
  mutate(log_MaxH = log(MaxH), 
         sqrt_MaxH = sqrt(MaxH)) 

h_null <- glmmTMB(log_MaxH ~ 1 + (1|SITE_ID), data = h_dat)

h_mod <- glmmTMB(log_MaxH ~ GRAZ + (1|SITE_ID), data = h_dat)
anova(h_null, h_mod) #p = 7.452e-05 ***
Anova(h_mod)


#maxLS~graz#
ls_dat <- FT |> 
  filter(!is.na(MaxLS)) |> 
  mutate(log_MaxLS = log(MaxLS), 
         sqrt_MaxLS = sqrt(MaxLS)) 

ls_null <- glmmTMB(log_MaxLS ~ 1 + (1|SITE_ID), data = ls_dat)

ls_mod <- glmmTMB(log_MaxLS ~ GRAZ + (1|SITE_ID), data = ls_dat)
anova(ls_null, ls_mod) #p = 0.003995 **
Anova(ls_mod)

##GRAZING MATTERS##

###Filling starts here####

#Change FT to long format
FT_long <- FT |> 
  rename(coverBiodesert100_FT = coverBiodesert100) |> #keep this column here so that we know what the cover was in the trait data
  select(!c(Genus, Species, coverDryfun20)) |> 
  pivot_longer(cols = c(MeanLL, MeanSLA, MeanLDMC, MeanLA, MaxH, MaxLS, percentN, percentC), 
                       names_to = "trait", values_to = "value")

FT_long$ID <- as.factor(FT_long$ID)
FT_long$SITE_ID <- as.factor(FT_long$SITE_ID)
FT_long$COU <- as.factor(FT_long$COU)
#traitstrap assumes the first level is the control, since we do not have a control, we make an empty first level
FT_long$GRAZ <- factor(FT_long$GRAZ, levels = c(100, 0,1,2,3))

#Create the comm table containing the cover values for each species
#Use CoverBiodesert100 because we cannot get cover values for the nurse species from the facilitation data
CovBio100 <- FT |> 
  select(c(ID, COU, SITE, SITE_ID, PLOT, taxon, coverBiodesert100, GRAZ, ARIDITY.v3)) |> 
  distinct() |> 
  filter(!is.na(coverBiodesert100)) #remove rows with NA cover values
CovBio100$ID <- as.factor(CovBio100$ID)
CovBio100$SITE_ID <- as.factor(CovBio100$SITE_ID)
CovBio100$COU <- as.factor(CovBio100$COU)
#traitstrap assumes the first level is the control, since we do not have a control, we make an empty first level
CovBio100$GRAZ <- factor(CovBio100$GRAZ, levels = c(100, 0,1,2,3)) 


FT_filled <- trait_fill(
  comm = CovBio100, #community data with abundance values
  traits = FT_long, #trait data
  abundance_col = "coverBiodesert100", #the column with the abundance data
  taxon_col = "taxon", #column with the speciesnames (must be the same in both datasets)
  value_col = "value",
  trait_col = "trait",
  
  global = FALSE, #do not calculate traits at the global scale
  
  keep_all = FALSE, #do not keep trait data at all availible levels, only on the finest scale availible
  
  # specifies sampling hierarchy
  scale_hierarchy = c("COU", "ID"), #in order to retrieve traits from only the same graz level, traits can only be filled from country.
                                    #you cannot find the same grazlevel from SITE_ID.
  #only fill with trait values from the same graz level
  treatment_col = "GRAZ",
  treatment_level = "COU", #hierarchy level at which you will find the same graz treatment level
  
  other_col = c("SITE_ID"),
  
  # min number of samples
  min_n_in_sample = 1 #if there is one trait value in the sample, do not search for values higher up in the hierarchy
)
FT_filled_sorted <- arrange(FT_filled, ID)
#make sure GRAZ_comm and GRAZ_trait are always the same
nrow(FT_filled[which(FT_filled$GRAZ_comm != FT_filled$GRAZ_comm) ,]) #GRAZ_comm is always equal to GRAZ_trait

write.csv(FT_filled, "Functional trait data\\Clean data\\FT_filled_all_sites_graz_conserved.csv")
#coverBiodesert100 is the cover value of that species in that ID
#coverBiodesert100_FT is the cover value of the "filler", i.e. it is the cover of the species entry that came from elsewhere to fill this gap


#show where the traits were filled from
autoplot(FT_filled, other_col_how = "ignore") +
  scale_fill_manual(values = c("grey","orange", "blue")) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5))

#look at missing traits
missing_traits <- trait_missing(
  filled_trait = FT_filled,
  comm = quad_summary)


###Subset FT_filled for the facilitation plots####
FT_filled <- read.csv("Functional trait data\\Clean data\\FT_filled_all_sites_graz_conserved.csv", row.names = 1)

#Import the country_v3 data
data_files <- list.files("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\Countriesv3")
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")
for(i in 1:length(data_files)) {                              
  assign(paste0(countrynames[i]),                                   
         read.csv2(paste0("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\Countriesv3\\",
                          data_files[i])))
}

#loop to get the IDs of plot in the fac data
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")

fac_IDs <- c(rep(NA, 3))

for (t in 1:length(countrynames)) {
  cou <- get(countrynames[t])
  
  if(t ==1) {
    
    fac_IDs <- c(unique(cou$ID))
    
  }else {
    
    temp_fac_IDs <- c(unique(cou$ID))
    fac_IDs <-c(fac_IDs, temp_fac_IDs)
  }
}

FT_filled_facplots <- FT_filled |> 
  filter(ID %in% c(fac_IDs)) |> 
  ungroup()
length(unique(FT_filled_facplots$ID))

#export
write.csv(FT_filled_facplots, "Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_graz_conserved.csv")



###SUBSET THE FT_filled_facplots FOR THE SPECIES IN THAT PLOT IN THE FACILITATION DATA####
FT_filled_facplots <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_graz_conserved.csv", row.names = 1)

###ALSO FIND OUT WHICH SPECIES ARE IN BOTH DATASETS OR ONLY ON ONE?
#import spnames in the facilitation data
fac_species <- read.csv("Functional trait data\\Clean data\\facilitation_species_and_positions.csv", row.names = 1)

#which sp are only in the facilitation data, but not in the trait data, this table is sp_matches
sp_matches <- data.frame(matrix(nrow = 10, ncol = 3))
colnames(sp_matches) <- c("ID", "spname", "position")
#Position = fac_only, if the sp is only in the facilitation data, FT_only if it is only in the FT data, match if the sp is present in both datasets

#We want to do this separately for each plot, because it will be a plot level analysis
fac_IDs <- c(unique(as.numeric(fac_species$ID))) #plot ID's from the facilitation data

##Loop to subset trait_subset for the facilitation sp and to get the names of species in both or only one dataset
##The table subsetted for the fac species in that plot is FT_sub
for(i in 1:length(fac_IDs)) {
  
  FT_plot <- FT_filled_facplots[which(FT_filled_facplots$ID == fac_IDs[i]) , ] #subset by plot in the trait data
  fac_sp_plot <- fac_species[which(fac_species$ID == fac_IDs[i]) , ] #subset by plot in fac_species
  
  fac_spnames <- c(unique(fac_sp_plot$spname)) #names of the species encountered in the facilitation dataset
  
  #species in both datasets
  matches <- c(unique(FT_plot$taxon[which(FT_plot$taxon %in% fac_spnames)]))
  #species in the facilitation but not in the trait data
  fac_only <- c(unique(fac_spnames[-which(fac_spnames %in% matches)]))
  #species in the trait but not in the facilitation data
  FT_only <- c(unique(FT_plot$taxon[-which(FT_plot$taxon %in% matches)]))
  
  if(i == 1) { 
    #only create FT sub for the first ID, afterwards we will just rbind to it
    FT_sub <- FT_plot[which(FT_plot$taxon %in% fac_spnames) , ] #subset the trait data to include only species encountered in the facilitation data
    
    #only create sp_matches for the first ID, and the first position (match) afterwards we will rbind to it
    sp_matches$ID <- fac_IDs[i]
    sp_matches$spname <- matches
    sp_matches$position <- "match"
    
  } else { #if i is not 1 we rbind to the table we created when i was 1
    temp_FT_sub <- FT_plot[which(FT_plot$taxon %in% fac_spnames) , ] #subset to inlcude only sp that are also in the facilitation data
    FT_sub <- rbind(FT_sub, temp_FT_sub) #rbind the temp_FT_sub files to the FT sub file
    
    #create table to rbind to the table we made when i was 1
    #species name matches
    temp_match <- cbind(rep(fac_IDs[i], length(matches)), matches , rep("match", length(matches)))
    colnames(temp_match) <- c("ID", "spname", "position")
    sp_matches <- rbind(sp_matches, temp_match) #rbind to sp_matches
  }
  #species only in facilitation data
  temp_match <- cbind(rep(fac_IDs[i], length(fac_only)) , fac_only, rep("fac_only", length(fac_only)))
  colnames(temp_match) <- c("ID", "spname", "position")
  sp_matches <- rbind(sp_matches, temp_match)
  
  #species only in FT data
  temp_match <- cbind(rep(fac_IDs[i], length(FT_only)), FT_only , rep("FT_only", length(FT_only)))
  colnames(temp_match) <- c("ID", "spname", "position")
  sp_matches <- rbind(sp_matches, temp_match)
  
}###end of loop
##FT sub did not try to fill any traits for species that are only in the fac data

write.csv(FT_sub, "Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species_graz_conserved.csv") #these species occur in that plot in the facilitation data

write.csv2(sp_matches, "Functional trait data\\Clean data\\filled_sp_matches.csv")


####Descriptive statistics####
#How many species were sampled in the facilitation survey but not in the trait survey?
#import spnames in the facilitation data
fac_species <- read.csv("Functional trait data\\Clean data\\facilitation_species_and_positions.csv", row.names = 1) |>
  rename(taxon = spname) |> 
  select(taxon) |> 
  distinct(taxon) 
fac_sp_total <- nrow(fac_species)

#import unfilled trait data
FT_species <- read.csv("Functional trait data\\Clean data\\FT_match_facilitation_plots.csv", row.names = 1) |> 
  select(taxon) |> 
  distinct(taxon)
FT_sp_total <- nrow(FT_species)

only_in_fac <- fac_species |> #species only sampled in the facilitation survey
  anti_join(FT_species, by = "taxon") |> 
  summarise(fac_only = n())
percent_only_in_fac <- only_in_fac/fac_sp_total *100

only_in_FT <- FT_species |> #species only sampled in the trait survey
  anti_join(fac_species, by = "taxon") |> 
  summarise(FT_only = n())
percent_only_in_FT <- only_in_FT/FT_sp_total *100


###Assess change in trait coverage###
#import unfilled trait data
FT_unfilled <- read.csv("Functional trait data\\Clean data\\FT_match_facilitation_plots_plotspecific_species.csv", row.names = 1) |> 
  pivot_longer(cols = c(MeanLL, MeanSLA, MeanLDMC, MeanLA, MaxH, MaxLS, percentN, percentC), 
               names_to = "trait", values_to = "value")

#import filled trait data
FT_filled <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species_graz_conserved.csv", row.names = 1)

###trait coverage in FT_unfilled
denom <- nrow(FT_unfilled) #this is the number of trait X species x plot combinations
#records that have trait values
unfilled_n_traits <- FT_unfilled |> 
                      filter(!is.na(value)) |>
                      summarise(unfilled_n_traits = n())
#percent of records that have trait values
unfilled_trait_coverage = (unfilled_n_traits$unfilled_n_traits/denom)*100

###trait coverage in FT_filled
filled_n_traits <- FT_filled |> 
  #trait_fill adds trait values from multiple places
  #ie one sp in a plot may have 3 different SLA values after filling
  #so we need to get the no of unique ID, species, trait combos first
  distinct(ID, taxon, trait, .keep_all = T) |> 
  summarise(filled_n_traits = n())

filled_trait_coverage = (filled_n_traits$filled_n_traits/denom)*100
#all records in FT_filled have values
#we use the same denominator because that is what we input into the trait filling function

###how much has trait coverage improved?
percent_improvement <- ((filled_n_traits$filled_n_traits - unfilled_n_traits$unfilled_n_traits)/denom)*100


