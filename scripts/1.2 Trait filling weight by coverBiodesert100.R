###TRAIT FILLING###
#This is the script for trait filling and subsetting the filled FT data for the facilitation plots##
#Weight the trait filling by CoverdryFun20
library(tidyverse)
library(tidylog)
library(traitstrap)

#we will import FT data for all sites that have FT data. 
#Then do trait filling
#Then subset for facilitation plots to have most complete trait data. 

#trait data 
FT <- read.csv("Functional trait data\\Clean data\\FT_all_sites.csv", row.names = 1) |> 
  #some values of CoverBiodesert 100 exceed 100, therefore I assume it is the sum of the cover of a species over all 100 quadrats
  #therefore divide by 100
  mutate(coverBiodesert100 = coverBiodesert100/100, 
         #make sure that there are no 0 cover values.
         #values must be >0 or NA
         coverBiodesert100 = case_when(coverBiodesert100 == 0 ~ NA, .default = as.numeric(coverBiodesert100)))



#Change FT to long format
FT_long <- FT |> 
  rename(coverBiodesert100_FT = coverBiodesert100) |> #keep this column here so that we know what the cover was in the trait data
  select(!c(Genus, Species, coverDryfun20)) |> 
  pivot_longer(cols = c(MeanLL, MeanSLA, MeanLDMC, MeanLA, MaxH, MaxLS, percentN, percentC), 
                       names_to = "trait", values_to = "value")

FT_long$ID <- as.factor(FT_long$ID)
FT_long$site_ID <- as.factor(FT_long$SITE_ID)
FT_long$GRAZ <- factor(FT_long$GRAZ, levels = c(0,1,2,3))

#Create the comm table containing the cover values for each species
#Use CoverBiodesert100
CovBio100 <- FT |> 
  select(c(ID, COU, SITE, SITE_ID, PLOT, taxon, coverBiodesert100, GRAZ, ARIDITY.v3)) |> 
  distinct() |> 
  filter(!is.na(coverBiodesert100)) #remove rows with NA cover values
CovBio100$ID <- as.factor(CovBio100$ID)
CovBio100$site_ID <- as.factor(CovBio100$SITE_ID)
CovBio100$GRAZ <- factor(CovBio100$GRAZ, levels = c(0,1,2,3)) 


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
  scale_hierarchy = c("COU", "SITE_ID", "ID"),
  
  other_col = "GRAZ",
  
  # min number of samples
  min_n_in_sample = 1 #if there is one trait value in the sample, do not search for values higher up in the hierarchy
)
FT_filled_sorted <- arrange(FT_filled, ID)
write.csv(FT_filled, "Functional trait data\\Clean data\\FT_filled_all_sites.csv")
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
FT_filled <- read.csv("Functional trait data\\Clean data\\FT_filled_all_sites.csv", row.names = 1)

#Import the country_v3 data
data_files <- list.files("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Countriesv3")
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")
for(i in 1:length(data_files)) {                              
  assign(paste0(countrynames[i]),                                   
         read.csv2(paste0("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Countriesv3\\",
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
write.csv(FT_filled_facplots, "Functional trait data\\Clean data\\FT_filled_match_facilitation_plots.csv")



###SUBSET THE FT_filled_facplots FOR THE SPECIES IN THAT PLOT IN THE FACILITATION DATA####
FT_filled_facplots <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots.csv", row.names = 1)

###ALSO FIND OUT WHICH SPECIES ARE IN BOTH DATASETS OR ONLY ON ONE?
#import spnames in the facilitation data
fac_species <- read.csv("Functional trait data\\Clean data\\facilitation_species_and_positions.csv", row.names = 1)

#which sp are only in the facilitation data, but not in the trait data, this table is sp_matches
sp_matches <- data.frame(matrix(nrow = 11, ncol = 3))
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

write.csv(FT_sub, "Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species.csv") #these species occur in that plot in the facilitation data

write.csv2(sp_matches, "Functional trait data\\Clean data\\filled_sp_matches.csv")


###Assess change in trait coverage###
#import unfilled trait data
FT_unfilled <- read.csv("Functional trait data\\Clean data\\FT_match_facilitation_plots.csv", row.names = 1) |> 
  pivot_longer(cols = c(MeanLL, MeanSLA, MeanLDMC, MeanLA, MaxH, MaxLS, percentN, percentC), 
               names_to = "trait", values_to = "value")

#import filled trait data
FT_filled <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots.csv", row.names = 1)

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


