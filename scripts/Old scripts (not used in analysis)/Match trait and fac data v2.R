###Subset the trait data for the sites and sp that are also in the facilitation data
library(tidyverse)
library(tidylog)
library(DescTools)


#Import the country_v3 data
data_files <- list.files("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Countriesv3")
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")
for(i in 1:length(data_files)) {                              
  assign(paste0(countrynames[i]),                                   
         read.csv2(paste0("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Countriesv3\\",
                          data_files[i])))
}

#Import the information about the facilitation sites
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\BIODESERT_sites_information.csv")

#Import the DRYPOP trait data
drypop <- read.csv("Functional trait data\\Raw data\\drypop_20May.csv")



###Subset drypop to include only the plots that are in the facilitation data####
#DRYPOP only has the site names, not ID or site_ID, so we have to work with the names
#Get a list of countries and sites in the facilitation dataset
countrylist <- list(algeria, argentina, australia, chile, chinachong, chinaxin, iranabedi, iranfarzam, israel, 
                    namibiablaum, namibiawang, southafrica, spainmaestre, spainrey)
#Get all the sites in each country
sitedf <- data.frame(matrix(nrow = length(100), ncol = 5))
colnames(sitedf) <- c("country", "site1", "site2", "site3", "site4")
l = 1
for (i in countrylist) {
  country <- unique(i$COU)
  site <- unique(i$SITE)
  
  sitedf[l,1] <- country[1]
  sitedf[l,2] <- site[1]
  sitedf[l,3] <- site[2]
  sitedf[l,4] <- site[3]
  sitedf[l,5] <- site[4]
  l = l+1
}

#Get a vector of all the sites in all the countries
fac_sites <- c(sitedf[,2], sitedf[,3], sitedf[,4])
fac_sites <- fac_sites[-which(is.na(fac_sites))] 

#subset drypop to include only the sites in sitedf
trait_subset_sites <- drypop |> 
  filter(Site %in% fac_sites) |> 
  select(Country, Site, Plot, Genus, Species, coverDryfun20,coverBiodesert100, Life_form, #keep only the columns we need
         MeanLL, MeanSLA, MeanLDMC, MeanLA, MaxH, MaxLS)

unique(trait_subset_sites$Site)



##Now add ID, site_ID and graz to trait_subset
#Add site_ID first
Site <- c(unique(siteinfo$SITE))
site_ID <- c(unique(siteinfo$SITE_ID))
add_site_ID <- data.frame(Site, site_ID) #create dataframe with the variables we want to merge

trait_subset_sites <- trait_subset_sites |> 
  left_join(add_site_ID, by = "Site")

##Now add ID, site_ID and graz to trait_subset_sites 
#Make a plot identification variable by pasting the site ID and plot number. We will use this variable to merge on 
#plotref of trait data
trait_subset_sites$plotref <- paste(trait_subset_sites$site_ID, trait_subset_sites$Plot , sep = "_")

#plotref of fac data
add_ID <- data.frame(plotref = c(rep(NA,3)), ID = c(rep(NA,3))) #create dataframe with plotref and ID of the fac data

#loop to get the plotrefs and IDs of each country
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")

for (t in 1:length(countrynames)) {
  cou <- get(countrynames[t])
  
  if(t ==1) {
  
    add_ID$plotref <- c(unique(paste(cou$SITE_ID, cou$PLOT, sep = "_")))
    add_ID$ID <- c(unique(cou$ID))
  
  }else {
    
    temp_cou_plotrefs <- c(unique(paste(cou$SITE_ID, cou$PLOT, sep = "_")))
    temp_cou_IDs <- c(unique(cou$ID))
    
    temp_rows <- cbind(plotref = temp_cou_plotrefs, ID = temp_cou_IDs)
    
    add_ID <- rbind(add_ID, temp_rows)
    
  }
}

trait_subset_plots <- trait_subset_sites |> 
  right_join(add_ID, by = "plotref") #do a right join, because we only want to keep ID's that are in the fac data
  


#Add aridity and graz
add_var <- siteinfo[, c(1, 7, 62)]
add_var <- add_var[- which(is.na(add_var$ID)) , ]
add_var$ID <- as.character(add_var$ID)

trait_subset_plots <- trait_subset_plots |> 
  left_join(add_var, by = "ID")


####NOW WE FIX THE TRAIT NAMES TO MATCH THE FACILITATION DATA###
##CLean up easy mistakes in names:
trait_subset_plots <- trait_subset_plots |> 
  mutate(taxon = str_c(Genus, Species, sep = " "), #spnames in one column
         taxon = str_squish(taxon), #remove spaces before and after string
         taxon = str_to_sentence(taxon)) |>  
  #remove infraspecies ranks
  separate_wider_delim(taxon, delim = " ",   
                       names = c("split1", "split2", "split3",  "split4"), too_few = "align_start") |> 
  mutate(taxon = str_c(split1, split2, sep = " ")) |> 
  select(!c(split1,split2,split3, split4))
  
##NOW RUN NAME TRAIL
#The sheet names_and_synonyms contains names and their synonyms accross the facilitation, trait and quadrat data.
#This list was made by comparing only_in_fac to the quadrat data in the script called trait name changing.
name_trail <- read_excel("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Name changes\\names_and_synonyms.xlsx") |> 
  mutate(correct_name = str_squish(correct_name), ##remove spaces before or after strings, and replace internal whitespace with a single space
         synonym1 = str_squish(synonym1),
         synonym2 = str_squish(synonym2))

#create synonyms variable which is all the synonyms concatenated
name_trail$synonyms <- paste(name_trail$synonym1, name_trail$synonym2, sep = "; ") 

name_trail <- name_trail %>% 
  select(correct_name, synonyms) # only keep "correct_name and "synonyms"


#create a dataframe that stores information about samples for which species names are modified by the code below
change_tracker <- data.frame( 
  old_spec = character(), 
  new_spec = character(), 
  stringsAsFactors=FALSE) 

###Now standardise each name to the name trail:

data_harmony <- trait_subset_plots
  
for (i in 1:nrow(data_harmony)) {
  old_sp <- data_harmony[i, which(colnames(data_harmony) == "taxon")]
  new_sp <- NA
  

    found <- FALSE
    for (j in 1:nrow(name_trail)) { # looks whether species name is a synonym and replaces it with the true_name if it is found to be a synonym
      found <- grepl(old_sp, name_trail[j, 2]) 
      
      if (found){ # only runs if the species is a synonym
        new_sp <- name_trail[j, 1] # finds the true name of the species and saves it 
        break
      }
    }
    
    if (found) { # replaces the species in the trait database with the saved true name if "found" is "TRUE"
      data_harmony[i, which(colnames(data_harmony) == "taxon")] <- new_sp 
      
      # add a new row with information about change to the change trackers dataset
      change_tracker[i, 1] <- old_sp
      change_tracker[i, 2] <- new_sp
    }
  } #close the loop through the names in data_harmony
head(change_tracker, 10)

#export:
write.csv(data_harmony, "Functional trait data\\FT_match_facilitation_plots.csv")
#this data has only the plots in the fac data and names are standardised


####SUBSET TRAIT DATA TO INCLUDE ONLY FACILITATION SPECIES####
trait_subset_plots <- read.csv("Functional trait data\\FT_match_facilitation_plots.csv", row.names = 1)

##GET A LIST OF SP IN THE FACILITATION DATA###
##Classify each species as a nurse, or as growing in a bare or open microsite

#Bind all the country data together, use only the columns we need
vars <- c("ID", "Microsite", "ID_Microsite", "Species.within.quadrat" )
data <- rbind(algeria[ , which(colnames(algeria) %in% vars)], argentina[ , which(colnames(argentina) %in% vars)], australia[ , which(colnames(australia) %in% vars)], chile[ , which(colnames(chile) %in% vars)], 
              chinachong[ , which(colnames(chinachong) %in% vars)], chinaxin[ , which(colnames(chinaxin) %in% vars)], iranabedi[ , which(colnames(iranabedi) %in% vars)], iranfarzam[ , which(colnames(iranfarzam) %in% vars)], 
              israel[ , which(colnames(israel) %in% vars)], namibiablaum[ , which(colnames(namibiablaum) %in% vars)], namibiawang[ , which(colnames(namibiawang) %in% vars)], southafrica[ , which(colnames(southafrica) %in% vars)],  
              spainmaestre[ , which(colnames(spainmaestre) %in% vars)], spainrey[ , which(colnames(spainrey) %in% vars)])

IDlist <- c(unique(data$ID)) #list of the unique ID's corresponding to each plot

#make a file to put the speciesnames in
fac_species <- data.frame(matrix(nrow = 1, ncol = 3))
colnames(fac_species) <- c("ID", "condition", "spname") #condition will be nurse if it is a nurse species
#condition will be micro_nurse if the species was growing in a nurse microsite
#condition will be micor_bare if the species was growing in a bare microsite

for(i in 1:length(IDlist)) {
  plot <- data[which(data$ID == IDlist[i]) , ] #isolate a single plot
  
  nurse_names <- unique(plot[which(plot$Microsite == 2) , which(colnames(plot) == "ID_Microsite")]) #get the nurse names in that plot
  
  
  if (i == 1) { #only make a new table for the first plot and nurse condition, we will just attach the results of subsequent plots to this table
  
    fac_species$ID <- IDlist[i]
    fac_species$condition <- "nurse"
    fac_species$spname <- nurse_names
  
  } else { #for subsequent plots, put the results in temp_mat and bind it to the result_prelim we make above
    temp_mat <- cbind(rep(IDlist[i], length(nurse_names)) , rep("nurse", length(nurse_names)) , nurse_names)
    colnames(temp_mat) <- c("ID", "condition", "spname")
    fac_species <- rbind(fac_species, temp_mat)
    }
  
 
   #now result_prelim already exists for us to bind the names of the other conditions to
  micro_nurse_names <- unique(plot[which(plot$Microsite == 2) , which(colnames(plot) == "Species.within.quadrat")])
  micro_nurse_names <- micro_nurse_names[!is.na(micro_nurse_names)] #remove NA
  temp_mat <- cbind(rep(IDlist[i], length(micro_nurse_names)) , rep("micro_nurse", length(micro_nurse_names)) , micro_nurse_names)
  colnames(temp_mat) <- c("ID", "condition", "spname")
  fac_species <- rbind(fac_species, temp_mat)
  
  micro_bare_names <- unique(plot[which(plot$Microsite == 1) , which(colnames(plot) == "Species.within.quadrat")])
  micro_bare_names <- micro_bare_names[!is.na(micro_bare_names)] #remove NA
  temp_mat <- cbind(rep(IDlist[i], length(micro_bare_names)) , rep("micro_bare", length(micro_bare_names)) , micro_bare_names)
  colnames(temp_mat) <- c("ID", "condition", "spname")
  fac_species <- rbind(fac_species, temp_mat)
  
  #Now move to next ID
  
} ###end of loop
##Remember the same sp my have multiple conditions
##write to a csv file
write.csv(fac_species, "Functional trait data\\facilitation_species_and_positions.csv")
fac_species <- read.csv("Functional trait data\\facilitation_species_and_positions.csv", row.names = 1)

###Subset trait_subset_plots to include only the species in the facilitation data####
##The species in the resulting dataset may not be present in that SPECIFIC PLOT in the facilitation data, they are just present someweher in the facilitatin data
FT_match_facilitation_plots_general_species <- trait_subset_plots |> 
  filter(taxon %in% fac_species$spname)
write.csv(FT_match_facilitation_plots_general_species, "Functional trait data\\FT_match_facilitation_plots_general_species_31jan.csv")

###SUBSET THE FT DATA FOR THE SPECIES IN THAT PLOT IN THE FACILITATION DATA####
#The resulting table will be called FT_sub
###ALSO FIND OUT WHICH SPECIES ARE IN BOTH DATASETS OR ONLY ON ONE?
#which sp are only in the facilitation data, but not in the trait data, this table is sp_matches
sp_matches <- data.frame(matrix(nrow = 11, ncol = 3))
colnames(sp_matches) <- c("ID", "spname", "position")
#Position = fac_only, if the sp is only in the facilitation data, FT_only if it is only in the FT data, match if the sp is present in both datasets

#We want to do this separately for each plot, because it will be a plot level analysis
fac_IDs <- c(unique(as.numeric(fac_species$ID))) #plot ID's from the facilitation data

##Loop to subset trait_subset for the facilitation sp and to get the names of species in both or only one dataset
for(i in 1:length(fac_IDs)) {
  
  FT_plot <- trait_subset_plots[which(trait_subset_plots$ID == fac_IDs[i]) , ] #subset by plot in the trait data
  fac_sp_plot <- fac_species[which(fac_species$ID == fac_IDs[i]) , ] #subset by plot in fac_species
  
  fac_spnames <- c(unique(fac_sp_plot$spname)) #names of the species encountered in the facilitation dataset
  
  #species in both datasets
  matches <- c(FT_plot$taxon[which(FT_plot$taxon %in% fac_spnames)])
  #species in the facilitation but not in the trait data
  fac_only <- c(fac_spnames[-which(fac_spnames %in% matches)])
  #species in the trait but not in the facilitation data
  FT_only <- c(FT_plot$taxon[-which(FT_plot$taxon %in% matches)])
  
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

#Now we need to add the species that are only present in the facilitation data, and give them NA trait values
#get the info associated with each ID
plotinfo <- FT_sub |> 
  select(ID, Site, Country, Plot, site_ID, GRAZ, ARIDITY.v3) |> 
  distinct(ID, .keep_all = TRUE) #we only need one record for each ID

#make a table of the fac_only species without traits
fac_only_df <- sp_matches |> 
  filter(position == "fac_only") |> 
  #make all the columns the same as in FT_sub:
  rename(taxon = spname) |> 
  add_column(MeanLL = NA, MeanSLA = NA, MeanLDMC = NA, 
             MeanLA = NA, MaxH = NA, MaxLS = NA) |> 
  select(!position) |>   #remove the position variable  
  mutate_at("ID", as.numeric) |>  #make ID a numeric variable  
  left_join(plotinfo, by = "ID") #merge the info associated with each plot


#now we can rbind fac_only_df to FT sub
FT_complete <- FT_sub |> 
  select(!c(Genus, Species, plotref)) |>  #remove these columns, they are uneccesary 
  bind_rows(fac_only_df)


write.csv(FT_complete, "Functional trait data\\FT_match_facilitation_plots_plotspecific_species_31jan.csv") #these species occur in that plot in the facilitation data
#Also, species that occur in the fac data only are in here with empty trait values.

write.csv2(sp_matches, "Functional trait data\\sp_matches_31jan.csv")







