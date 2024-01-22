###Subset the trait data for the sites and sp that are also in the facilitation data
library(tidyverse)
library(tidylog)

wd <- "C:\\Users\\imke6\\Documents\\Msc Projek"
setwd(wd)

#Import the results from the facilitation analysis
fac_result <- read.csv("Facilitation data\\results\\NIntc_results_allcountries_26Sep.csv")

#Import the country_v2 data
data_files <- list.files("Facilitation data\\Countriesv2")
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")
for(i in 1:length(data_files)) {                              
  assign(paste0(countrynames[i]),                                   
         read.csv2(paste0("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation data\\Countriesv2\\",
                          data_files[i])))
}

#Import the information about the facilitation sites
siteinfo <- read.csv("Facilitation data\\BIODESERT_sites_information.csv")

#Import the DRYPOP trait data
drypop <- read.csv("Functional trait data\\drypop_20May.csv")



###Subset drypop to include only the plots that are in the facilitation results####
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
trait_subset <- drypop[which(drypop$Site %in% fac_sites) ,]

#Drop the columns that we don't need
colnames(trait_subset)
trait_subset <- trait_subset[ , which((colnames(trait_subset)) %in% c("Country", "Site", "Plot", "Genus", "Species", 
                                                       "MeanLL", "MeanSLA", "MeanLDMC", "MeanLA", "MaxH", "MaxLS"))]


##Now add ID, site_ID and graz to trait_subset
#Add site_ID first
sitename <- c(unique(siteinfo$SITE))
site_ID <- c(unique(siteinfo$SITE_ID))
add_site_ID <- data.frame(sitename, site_ID) #create dataframe with the variables we want to merge

trait_subset <- merge(x = trait_subset, y = add_site_ID, by.x = "Site", by.y = "sitename", all.x = TRUE, all.y = FALSE) #want to show all rows in x even if there is o match, but do not want to add rows from y that have no match

##Add ID 
#Make a plot identification variable by pasting the site ID and plot number. We will use this variable to merge on 
trait_subset$plotref <- paste(trait_subset$site_ID, trait_subset$Plot , sep = "_")
fac_result$plotref <- paste(fac_result$site_ID, fac_result$plot, sep = "_")

plotref <- unique(fac_result$plotref)
ID <- unique(fac_result$ID)
add_ID <- data.frame(plotref, ID) #create dataframe with the variables we want to merge

trait_subset <- merge(x = trait_subset, y = add_ID, by = "plotref", all.x = TRUE, all.y = FALSE) #want to show all rows in x even if there is o match, but do not want to add rows from y that have no match

#Add aridity and graz
add_var <- siteinfo[, c(1, 7, 62)]
add_var <- add_var[- which(is.na(add_var$ID)) , ]

trait_subset <- merge(x = trait_subset, y = add_var, by = "ID", all.x = TRUE, all.y = FALSE)
#export:
write.csv(trait_subset, "Functional trait data\\FT_match_facilitation_plots.csv")

####SUBSET TRAIT DATA TO INCLUDE ONLY FACILITATION SPECIES####
trait_subset <- read.csv("Functional trait data\\FT_match_facilitation_plots.csv", row.names = 1)

###GET A LIST OF SP IN THE FACILITATION DATA####
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
#write.csv(fac_species, "Functional trait data\\facilitation_species_and_positions.csv")
fac_species <- read.csv("Functional trait data\\facilitation_species_and_positions.csv", row.names = 1)

###Subset trait_subset to include only the species in the facilitation data####
##The species in the resulting dataset may not be present in that SPECIFIC PLOT in the facilitation data, they are just present someweher in the facilitatin data
#First we need to concatenate the genus and species names in trait_subset
trait_subset$sp_fullname <- paste(trait_subset$Genus, trait_subset$Species, sep = " ")

#Some names need to be changed in the trait data to match the facilitation data
trait_subset[trait_subset$sp_fullname %like% "Prosopis flex" , which(colnames(trait_subset) == "sp_fullname")] <- "Prosopis flexuosa"
trait_subset[trait_subset$sp_fullname %like% "Chuquiraga erinacea" , which(colnames(trait_subset) == "sp_fullname")] <- "Chuquiraga erinacea"

#subset to include only sp in the facilitation data
FT_match_facilitation_plots_general_species <- trait_subset[which(trait_subset$sp_fullname %in% c(fac_species$spname)) , ]
#write.csv(FT_match_facilitation_plots_general_species, "Functional trait data\\FT_match_facilitation_plots_general_species.csv")
trait_subset <- read.csv("Functional trait data\\FT_match_facilitation_plots_general_species.csv", row.names = 1)

###SUBSET THE FT DATA FOR THE SPECIES IN THAT PLOT IN THE FACILITATION DATA
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
  
  FT_plot <- trait_subset[which(trait_subset$ID == fac_IDs[i]) , ] #subset by plot in the trait data
  fac_sp_plot <- fac_species[which(fac_species$ID == fac_IDs[i]) , ] #subset by plot in fac_species
  
  fac_spnames <- c(unique(fac_sp_plot$spname)) #names of the species encountered in the facilitation dataset
  
  #species in both datasets
  matches <- c(FT_plot$sp_fullname[which(FT_plot$sp_fullname %in% fac_spnames)])
  #species in the facilitation but not in the trait data
  fac_only <- c(fac_spnames[-which(fac_spnames %in% matches)])
  #species in the trait but not in the facilitation data
  FT_only <- c(FT_plot$sp_fullname[-which(FT_plot$sp_fullname %in% matches)])
  
  if(i == 1) { 
    #only create FT sub for the first ID, afterwards we will just rbind to it
    FT_sub <- FT_plot[which(FT_plot$sp_fullname %in% fac_spnames) , ] #subset the trait data to include only species encountered in the facilitation data
    
    #only create sp_matches for the first ID, and the first position (match) afterwards we will rbind to it
    sp_matches$ID <- fac_IDs[i]
    sp_matches$spname <- matches
    sp_matches$position <- "match"
    
  } else { #if I is not 1 we rbind to the table we created when i was 1
    temp_FT_sub <- FT_plot[which(FT_plot$sp_fullname %in% fac_spnames) , ] #subset to inlcude only sp that are also in the facilitation data
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

#make a of the fac_only species without traits
fac_only_df <- sp_matches |> 
  filter(position == "fac_only") |> 
  #make all the columns the same as in FT_sub:
  rename(sp_fullname = spname) |> 
  add_column(MeanLL = NA, MeanSLA = NA, MeanLDMC = NA, 
             MeanLA = NA, MaxH = NA, MaxLS = NA) |> 
  select(!position) |>   #remove the position variable  
  mutate_at("ID", as.numeric) |>  #make ID a numeric variable  
  left_join(plotinfo, by = "ID") #merge the info associated with each plot


#now we can rbind fac_only_df to FT sub
FT_complete <- FT_sub |> 
  select(!c(Genus, Species, plotref)) |>  #remove these columns, they are uneccesary 
  bind_rows(fac_only_df)


write.csv(FT_complete, "Functional trait data\\FT_match_facilitation_plots_plotspecific_species_20jan.csv") #these species occur in that plot in the facilitation data
#Also, species that occur in the fac data only are in here with empty trait values.

write.csv2(sp_matches, "Functional trait data\\sp_matches_23Oct.csv")







