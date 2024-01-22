###Subset the trait data for the sites and sp that are also in the facilitation data
#And create a file with informatio about all the sites for which we have facilitation data
#install.packages("stringr")
library(stringr)
#install.packages("ggplot2")
library(ggplot2)

#Which are the species in the facilitation dataset (nurses and targets)
wd <- "C:\\Users\\user\\OneDrive\\Documents\\Msc Projek"
setwd(wd)

algeria <- read.csv("Facilitation data\\Countries\\Algeria_facilitation.csv")
argentina <- read.csv("Facilitation data\\Countries\\Argentina_facilitation.csv")
australia <- read.csv("Facilitation data\\Countries\\Australia_facilitation.csv")
chile <- read.csv("Facilitation data\\Countries\\CHile_facilitation.csv")
chinachong <- read.csv("Facilitation data\\Countries\\China_Chong_facilitation.csv")
chinaxin <- read.csv("Facilitation data\\Countries\\China_Xin_facilitation.csv")
iranabedi <- read.csv("Facilitation data\\Countries\\Iran_Abedi_facilitation.csv")
iranfarzam <- read.csv("Facilitation data\\Countries\\Iran_Farzam_facilitation.csv")
israel <- read.csv("Facilitation data\\Countries\\Israel_facilitation.csv")
namibiablaum <- read.csv("Facilitation data\\Countries\\Namibia_Blaum_facilitation.csv")
namibiawang <- read.csv("Facilitation data\\Countries\\Namibia_Wang_facilitation.csv")
southafrica <- read.csv("Facilitation data\\Countries\\SouthAfrica_facilitation.csv")
spainmaestre <- read.csv("Facilitation data\\Countries\\Spain_Maestre_facilitation.csv")
spainrey <- read.csv("Facilitation data\\Countries\\Spain_Rey_facilitation.csv")

#remove all the empty rows from each country dataset
algeria <- algeria[which(!is.na(algeria[,1])),]
southafrica <- southafrica[which(!is.na(southafrica[,1])),]
argentina <- argentina[which(!is.na(argentina[,1])),] #has no NA values
australia <- australia[which(!is.na(australia[,1])),]
chile <- chile[which(!is.na(chile[,1])),]
chinachong <- chinachong[which(!is.na(chinachong[,1])),]
chinaxin <- chinaxin[which(!is.na(chinaxin[,1])),]
iranabedi <- iranabedi[which(!is.na(iranabedi[,1])),]
iranfarzam<- iranfarzam[which(!is.na(iranfarzam[,1])),]
israel <- israel[which(!is.na(israel[,1])),]
namibiablaum <- namibiablaum[which(!is.na(namibiablaum[,1])), ]
namibiawang <- namibiawang[which(!is.na(namibiawang[,1])),]
spainmaestre <- spainmaestre[which(!is.na(spainmaestre[,1])),]
spainrey <- spainrey[which(!is.na(spainrey[,1])),]

#the species names in southafrica and chinachong are formatted incorrectly
#remove the underscore separating genus and species names in southafrica
southafrica[c("target genus", "target species")] <- str_split_fixed(southafrica$Species.within.quadrat, "_", 2)
southafrica$Species.within.quadrat <- paste(southafrica$`target genus`, southafrica$`target species`, sep = " " )
southafrica[which(southafrica$Species.within.quadrat == " ") , 16] <- NA

#remove the author names in chinachong
chinachong[c("target genus", "target species", "author", "author", "author")] <- str_split_fixed(chinachong$Species.within.quadrat, " ", 5)
chinachong$Species.within.quadrat <- paste(chinachong$`target genus`, chinachong$`target species`, sep = " " )

#Get a list of countries and sites in the facilitation dataset
countrylist <- list(algeria, southafrica, argentina, australia, chile, chinachong, chinaxin, iranabedi, iranfarzam, israel, namibiablaum,
                namibiawang, spainmaestre, spainrey)
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
fac_sites <- fac_sites[-which(is.na(fac_sites))] ##10 countries, 27 sites

####Import the drypop trait data and subset it to only include the sites also in the facilitation data
drypop <- read.csv("Functional trait data\\drypop_20May.csv")
trait_subset <- drypop[which(drypop$Site %in% fac_sites) ,]
#subset again to include only the species (only in the species within quadrat column for now) that is also in the facilitation data
common_sp <- c(southafrica$Species.within.quadrat, argentina$Species.within.quadrat, australia$Species, chile$Species.within.quadrat, 
               chinachong$Species.within.quadrat, chinaxin$Species.within.quadrat, iranabedi$Species.within.quadrat, iranfarzam$Species.within.quadrat, 
               israel$Species.within.quadrat, namibiablaum$Species.within.quadrat, namibiawang$Species.within.quadrat, 
               spainmaestre$Species.within.quadrat, spainrey$Species.within.quadrat)
trait_subset$speciesname <- paste(trait_subset$Genus, trait_subset$Species, sep = " ")
trait_subset <- trait_subset[which(trait_subset$speciesname %in% common_sp) ,]


#write to a csv file to explore in excel
write.csv(trait_subset, "Functional trait data\\trait_subset.csv")

###Create a table with the coordinates of each facilitation site to plot in qgis
sites_info <- read.csv("Facilitation data\\BIODESERT_sites_information.csv")
sites_info <- sites_info[-which(is.na(sites_info[,1])),]
#subset to include only the sites in the facilitation dataset
fac_sites <- sites_info[which(sites_info$SITE %in% fac_sites) ,]
#remove columns we dont need
fac_sites <- fac_sites[ , c(2,4,5,6,7,20,21,59, 60,61,62)]
write.csv(fac_sites, "Spatial data\\facilitation_sites_for_mapping.csv")

###make a graph showing the aridity of the first plot of each site
plots1 <- fac_sites[which(fac_sites$PLOT == 1) ,]
plots1$SITE <- as.factor(plots1$SITE)
plot(AI ~ SITE, data = plots1)
#test <- plots1[order(plots1$AI, decreasing = TRUE ), ]
#plot(AI ~ SITE, data = test)

plots1$SITE = with(plots1, reorder(SITE, AI, median))
mid <- median(plots1$AI)

library(RColorBrewer)
#display.brewer.all()
#display.brewer.pal(6, "BrBG")
pal <- brewer.pal(6, "BrBG")
brown <- pal[1]
sand <- pal[2]
turq <- pal[5]
teal <- pal[6]

vis <- 
  ggplot(plots1, aes(x=SITE, y=AI, fill = AI)) + 
  geom_bar(stat = "identity") +
  scale_fill_stepsn(n.breaks = 4, nice.breaks = F, breaks = c(0.16, 0.33, 0.49), colors = c(brown,  sand, turq , teal)) +
  coord_flip() + 
  theme_classic() +
  xlab("Site") +
  ylab("AI") +
  theme(legend.position = c(0.8, 0.5))
