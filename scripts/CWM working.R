library(tidyverse)
library(tidylog)
library(traitstrap)
library(RColorBrewer)

#trait data 
FT <- read.csv("Functional trait data\\Clean data\\FT_all_sites.csv", row.names = 1) |> 
  #some values of CoverBiodesert 100 exceed 100, therefore I assume it is the sum of the cover of a species over all 100 quadrats
  #therefore divide by 100
  mutate(coverBiodesert100 = coverBiodesert100/100, 
         #make sure that there are no 0 cover values.
         #values must be >0 or NA
         coverBiodesert100 = case_when(coverBiodesert100 == 0 ~ NA, .default = as.numeric(coverBiodesert100))) |> 
  #calculate C:N ratio
  mutate(C_N_ratio = percentC/percentN) |> 
  select(!c(percentC, percentN))



#Change FT to long format
FT_long <- FT |> 
  rename(coverBiodesert100_FT = coverBiodesert100) |> #keep this column here so that we know what the cover was in the trait data
  select(!c(Genus, Species, coverDryfun20)) |> 
  pivot_longer(cols = c(MeanLL, MeanSLA, MeanLDMC, MeanLA, MaxH, MaxLS, C_N_ratio), 
               names_to = "trait", values_to = "value")

FT_long$ID <- as.factor(FT_long$ID)
FT_long$site_ID <- as.factor(FT_long$SITE_ID)
FT_long$GRAZ <- factor(FT_long$GRAZ, levels = c(100, 0,1,2,3))

#Create the comm table containing the cover values for each species
#Use CoverBiodesert100
CovBio100 <- FT |> 
  select(c(ID, COU, SITE, SITE_ID, PLOT, taxon, coverBiodesert100, GRAZ, ARIDITY.v3)) |> 
  distinct() |> 
  filter(!is.na(coverBiodesert100)) #remove rows with NA cover values
CovBio100$ID <- as.factor(CovBio100$ID)
CovBio100$site_ID <- as.factor(CovBio100$SITE_ID)
CovBio100$GRAZ <- factor(CovBio100$GRAZ, levels = c(100, 0,1,2,3)) 


###Classify species as nurses and targets
#all species in each graz level
graz_allsp <- FT_long |> 
  group_by(GRAZ) |> 
  distinct(taxon) |> 
  ungroup()

#Get the nurses at each graz level
#import the nint results
nint_result <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\results\\NIntc_results_allcountries_6Feb2024.csv", row.names =1)   
graz_nurses <- nint_result |> 
  group_by(graz) |> 
  distinct(nurse) |> 
  ungroup()

##Loop to classify species in trait data##
FT_long$status <- NA
grazlevels <- c(0,1,2,3)

for(g in 1:length(grazlevels)) {
  nurses <- graz_nurses |> 
    filter(graz == grazlevels[g])
  
  allsp <- graz_allsp |> 
    filter(GRAZ == grazlevels[g])
  
  targets <- allsp[-which(allsp$taxon %in% nurses$nurse) , ]
  
    if(g == 1) {
    
    FT_graz <- FT_long[which(FT_long$GRAZ == grazlevels[g]) , ]
    
    for(i in 1:nrow(FT_graz)) {
      if(FT_graz$taxon[i] %in% nurses$nurse) {
        FT_graz[i, which(colnames(FT_graz) == "status")] <- "nurse" }
      
      if(FT_graz$taxon[i] %in% targets$taxon) {
        FT_graz[i, which(colnames(FT_graz) == "status")] <- "target" }} 
    
      }else {
      temp_FT_graz <- FT_long[which(FT_long$GRAZ == grazlevels[g]) , ]
      
      for(i in 1:nrow(temp_FT_graz)) {
        if(temp_FT_graz$taxon[i] %in% nurses$nurse) {
          temp_FT_graz[i, which(colnames(temp_FT_graz) == "status")] <- "nurse" }
        
        if(temp_FT_graz$taxon[i] %in% targets$taxon) {
          temp_FT_graz[i, which(colnames(temp_FT_graz) == "status")] <- "target" }}
    
      FT_graz <- rbind(FT_graz, temp_FT_graz)
    }} ##loop ends here##



##Loop to classify species in community data##
CovBio100$status <- NA
grazlevels <- c(0,1,2,3)

for(g in 1:length(grazlevels)) {
  nurses <- graz_nurses |> 
    filter(graz == grazlevels[g])
  
  allsp <- graz_allsp |> 
    filter(GRAZ == grazlevels[g])
  
  targets <- allsp[-which(allsp$taxon %in% nurses$nurse) , ]
  
  if(g == 1) {
    
    comm_graz <- CovBio100[which(CovBio100$GRAZ == grazlevels[g]) , ]
    
    for(i in 1:nrow(comm_graz)) {
      if(comm_graz$taxon[i] %in% nurses$nurse) {
        comm_graz[i, which(colnames(comm_graz) == "status")] <- "nurse" }
      
      if(comm_graz$taxon[i] %in% targets$taxon) {
        comm_graz[i, which(colnames(comm_graz) == "status")] <- "target" }} 
    
  }else {
    temp_comm_graz <- CovBio100[which(CovBio100$GRAZ == grazlevels[g]) , ]
    
    for(i in 1:nrow(temp_comm_graz)) {
      if(temp_comm_graz$taxon[i] %in% nurses$nurse) {
        temp_comm_graz[i, which(colnames(temp_comm_graz) == "status")] <- "nurse" }
      
      if(temp_comm_graz$taxon[i] %in% targets$taxon) {
        temp_comm_graz[i, which(colnames(temp_comm_graz) == "status")] <- "target" }}
    
    comm_graz <- rbind(comm_graz, temp_comm_graz)
  }} ##loop ends here##



FT_filled <- trait_fill(
  comm = comm_graz, #community data with abundance values
  traits = FT_graz, #trait data
  abundance_col = "coverBiodesert100", #the column with the abundance data
  taxon_col = "taxon", #column with the speciesnames (must be the same in both datasets)
  value_col = "value",
  trait_col = "trait",
  
  global = FALSE, #do not calculate traits at the global scale
  
  keep_all = FALSE, #do not keep trait data at all availible levels, only on the finest scale availible
  
  # specifies sampling hierarchy
  scale_hierarchy = c("COU", "SITE_ID", "ID"),
  
  treatment_col = "GRAZ",
  treatment_level = "COU", #hierarchy level at which you will find the same graz treatment level
  
  other_col = "status",
  
  # min number of samples
  min_n_in_sample = 1 #if there is one trait value in the sample, do not search for values higher up in the hierarchy
)




#Get bootstrapped distributions
FT_bs <- trait_np_bootstrap(FT_filled, sample_size = 100, nrep = 50)

#Look at BS moments
FT_moments <- trait_summarise_boot_moments(
  FT_bs, 
  parametric = TRUE) #uses moment +- 1sd

#plot them
pal <- brewer.pal(3, "Dark2")[1:2]

cwm <- ggplot(FT_moments, aes(x = GRAZ_comm, y = mean, fill = status)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  #geom_errorbar(aes(ymin = ci_low_mean, ymax = ci_high_mean, group = trait)) +
  facet_wrap(~trait, scale = "free_y") +
  scale_fill_manual(labels = c("dominant", "target"), 
                    values = pal) +
  labs(x = " ", y = "CWM", fill = " ") +
  scale_x_discrete(labels = c("Ungrazed", "Low", "Medium", "High")) +
  theme_classic() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 90))

ggsave("cwm_bygraz.png", cwm, height = 1300, width = 1600, units = "px",
       path = "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Figures")
