##Functional distance analysis
library(tidyverse)
library(tidylog)
library(ggplot2)
library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)
library(MuMIn)
#library(multcomp)
#library(multcompView)

#This script gets th euclidean distance between the nurse species and every  other sp in the replicate. 
#It doesn't average those distances, if there are 3 target sp in the replicate there will be 3 distances.
#It uses standardised trait values
#The resulting data is called twosp_dist

####Create function to standardise sp X trait matrices####
standard_trait_matrix <- function(trait_matrix, traitlist) {
  std_trait_matrix <- trait_matrix
  for(t in 1:length(traitlist)) {
    
    #get the grand mean and sd for each trait
    m <- trait_matrix |> 
      summarise(m = mean(trait_matrix[ , which(colnames(trait_matrix) == traitlist[t])]))
    
    sd <- trait_matrix |> 
      summarise(m = sd(trait_matrix[ , which(colnames(trait_matrix) == traitlist[t])]))
    
    for (i in 1:nrow(trait_matrix)) {
      
      raw_value <- trait_matrix[i , which(colnames(trait_matrix) == traitlist[t])]
      #replace the raw value with the standardised value using m and sd
      std_trait_matrix[i , which(colnames(std_trait_matrix) == traitlist[t])] <- (raw_value - m)/sd
    }
  }
    return(std_trait_matrix)
}


###Create function to retrieve pairwise distances between nurse and target plants####
#This function gets the distance between the dominant species and every other species in the replicate. And then classifies that distance as 
#nurse_bare_only = dist between nurse and a species that occurs only in bare microsite
#nurse_nurse_only = dist between nurse and a species that occurs only in nurse microsite
#nurse_both = dist between nurse and a species that occurs in both microsites

pairwise_fdist <- function(distmat, #euclidean distances between species pairs based on functional traits
                           sp_positions) #whther species are nurses, or grow only in a specific microsite
{ #function returns a dataframe with euclidean distance between nurse and target species and the microsite affinity of the target.
  
  IDlist <- unique(sp_positions$ID)
  
  twosp_dist <- cbind("ID" = character(), "replicate" = character(), 
                      "euclidean_dist" = numeric(), "grouping" = character(), "nurse" = character(), "target" = character())
  
  for(i in 1:length(IDlist)) {
    plot <- sp_positions |> 
      filter(ID == IDlist[i])
    
    replist <- unique(plot$replicate) 
    
    for(r in 1:length(replist)) {
      #Get the names of nurse, bare only, nurse only and both species
      NURSE <- sp_positions |> 
        filter(ID == IDlist[i], replicate == replist[r], position == "nurse_species") |> 
        select(taxon)
      NURSE <- NURSE$taxon
      
      bare_only <- sp_positions |> 
        filter(ID == IDlist[i], replicate == replist[r], position == "bare_only") |> 
        select(taxon)
      bare_only <- c(bare_only$taxon)
      
      nurse_only <- sp_positions |> 
        filter(ID == IDlist[i], replicate == replist[r], position == "nurse_only") |> 
        select(taxon)
      nurse_only <- c(nurse_only$taxon)
      
      both <- sp_positions |> 
        filter(ID == IDlist[i], replicate == replist[r], position == "both") |> 
        select(taxon)
      both <- c(both$taxon)
      
      ###
      for(b in 1:length(bare_only)) {
        D_bare_only <- 
          mean(distmat[which(rownames(distmat) == NURSE), which(colnames(distmat) == bare_only[b])])
        
        if(b == 1){
          twosp_bare_only <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                   "euclidean_dist" = D_bare_only, "grouping" = "nurse_bare_only", 
                                   "nurse" = NURSE, "target" = bare_only[b])
        } else {
          temp_twosp_bare_only <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                        "euclidean_dist" = D_bare_only, "grouping" = "nurse_bare_only", 
                                        "nurse" = NURSE, "target" = bare_only[b])
          twosp_bare_only <- rbind(twosp_bare_only, temp_twosp_bare_only)
        }}
      
      ###
      for(n in 1:length(nurse_only)) {
        D_nurse_only <- 
          mean(distmat[which(rownames(distmat) == NURSE), which(colnames(distmat) == nurse_only[n])])
        
        if(n == 1){
          twosp_nurse_only <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                    "euclidean_dist" = D_nurse_only, "grouping" = "nurse_nurse_only", 
                                    "nurse" = NURSE, "target" = nurse_only[n])
        } else {
          temp_twosp_nurse_only <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                         "euclidean_dist" = D_nurse_only, "grouping" = "nurse_nurse_only", 
                                         "nurse" = NURSE, "target" = nurse_only[n])
          twosp_nurse_only <- rbind(twosp_nurse_only, temp_twosp_nurse_only)
        }}
      
      ###
      for(z in 1:length(both)) {
        D_both <- 
          mean(distmat[which(rownames(distmat) == NURSE), which(colnames(distmat) == both[z])])
        
        if(z == 1){
          twosp_both <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                              "euclidean_dist" = D_both, "grouping" = "nurse_both", 
                              "nurse" = NURSE, "target" = both[z])
        } else {
          temp_twosp_both <- twosp_both <- cbind("ID" = IDlist[i], "replicate" = replist[r], 
                                                 "euclidean_dist" = D_both, "grouping" = "nurse_both", 
                                                 "nurse" = NURSE, "target" = both[z])
          twosp_both <- rbind(twosp_both, temp_twosp_both)
        }}
      
      twosp_dist <- rbind(twosp_dist, twosp_bare_only,twosp_nurse_only, twosp_both)
    }
  }
  #The distance is NaN if nurse_only, bare_only or both is NA
  
  return(twosp_dist)
}


##Get the species that are nurses, growing with nurses, or in the bare microsite in each country####
#We require raw country data
data_files <- list.files("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\Countriesv3")
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")
for(i in 1:length(data_files)) {                              
  assign(paste0(countrynames[i]),                                   
         read.csv2(paste0("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Countriesv3\\",
                          data_files[i])))
}

countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")

##LOOP starts here
l = 1
for (t in 1:length(countrynames)) {
  cou <- get(countrynames[t])
  
  IDlist <- unique(cou$ID)
  
    for(i in 1:length(IDlist)) {
    plot <- cou |> 
      filter(ID == IDlist[i])
    
    replist <- unique(plot$Number.of.replicate)
  
      for(r in 1:length(replist)) {
        one_rep <- plot |> 
          filter(Number.of.replicate == replist[r])
        
        NURSE <- one_rep |> 
          filter(Microsite == 2) |> 
          distinct(ID_Microsite)
        
        bare_sp <- one_rep |> #sp in bare microsite
          filter(Microsite == 1) |> 
          select(Species.within.quadrat)
        
        nurse_sp <- one_rep |> #sp in nurse microsite
          filter(Microsite == 2) |> 
          select(Species.within.quadrat)
        
        #get the species that are in both nurse and bare microsites
        both_sp <- c(bare_sp[match(bare_sp$Species.within.quadrat, nurse_sp$Species.within.quadrat) , ])
        both_sp <- both_sp[!is.na(both_sp)]
        #if there are no species in both microsites, assign NA
        if(length(both_sp) == 0) {
          both_sp <- NA
        }
        
        #only if both_sp is not NA:
        if(FALSE %in% is.na(both_sp)) {
        
        #species ONLY in nurse microsite:
        nurse_only <- nurse_sp[-which(nurse_sp$Species.within.quadrat %in% both_sp), ]
        #if there are no species, assign NA
        if(length(nurse_only) == 0) {
          nurse_only <- NA
        }
        
        #species ONLY in bare microsites
        bare_only <- bare_sp[-which(bare_sp$Species.within.quadrat %in% both_sp) , ]
        #if there are no species, assign NA
        if(length(bare_only) == 0) {
          bare_only <- NA
        }
        
        #if both sp is NA:
        }else {
          nurse_only <- nurse_sp$Species.within.quadrat
          bare_only <- bare_sp$Species.within.quadrat
        }
        
        #only for the first run of the loop:
        if(l == 1){
        #put species and their classifications in the dataframe
        nurse_identity_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = NURSE$ID_Microsite, position = "nurse_species")
        nurse_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = nurse_only, position = "nurse_only")
        bare_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = bare_only, position = "bare_only")
        both_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = both_sp, position = "both")
        
        sp_positions <- rbind(nurse_identity_df, nurse_df, bare_df, both_df)
        
        } else { #for subsequent runs we rbind to the df we made whn l was 1
          nurse_identity_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = NURSE$ID_Microsite, position = "nurse_species")
          nurse_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = nurse_only, position = "nurse_only")
          bare_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = bare_only, position = "bare_only")
          both_df <- cbind(ID = IDlist[i], replicate = replist[r], taxon = both_sp, position = "both")
          
          sp_positions <- rbind(sp_positions, nurse_identity_df, nurse_df, bare_df, both_df)
        }
        
        l = l+1
      }#close loop through replicates
  }#close loop through ID's
}#close loop through countries

sp_positions <- as.data.frame(sp_positions)
#save the table 
write.csv(sp_positions, "Functional trait data\\results\\sp_positions.csv")


####Functional distance in 7 dimensional space####
sp_positions <- read.csv("Functional trait data\\results\\sp_positions.csv", row.names = 1) 
#From the filled trait data for plotspecific species
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species_graz_conserved.csv", row.names = 1)

complete_FT_mean <- FT |> #get the mean trait value for each sp accross the whole dataset. we will use this to fill missing data values
  group_by(taxon,trait) |> 
  summarise(mean_value = mean(value))

IDlist <- c(unique(sp_positions$ID))[10]

for(p in 1:length(IDlist)) {
  positions_plot <- sp_positions[which(sp_positions$ID == IDlist[p]) , ]
  FT_plot <- FT[which(FT$ID == IDlist[p]) , ]
  
  #now make FT wide and standardise it
  #Get the mean of each trait for each sp
  FT_mean <- FT_plot |> 
    group_by(taxon,trait) |> 
    summarise(mean_value = mean(value))
  
  #Find out if traits are missing:
  splist <- c(unique(FT_mean$taxon)) #species in plot
  traits_required <- c("MaxH","MaxLS","MeanLA","MeanLDMC","MeanLL","MeanSLA","percentC", "percentN")
  
  for (s in 1:length(splist)) {
  traits_present <- FT_mean[which(FT_mean$taxon == splist[s]) , ]$trait
  #get missing traits
  traits_tofill <- c(traits_required[which(is.na(match(traits_required, traits_present)))])
  
  #if there are missing traits:
  if(length(traits_tofill) > 0) {
  for(f in 1:length(traits_tofill)) {
  #get the filler value
  filler_value <- complete_FT_mean |> 
    filter(taxon == splist[s], trait == traits_tofill[f])
  #if there is a filler available:
  if(nrow(filler_value) >0) {
  filler_row = data.frame(taxon = splist[s], trait = traits_tofill[f], mean_value = filler_value$mean_value)
  FT_mean <- rbind(FT_mean, filler_row)
  } else {
    #if there is no filler, just insert an NA value
    filler_row = data.frame(taxon = splist[s], trait = traits_tofill[f], mean_value = NA)
    FT_mean <- rbind(FT_mean, filler_row)
  }
  }}}
  
  #Get it into wide format
  FT_wide <- FT_mean |>
    pivot_wider(names_from = trait, values_from = mean_value) |> 
    column_to_rownames(var = "taxon") |> 
    filter(!is.na(MaxH), #remove species that do not have all seven traits
           !is.na(MaxLS), 
           !is.na(MeanLA),
           !is.na(MeanLDMC),
           !is.na(MeanLL),
           !is.na(MeanSLA), 
           !is.na(percentC),
           !is.na(percentN)) |> 
    mutate(C_N_ratio = percentC/percentN) |> 
    select(!c(percentC, percentN))
  
  if(nrow(FT_wide >0)){ #only do the following if FT_wide has entries:
  
    #standardise trait values
    std_FT_wide <- standard_trait_matrix(trait_matrix = FT_wide, traitlist = c(colnames(FT_wide)))
    
    #now get the distance between each pair of species
    distmat <- as.matrix(dist(std_FT_wide, method = "euclidean"))
    
    if (p == 1) {
      twosp_dist <- as.data.frame(pairwise_fdist(distmat = distmat, sp_positions = positions_plot))
    }else{
      twosp_dist_temp <- as.data.frame(pairwise_fdist(distmat = distmat, sp_positions = positions_plot))
      twosp_dist <- rbind(twosp_dist, twosp_dist_temp)
    } 
  
  }else { #if the plot has no species with complete traits, do the following:
        twosp_dist_temp <- cbind(ID = IDlist[p], replicate = "no complete traits in plot", euclidean_dist = NA, grouping = NA, nurse = NA, target = NA)
        twosp_dist <- rbind(twosp_dist, twosp_dist_temp)
  }
}
##How many plots were thrown out?
nrow(twosp_dist[which(twosp_dist$replicate %in% c("no complete traits in plot", "not all traits measured in plot")), ])
#10 eish

##add siteinfo
twosp_dist <- as.data.frame(twosp_dist)
twosp_dist$ID <- as.numeric(twosp_dist$ID)
#import siteinfo
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  select(ID,COU,SITE,SITE_ID,PLOT,GRAZ, ARIDITY.v3) |> 
  distinct() |> 
  filter(!is.na(ID))

#do the join
twosp_dist <- twosp_dist |> 
  filter(!is.na(euclidean_dist), !euclidean_dist == "NaN") |> 
  inner_join(siteinfo, by = "ID") 

#save the output
write.csv(twosp_dist, "Functional trait data\\results\\Functional_distances_between_2sp_traits_varying.csv")

####Models of dist ~ association####
twosp_dist <- read.csv("Functional trait data\\results\\Functional_distances_between_2sp_traits_varying.csv", row.names = 1) 
twosp_dist$GRAZ <- as.factor(twosp_dist$GRAZ)
twosp_dist$SITE_ID <- as.factor(twosp_dist$SITE_ID)
twosp_dist$grouping <- as.factor(twosp_dist$grouping)
twosp_dist$arid_sq <- (twosp_dist$ARIDITY.v3)^2

###Lets join the results of the CHi2 tests to twosp-dist###
ass <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\Chisq_results_6Feb2024.csv", row.names = 1) |> 
  select(ID, species, association) |> 
  rename(target = species)
#remember that these associations were calculated were calculated at the plot scale. Eg in a specific plot, a species has a significant association with nurse microsites

dist_ass_join <- twosp_dist |> 
  left_join(ass, by = c("target", "ID")) |> 
  filter(association %in% c("nurse", "bare")) #only work with these associations
dist_ass_join$association <- as.factor(dist_ass_join$association)
dist_ass_join$nurse <- as.factor(dist_ass_join$nurse)
dist_ass_join$SITE_ID <- as.factor(dist_ass_join$SITE_ID)
dist_ass_join$ID <- as.factor(dist_ass_join$ID)
dist_ass_join$euclidean_dist <- as.numeric(dist_ass_join$euclidean_dist)

#how many of each association?
dist_ass_join |> 
  group_by(association) |> 
  summarise(obs = n()) #618 bare, 1791 nurse
#how many species?
length(unique(dist_ass_join$target)) #71, many species are lost when we filter out the other associations
length(unique(dist_ass_join$nurse)) #45
length(unique(dist_ass_join$ID)) #41

##does the distance between nurses and bare associated species differ from the distance between nurses and nurse associated species?
dist_ass_null <- glmmTMB(euclidean_dist ~ 1 + (1|nurse) + (1|SITE_ID/ID), data = dist_ass_join)

dist_ass_mod <- glmmTMB(euclidean_dist ~ association + (1|nurse) + (1|SITE_ID/ID), data = dist_ass_join)

summary(dist_ass_mod)
Anova(dist_ass_mod)
anova(dist_ass_null, dist_ass_mod) #p = <0.001

emmeans(dist_ass_mod, specs = "association")
cld(glht(model = dist_ass_mod, mcp(association = "Tukey")))

#model diagnostics
simres <- simulateResiduals(dist_ass_mod)
plot(simres)
#a little underdispersed
#HOV violated

##t tests to see if traits are different from the nurse
dist_nurse_ass_test <- t.test(dist_ass_join[which(dist_ass_join$association == "nurse") , ]$euclidean_dist, 
                        mu = 0, alternative = "greater")

dist_bare_ass_test <- t.test(dist_ass_join[which(dist_ass_join$association == "bare") , ]$euclidean_dist, 
                              mu = 0, alternative = "greater")


###One-dimensional (trait) difference between species####
#This is just the difference in (unstandardised) trait values between species

traitlist <- c("MaxH","MaxLS","MeanLA","MeanLDMC","MeanLL","MeanSLA","C_N_ratio")
sp_positions <- read.csv("Functional trait data\\results\\sp_positions.csv", row.names = 1) 

FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species_graz_conserved.csv", row.names = 1)

complete_FT_mean <- FT |> #get the mean trait value for each sp accross the whole dataset. we will use this to fill missing data values
  group_by(taxon,trait) |> 
  summarise(mean_value = mean(value))

IDlist <- c(unique(sp_positions$ID))

for(p in 1:length(IDlist)) {
  #isolate one plot
  positions_plot <- sp_positions[which(sp_positions$ID == IDlist[p]) , ]
  FT_plot <- FT[which(FT$ID == IDlist[p]) , ]
  
  #Get the mean of each trait for each sp
  FT_mean <- FT_plot |> 
    group_by(taxon,trait) |> 
    summarise(mean_value = mean(value))
  
  #Find out if traits are missing:
  splist <- c(unique(FT_mean$taxon)) #species in plot
  traits_required <- c("MaxH","MaxLS","MeanLA","MeanLDMC","MeanLL","MeanSLA","percentC", "percentN") #traits that should be there
  
  for (s in 1:length(splist)) {
    #traits present
    traits_present <- FT_mean[which(FT_mean$taxon == splist[s]) , ]$trait
    #missing traits
    traits_tofill <- c(traits_required[which(is.na(match(traits_required, traits_present)))])
    
    #if there are missing traits:
    if(length(traits_tofill) > 0) {
      for(f in 1:length(traits_tofill)) {
        #get the filler value from complete_FT_mean
        filler_value <- complete_FT_mean |> 
          filter(taxon == splist[s], trait == traits_tofill[f])
        #if a filler exists, rbind it to FT_mean:
        if(nrow(filler_value) >0) {
          filler_row = data.frame(taxon = splist[s], trait = traits_tofill[f], mean_value = filler_value$mean_value)
          FT_mean <- rbind(FT_mean, filler_row)
          #if there is no filler availible, add an NA
        } else { 
          filler_row = data.frame(taxon = splist[s], trait = traits_tofill[f], mean_value = NA)
          FT_mean <- rbind(FT_mean, filler_row)
        }
      }}} #now FT_mean is as filled as possible
    
    #Get it into wide format
    FT_wide <- FT_mean |>
      pivot_wider(names_from = trait, values_from = mean_value) |> 
      column_to_rownames(var = "taxon") |> 
      mutate(C_N_ratio = percentC/percentN) |> 
      select(!c(percentC, percentN))
  
      for(t in 1:length(traitlist)) {
        #isolate one trait from the sp x trait matrix
        #remove sp with NA values for the trait
        one_trait <- FT_wide[, which(colnames(FT_wide) == traitlist[t])]
        names(one_trait) <- rownames(FT_wide)
        one_trait <- one_trait[which(!is.na(one_trait))]
        
        if(length(one_trait) >0) { #only get the differences if there are actually trait values
          
        #get the difference in trait values between every two species
        trait_diff_matrix <- outer(one_trait, one_trait, FUN = "-")
        
          if(t == 1 & p == 1) { #only for the first run of the loop
          
            trait_diff <- as.data.frame(pairwise_fdist(distmat = trait_diff_matrix, sp_positions = positions_plot)) #add the positions/associations of species
            trait_diff$trait <- traitlist[t]
          
            } else {
              temp_trait_diff <- as.data.frame(pairwise_fdist(distmat = trait_diff_matrix, sp_positions = positions_plot)) #add the positions/associations of species
              temp_trait_diff$trait <- traitlist[t]
            
              trait_diff <- rbind(trait_diff, temp_trait_diff)
            }
        
        #if there are no trait values, just put an NA in
        }else {
          temp_trait_diff <- data.frame(ID = IDlist[p], 
                                        replicate = "No trait values", 
                                        euclidean_dist = NA, grouping = NA, 
                                        nurse = NA, target = NA, trait = traitlist[t])
          trait_diff <- rbind(trait_diff, temp_trait_diff)
        }
  }#end loop through traits
}#end loop through IDlist

#!! there are differnces = 0 because sometimes the dominant species also occurs in the bare or nurse microsite. Thus the dominant and target sp can be the same

#Add ardidty and graz
trait_diff$ID <- as.numeric(trait_diff$ID)
#import siteinfo
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  select(ID,COU,SITE,SITE_ID,PLOT,GRAZ, ARIDITY.v3) |> 
  distinct() |> 
  filter(!is.na(ID))

#do the join
trait_diff <- trait_diff |> 
  rename(trait_difference = euclidean_dist) |> 
  filter(!is.na(trait_difference), !trait_difference == "NaN") |> 
  inner_join(siteinfo, by = "ID") 

write.csv(trait_diff, "Functional trait data\\results\\trait_differences_between_2sp_traits_vary.csv")


###models of trait difference ~association####
trait_fdist <- read.csv("Functional trait data\\results\\trait_differences_between_2sp_traits_vary.csv", row.names = 1)
trait_fdist$SITE_ID <- as.factor(trait_fdist$SITE_ID)
##Lets join the results of the CHi2 tests to sla_fdist###
ass <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis clone\\Facilitation data\\results\\Chisq_results_6Feb2024.csv", row.names = 1) |> 
  select(ID, species, association) |> 
  rename(target = species)
#remember that these associations were calculated were calculated at the plot scale. Eg in a specific plot, a species has a significant association with nurse microsites

trait_ass_join <- trait_fdist |> 
  left_join(ass, by = c("target", "ID")) |> 
  filter(association %in% c("nurse", "bare")) #only work with these associations
trait_ass_join$association <- as.factor(trait_ass_join$association)
trait_ass_join$nurse <- as.factor(trait_ass_join$nurse)
trait_ass_join$SITE_ID <- as.factor(trait_ass_join$SITE_ID)
trait_ass_join$ID <- as.factor(trait_ass_join$ID)


#MaxH model#
maxh_data <- trait_ass_join |> 
  filter(trait == "MaxH") |> 
  mutate(sqrt_trait_difference = sqrt(trait_difference), 
         log_trait_difference = log(trait_difference), 
         neginv_trait_difference = -1/(1+trait_difference))
hist(maxh_data$neginv_trait_difference)
hist(maxh_data$trait_difference)

#how many of each association?
maxh_data |> 
  group_by(association) |> 
  summarise(obs = n()) #751 bare, 1963 nurse

#null model
maxh_null <- glmmTMB(trait_difference ~ 1 + (1|nurse) + (1|SITE_ID/ID), data = maxh_data) #cannot use transformations because of zeroes and negative values
#alternative model
maxh_mod <- glmmTMB(trait_difference ~ association + (1|nurse) + (1|SITE_ID/ID), data = maxh_data)

summary(maxh_mod)
Anova(maxh_mod) #significant effect
anova(maxh_null, maxh_mod) #p = 1.113e-07 ***
emmeans(maxh_mod, specs = "association")
r.squaredGLMM(maxh_mod)

#model diagnostics
simres <- simulateResiduals(maxh_mod)
plot(simres)#underispersed, HOV violated

ggplot(maxh_data, aes(x = association, y = trait_difference)) +
  geom_boxplot()

##Are the means of each association group different from 0?
#Do Whelch's t.test, which does not assume equal variances
maxh_test_nurse <- t.test(maxh_data[which(maxh_data$association == "nurse") , ]$trait_difference, 
                          mu = 0, alternative = "two.sided")
#p <0.001, true mean greater than 0

maxh_test_bare <- t.test(maxh_data[which(maxh_data$association == "bare") , ]$trait_difference, 
                          mu = 0, alternative = "two.sided")
#p <0.001, true mean greater than 0


#MaxLS#
maxls_data <- trait_ass_join |> 
  filter(trait == "MaxLS") |> 
  mutate(sqrt_trait_difference = sqrt(trait_difference), 
         log_trait_difference = log(trait_difference), 
         neginv_trait_difference = -1/(1+trait_difference))
hist(maxls_data$trait_difference)
hist(maxls_data$neginv_trait_difference)
hist(maxls_data$sqrt_trait_difference)

#how many of each association?
maxls_data |> 
  group_by(association) |> 
  summarise(obs = n()) #751 bare, 1963 nurse

#null model
maxls_null <- glmmTMB(sqrt_trait_difference ~ 1 + (1|nurse) + (1|SITE_ID/ID), data = maxls_data) #cannot use transformations because of 0 and - values
#alternative model
maxls_mod <- glmmTMB(sqrt_trait_difference ~ association + (1|nurse) + (1|SITE_ID/ID), data = maxls_data)
#models do not converge for untransformed response, forced to use sqrt despite NA values

summary(maxls_mod)
Anova(maxls_mod) 
anova(maxls_null, maxls_mod) #p = 0.6987
emmeans(maxls_mod, specs = "association") #use emmeans for the post hoc tests and report those results in table!
r.squaredGLMM(maxls_mod)

#model diagnostics
maxls_simres <- simulateResiduals(maxls_mod)
plot(maxls_simres)#unerdispersed, HOV violated

##Are the means of each association group different from 0?
#Do Whelch's t.test, which does not assume equal variances
maxls_test_nurse <- t.test(maxls_data[which(maxls_data$association == "nurse") , ]$trait_difference, 
                          mu = 0, alternative = "two.sided")
#p <0.001, true mean greater than 0

maxls_test_bare <- t.test(maxls_data[which(maxls_data$association == "bare") , ]$trait_difference, 
                         mu = 0, alternative = "two.sided")
#p <0.001, true mean greater than 0


#MeanLA#
meanla_data <- trait_ass_join |> 
  filter(trait == "MeanLA") |> 
  mutate(sqrt_trait_difference = sqrt(trait_difference), 
         log_trait_difference = log(trait_difference),
         neginv_trait_difference = -1/(1+trait_difference))
hist(meanla_data$trait_difference)
hist(meanla_data$neginv_trait_difference)

#null model
meanla_null <- glmmTMB(trait_difference ~ 1 + (1|nurse) + (1|SITE_ID/ID), data = meanla_data) #cannot use transformed respinse because of 0 and - values
#alternative model
meanla_mod <- glmmTMB(trait_difference ~ association + (1|nurse) + (1|SITE_ID/ID), data = meanla_data)

summary(meanla_mod)
Anova(meanla_mod) 
anova(meanla_null, meanla_mod) #0.4486
emmeans(meanla_mod, specs = "association")
r.squaredGLMM(meanla_mod)

#model diagnostics
meanla_simres <- simulateResiduals(meanla_mod)
plot(meanla_simres)#underdispersed, HOV violated

##Are the means of each association group different from 0?
#Do Whelch's t.test, which does not assume equal variances
meanla_test_nurse <- t.test(meanla_data[which(meanla_data$association == "nurse") , ]$trait_difference, 
                          mu = 0, alternative = "two.sided")
#p = 0.6672, mean equal to 0

meanla_test_bare <- t.test(meanla_data[which(meanla_data$association == "bare") , ]$trait_difference, 
                         mu = 0, alternative = "two.sided")
#p <0.001, true mean less than 0


#MeanLDMC#
meanldmc_data <- trait_ass_join |> 
  filter(trait == "MeanLDMC") |> 
  mutate(sqrt_trait_difference = sqrt(trait_difference), 
         log_trait_difference = log(trait_difference), 
         neginv_trait_difference = -1/(1+trait_difference))
hist(meanldmc_data$trait_difference)

#null model
meanldmc_null <- glmmTMB(trait_difference ~ 1 + (1|nurse) + (1|SITE_ID/ID), data = meanldmc_data)
#alternative model
meanldmc_mod <- glmmTMB(trait_difference ~ association + (1|nurse) + (1|SITE_ID/ID), data = meanldmc_data)

summary(meanldmc_mod)
Anova(meanldmc_mod) #significant effect
anova(meanldmc_null, meanldmc_mod) #p < 2.2e-16 ***
emmeans(meanldmc_mod, specs = "association")
r.squaredGLMM(meanldmc_mod)

#model diagnostics
meanldmc_simres <- simulateResiduals(meanldmc_mod)
plot(meanldmc_simres)#residuals normal, HOV violated

ggplot(meanldmc_data, aes(x = association, y = trait_difference))+
  geom_boxplot()

##Are the means of each association group different from 0?
#Do Whelch's t.test, which does not assume equal variances
meanldmc_test_nurse <- t.test(meanldmc_data[which(meanldmc_data$association == "nurse") , ]$trait_difference, 
                            mu = 0, alternative = "two.sided")
#p <0.001, mean greater than 0

meanldmc_test_bare <- t.test(meanldmc_data[which(meanldmc_data$association == "bare") , ]$trait_difference, 
                           mu = 0, alternative = "two.sided")
#p <0.001, mean greater than 0


#MeanLL#
meanll_data <- trait_ass_join |> 
  filter(trait == "MeanLL") |> 
  mutate(sqrt_trait_difference = sqrt(trait_difference), 
         log_trait_difference = log(trait_difference), 
         neginv_trait_difference = -1/(1+trait_difference))
hist(meanll_data$trait_difference)
hist(meanll_data$neginv_trait_difference)

#null model
meanll_null <- glmmTMB(trait_difference ~ 1 + (1|nurse) + (1|SITE_ID/ID), data = meanll_data) #cannot use transfromations because of 0 and - values
#alternative model
meanll_mod <- glmmTMB(trait_difference ~ association + (1|nurse) + (1|SITE_ID/ID), data = meanll_data)

summary(meanll_mod)
Anova(meanll_mod) #no significant effect
anova(meanll_null, meanll_mod) #p = 0.5138
emmeans(meanll_mod, specs = "association")
r.squaredGLMM(meanll_mod)

#model diagnostics
meanll_simres <- simulateResiduals(meanll_mod)
plot(meanll_simres)#a little overdispersed, HOV violated

##Are the means of each association group different from 0?
#Do Whelch's t.test, which does not assume equal variances
meanll_test_nurse <- t.test(meanll_data[which(meanll_data$association == "nurse") , ]$trait_difference, 
                              mu = 0, alternative = "two.sided")
#p <0.001, mean less than 0

meanll_test_bare <- t.test(meanll_data[which(meanll_data$association == "bare") , ]$trait_difference, 
                             mu = 0, alternative = "two.sided")
#p <0.001, mean greater than 0


#MeanSLA#
meansla_data <- trait_ass_join |> 
  filter(trait == "MeanSLA") |> 
  mutate(sqrt_trait_difference = sqrt(trait_difference), 
         log_trait_difference = log(trait_difference), 
         neginv_trait_difference = -1/(1+trait_difference))
hist(meansla_data$trait_difference)
hist(meansla_data$neginv_trait_difference)

#null model
meansla_null <- glmmTMB(trait_difference ~ 1 + (1|nurse) + (1|SITE_ID), data = meansla_data)
#alternative model
meansla_mod <- glmmTMB(trait_difference ~ association + (1|nurse) + (1|SITE_ID), data = meansla_data)
#model doesnt converge with site_ID/ID RE

summary(meansla_mod)
Anova(meansla_mod) #
anova(meansla_null, meansla_mod) #p = 0.0002708 ***
emmeans(meansla_mod, specs = "association")
r.squaredGLMM(meansla_null)

#model diagnostics
meansla_simres <- simulateResiduals(meansla_mod)
plot(meansla_simres)#a little underdispersed, HOV violated

ggplot(meansla_data, aes(x = association, y = trait_difference))+
  geom_boxplot()

##Are the means of each association group different from 0?
#Do Whelch's t.test, which does not assume equal variances
meansla_test_nurse <- t.test(meansla_data[which(meansla_data$association == "nurse") , ]$trait_difference, 
                            mu = 0, alternative = "two.sided")
#p <0.001, mean less than 0

meansla_test_bare <- t.test(meansla_data[which(meansla_data$association == "bare") , ]$trait_difference, 
                           mu = 0, alternative = "two.sided")
#p <0.001, mean less than 0


#C_N_ratio#
cn_data <- trait_ass_join |> 
  filter(trait == "C_N_ratio") |> 
  mutate(sqrt_trait_difference = sqrt(trait_difference), 
         log_trait_difference = log(trait_difference), 
         neginv_trait_difference = -1/(1+trait_difference))
hist(cn_data$trait_difference)

#null model
cn_null <- glmmTMB(trait_difference ~ 1 + (1|nurse) + (1|SITE_ID/ID), data = cn_data)
#alternative model
cn_mod <- glmmTMB(trait_difference ~ association + (1|nurse) + (1|SITE_ID/ID), data = cn_data)

summary(cn_mod)
Anova(cn_mod) 
anova(cn_null, cn_mod) #p = 9.832e-05 ***
emmeans(cn_mod, specs = "association")
r.squaredGLMM(cn_mod)

#model diagnostics
cn_simres <- simulateResiduals(cn_mod)
plot(cn_simres)#resiuals normal enough, HOV ok

ggplot(cn_data, aes(x = association, y = trait_difference))+
  geom_boxplot()


##Are the means of each association group different from 0?
#Do Whelch's t.test, which does not assume equal variances
cn_test_nurse <- t.test(cn_data[which(cn_data$association == "nurse") , ]$trait_difference, 
                             mu = 0, alternative = "two.sided")
#p <0.001, mean less than 0

cn_test_bare <- t.test(cn_data[which(cn_data$association == "bare") , ]$trait_difference, 
                            mu = 0, alternative = "two.sided")
#p <0.001, mean greater than 0
