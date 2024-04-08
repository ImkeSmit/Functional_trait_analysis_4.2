###RDA:###
###Do the response traits of associated communities depend on the effect traits of dominant plants###
library(tidyverse)
library(tidylog)
library(vegan)

##Import raw countries data
data_files <- list.files("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Countriesv3")
countrynames <- c("algeria", "argentina", "australia", "chile", "chinachong", "chinaxin", "iranabedi", "iranfarzam", 
                  "israel", "namibiablaum", "namibiawang", "southafrica",  "spainmaestre", "spainrey")
for(i in 1:length(data_files)) {                              
  assign(paste0(countrynames[i]),                                   
         read.csv2(paste0("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\Countriesv3\\",
                          data_files[i])))
}

##lets only work with Southafrica for now
#only work with nurse microsites
sa <- southafrica |> 
  filter(Microsite == 2, 
         !is.na(Species.within.quadrat))


##Import trait data
#From the filled trait data for plotspecific species
##!remember that these species were not necessarily filled from the same graz level
FT <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species.csv",
               row.names = 1) |> 
  filter(COU == "South Africa")

###Loop to make sp x trait matrices for each nurse microsite####

#names of the traits collectd
traits_collected <- c(unique(FT$trait))
#plot Id's
IDlist <- c(unique(sa$ID))
#empty list of matrices
matlist <- list()

l = 1
for (i in 1:length(IDlist)) {
  
  #isolate a plot in the facilitation and FT data
  fac_plot <- sa[which(sa$ID == IDlist[i]) , ]
  FT_plot <- FT[which(FT$ID == IDlist[i]) , ]
  
  #list of reps in this plot
  replist <- c(unique(fac_plot$Number.of.replicate))
  
  #get the species that occur in the nurse microsite of every rep
  for (r in 1:length(replist)) {
    one_rep <- fac_plot[which(fac_plot$Number.of.replicate == replist[r]) , ]
    comm_sp <- c(one_rep$Species.within.quadrat) #species names
    
    #create a vector of ID_rep to identify each matrix. we will use this to name the list at the end of the loop
    if(l == 1) {comm_index <- paste(IDlist[i], replist[r], sep = "_")}
    else{comm_index_temp <- paste(IDlist[i], replist[r], sep = "_")
        comm_index <- c(comm_index, comm_index_temp)}
    
    #create an empty sp x trait matrix
    comm_mat <- matrix(nrow = length(comm_sp), ncol = length(traits_collected))
    row.names(comm_mat) <- comm_sp
    colnames(comm_mat) <- traits_collected
    
      #get the trait vals of species in comm_mat
      for(s in 1:length(comm_sp)) { #for each species
        for(t in 1:length(traits_collected)) { #and each trait
        
        #isolate the trait value(s) for each species and trait  
        val <- FT_plot |> 
          filter(taxon == comm_sp[s], 
                 trait == traits_collected[t]) |> 
          select(value) 
      
        #if there are no values for this trait, put NA in the matrix
        if(nrow(val) == 0) {comm_mat[s,t] <- NA} 
          #if there are 2 or more trait measurements, get the mean and put that in the matrix
          else if (nrow(val) > 1) {mean_val <- mean(val$value)
            comm_mat[s,t] <- mean_val} 
            #if there is only on trait value, put that in the matrix
            else {comm_mat[s,t] <- val$value} 
        }#end loop through traits
      }#end loop through species
    
      matlist[[l]] <- comm_mat #put filled matrix in the list
  
      l = l+1
  }#end loop through reps
}#end loop through plots
names(matlist) <- comm_index #name the list elements


###Now we need to get the mean trait values for each community####
comm_means <- data.frame(ID_rep = comm_index, mean_C_N_ratio = NA, meanLL = NA, meanSLA = NA, meanLDMC = NA, meanLA = NA, 
                         mean_H = NA, mean_LS = NA)

for (m in 1:length(matlist)) {
  comm <- matlist[[m]]
  ID_rep <- names(matlist[m])
  
  comm_means[m, 1] <- ID_rep
  comm_means[m, 2] <- mean(c(comm[, 2]/comm[, 1]), na.rm = T)
  comm_means[m, 3] <- mean(comm[, 3], na.rm = T)
  comm_means[m, 4] <- mean(comm[, 4], na.rm = T)
  comm_means[m, 5] <- mean(comm[, 5], na.rm = T)
  comm_means[m, 6] <- mean(comm[, 6], na.rm = T)
  comm_means[m, 7] <- mean(comm[, 7], na.rm = T)
  comm_means[m, 8] <- mean(comm[, 8], na.rm = T)
}


###Now we need to get the explanatory variables together####
#names of the traits collectd
traits_collected <- c(unique(FT$trait))
#plot Id's
IDlist <- c(unique(sa$ID))

exp_var <- data.frame(ID_rep = NA, nurse_sp = NA, aridity = NA, graz = NA, mean_C_N_ratio = NA, meanLL = NA, 
                      meanSLA = NA, meanLDMC = NA, meanLA = NA, mean_H = NA, mean_LS = NA)




l = 1
for (i in 1:length(IDlist)) {
  
  #isolate a plot in the facilitation and FT data
  fac_plot <- sa[which(sa$ID == IDlist[i]) , ]
  FT_plot <- FT[which(FT$ID == IDlist[i]) , ]
  
  #list of reps in this plot
  replist <- c(unique(fac_plot$Number.of.replicate))
  
  #get the nurse species of each rep
  for (r in 1:length(replist)) {
    one_rep <- fac_plot[which(fac_plot$Number.of.replicate == replist[r]) , ]
    nurse_sp <- c(one_rep$ID_Microsite) 
    
    #concatenate ID and repto make an identifier
    comm_index <- paste(IDlist[i], replist[r], sep = "_")
    
    exp_var[l,1] <- comm_index
    exp_var[l,2] <- nurse_sp
    
    #get the traits of that nurse
    for(t in 1:length(traits_collected)) { #and each trait
      
      #isolate the trait value(s) for each species and trait  
      val <- FT_plot |> 
        filter(taxon == nurse_sp, 
               trait == traits_collected[t]) |> 
        select(value) 
      
      #if there are no values for this trait, put NA in the matrix
      if(nrow(val) == 0) {exp_var[l,t] <- NA} 
      #if there are 2 or more trait measurements, get the mean and put that in the matrix
      else if (nrow(val) > 1) {mean_val <- mean(val$value)
      comm_mat[s,t] <- mean_val} 
      #if there is only on trait value, put that in the matrix
      else {comm_mat[s,t] <- val$value} 
    }#end loop through traits
    
    
  }
}




