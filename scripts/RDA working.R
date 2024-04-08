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

traits_collected <- c(unique(FT$trait))

#Now we need to assemble species by trait matrices for every nurse microsite in southafrica
IDlist <- c(unique(sa$ID))
matlist <- list()

l = 1
for (i in 1:length(IDlist)) {
  one_plot <- sa[which(sa$ID == IDlist[i]) , ]
  
  replist <- c(unique(one_plot$Number.of.replicate))
  
  for (r in 1:length(replist)) {
    one_rep <- one_plot[which(one_plot$Number.of.replicate == replist[r]) , ]
    
    comm_sp <- c(one_rep$Species.within.quadrat)
    
    if(l == 1) {comm_index <- paste(IDlist[i], replist[r], sep = "_")}
    else{comm_index_temp <- paste(IDlist[i], replist[r], sep = "_")
        comm_index <- c(comm_index, comm_index_temp)}
    
    comm_mat <- matrix(nrow = length(comm_sp), ncol = length(traits_collected))
    row.names(comm_mat) <- comm_sp
    colnames(comm_mat) <- traits_collected
    
    matlist[[l]] <- comm_mat

    l = l+1
  }
}
names(matlist) <- comm_index
