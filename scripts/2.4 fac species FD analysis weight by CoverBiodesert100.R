###FUNCTIONAL DIVERSITY ANALYSIS###
#Get the functional richness of each plot with FD
#Use FT match facilitation plots plotspecific species

library(tidyverse)
library(tidylog)
library(traitstrap)
library(FD)
library(glmmTMB)
library(car)
library(DescTools)
library(corrplot)

##Import filled FT data for the facilitation plots
FT_raw <- read.csv("Functional trait data\\Clean data\\FT_filled_match_facilitation_plots_plotspecific_species_graz_conserved.csv", row.names = 1)

FT <- FT_raw|> 
  mutate(taxon = str_replace(taxon, " ", "_")) |> 
  pivot_wider(names_from = trait, values_from = value) |> 
  #calculate the C:N ratio
  mutate(C_N_ratio = percentC/percentN) |> 
  select(!c(percentC, percentN)) |> 
  pivot_longer(cols = c("MeanLL","MeanSLA","MeanLDMC","MeanLA","MaxH","MaxLS","C_N_ratio"),
               names_to = "trait", values_to = "value")

#get species in the FT data
FT_sp <- FT |> 
  distinct(taxon)


#make a cover matrix using coverBiodesert100
#remember coverBiodesert100_FT is the cover of the"filler"
covermat <- FT |> 
  select(c(ID, taxon, coverBiodesert100)) |> 
  distinct() |> 
  filter(!is.na(coverBiodesert100)) |> 
  pivot_wider(names_from = taxon, values_from = coverBiodesert100) |> 
  column_to_rownames(var = "ID")
#open records are species for which we do not have trait data

#replace NA's in covermat with 0
for (r in 1:nrow(covermat)) {
  for(c in 1:ncol(covermat)) {
    
    record <- covermat[r , c]
    
    if(is.na(record)) {
      covermat[r , c] <- 0
    }
  }
}
sum(covermat[ , 1])


###Create separate species x trait matrices, and sp x plot cover matrices for each plot, put them in a list####
#get the mean trait value where there are multiple sp entries in a plot
FT_mean <- FT |> 
  filter(!is.na(value)) |> 
  group_by(ID, taxon, trait) |> 
  summarise(mean_value = mean(value)) |> 
  ungroup()

IDlist <- c(unique(FT_mean$ID))

#list of trait matrices
xlist <- vector(mode='list', length = length(IDlist))
names(xlist) <- IDlist

#list of cover matrices
alist <- vector(mode='list', length = length(IDlist))
names(alist) <- IDlist

for(i in 1:length(IDlist)) {
  plot_wide <- FT_mean |> 
    filter(ID == IDlist[i]) |> 
    select(!ID) |> 
    pivot_wider(names_from = trait, values_from = mean_value) |> 
    column_to_rownames(var = "taxon")
  
  FT_plot_names <- c(row.names(plot_wide))
  
  xlist[[i]] <- plot_wide
  
  #subset covermat for the  ID and the sp in plot_wide
  #because dbfd needs the sp between the traits and the community to correspond exactly
  covermat_sub <- covermat |> 
    filter(row.names(covermat) == IDlist[i]) |> 
    select(any_of(FT_plot_names))
  
  alist[[i]] <- covermat_sub
}


###Quality assesment - do we have enough species?####
#compare the species and cover in xlist to the species and cover in the quadrat data
quad <- read.csv("Functional trait data\\Clean data\\quadrat_survey_all_plots.csv", row.names = 1)

IDlist <- c(row.names(covermat))

nsp_dat <- data.frame(ID = IDlist , trait_nsp = c(rep(NA, length(IDlist))) , 
                      quad_nsp = c(rep(NA, length(IDlist))), 
                      trait_percent_cover = c(rep(NA, length(IDlist))), 
                      quad_percent_cover = c(rep(NA, length(IDlist))))
l = 1
for (t in 1:length(IDlist)) {
  
  ##number of species in FT data
  trait_plot <- xlist[[name = IDlist[t]]]
  trait_nsp <- nrow(trait_plot)
  
  ##number of sp in quad data
  quad_nsp <- quad |> #use quad before we made it wide and removed extra species
    filter(ID == IDlist[t]) |> 
    distinct(taxon) |> 
    nrow()
  
  ##added cover of all species in FT data
  trait_cov <- FT |> 
    filter(ID == IDlist[t]) |> 
    distinct(taxon, coverBiodesert100) |> 
    summarise(trait_percent_cover = sum(coverBiodesert100))
  
  ##added cover of all species in quad data
  quad_cov <- quad |> #use quad before we made it wide and removed extra species
    filter(ID == IDlist[t]) |> 
    group_by(quadrat) |> 
    mutate(sum_quadrat_cover = sum(percent_cover)) |> 
    ungroup() |>
    distinct(quadrat, sum_quadrat_cover) |> 
    summarise(sum_plot_cover = sum(sum_quadrat_cover), max_quadrats = max(quadrat), 
              plot_percent_cover = sum_plot_cover/max_quadrats)
  
  nsp_dat[l, 2] <- trait_nsp
  nsp_dat[l, 3] <- quad_nsp
  nsp_dat[l, 4] <- trait_cov
  nsp_dat[l, 5] <- quad_cov[, 3]
  
  l = l+1
}

nsp_dat <- nsp_dat |> 
  mutate(percent_sp_in_FT = (trait_nsp/quad_nsp)*100, 
         trait_cover_div_quad_cover = (trait_percent_cover/quad_percent_cover)*100)


#plots where we have data for more than 80% cover
plots_for_FD <- nsp_dat |> 
  filter(trait_cover_div_quad_cover >= 80) |> 
  distinct(ID) #resulting in 75 plots

##However to make this comparable with the full version of this analysis (with FT-match facilitation plots), we need to run it for exactly the same plots
full_FD_IDlist <- read.csv("Functional trait data\\results\\FD_results_20May2024.csv") |> 
  distinct(ID)
#subset for the plots in full_FD_IDlist
xlist_red <- xlist[which(names(xlist) %in% c(as.character(full_FD_IDlist$ID)))]

alist_red <- alist[which(names(alist) %in% c(as.character(full_FD_IDlist$ID)))]


##plot 115 breaks  FD because Ferula sink has many NA's. Let's remove this species.
new_115 <- xlist_red[[which(names(xlist_red) == "115")]] |> 
  filter(!row.names(xlist_red[[which(names(xlist_red) == "115")]]) == "Ferula_sinkiangensis")
xlist_red[[which(names(xlist_red) == "115")]] <- new_115

new_cov_115 <- alist_red[[which(names(alist_red) == "115")]] |> 
  select(!Ferula_sinkiangensis)
alist_red[[which(names(alist_red) == "115")]]  <- new_cov_115


###Now calculate FD for each plot
#Table to put results in 
column_headers <- c("ID", "nsp", "FRic", "qual.FRic", "FEve", "FDiv", "FDis", "RaoQ")
FD_results <- data.frame(matrix(nrow = nrow(full_FD_IDlist), ncol = length(column_headers)))
colnames(FD_results) <- column_headers

IDlist <- c(full_FD_IDlist$ID)

for (k in 1:length(xlist_red)) {
  
  nsp <- nrow(xlist_red[[k]])
  
  if(nsp > 1) {
    
    result <- dbFD(x = xlist_red[[k]], #trait matrix
                   a = alist_red[[k]], #community cover matrix
                   w.abun = T,     #weight metrics by relative abundance of species
                   stand.x = T,    #standardise trait to mean 0 and unit variance
                   corr = "cailliez",  #if distance matrix cannot be represented in euclidean space, take the sqrt of the distance
                   calc.FRic = T, stand.FRic = F, #do not standardise FRic to the global FRic (which includes all sp)
                   m = 2,          #the number of PCoA axes to use when calculating FRic
                   calc.FDiv = T, 
                   print.pco = T)  #return eigenvalues and PCoA axes
    
    FD_results[k , ]$ID <- IDlist[k]
    FD_results[k , ]$nsp <- result$nbsp
    FD_results[k , ]$FRic <- result$FRic
    FD_results[k , ]$qual.FRic <- result$qual.FRic
    FD_results[k , ]$FEve <- result$FEve
    FD_results[k , ]$FDiv <- result$FDiv
    FD_results[k , ]$FDis <- result$FDis
    FD_results[k , ]$RaoQ <- result$RaoQ
    
  }else {
    FD_results[k , ]$ID <- IDlist[k]
    FD_results[k , c(2:8)] <- "only one sp"
  }
}


write.csv(FD_results, "Functional trait data\\results\\FD_results_20May2024_for_facilitation_species.csv")
