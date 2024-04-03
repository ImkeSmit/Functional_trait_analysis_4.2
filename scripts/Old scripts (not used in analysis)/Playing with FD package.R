####Playing with the FD package
#install.packages("FD")
library(FD)
library(dplyr)
library(glmmTMB)
library(stats)

wd <- "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait data"
setwd(wd)

#This data contains species without trait values
trait_dat_complete <- read.csv("FT_match_facilitation_plots_plotspecific_species.csv", row.names = 1)
#remove sp without trait values
trait_dat <- trait_dat_complete[which(!is.na(trait_dat_complete$MeanSLA)) , ] 


#How many species do we have data for in each plot?
sp_perplot <- as.data.frame(tapply(trait_dat$sp_fullname, trait_dat$ID, FUN = length))
colnames(sp_perplot) <- "nsp"
sp_perplot$ID <- row.names(sp_perplot)
row.names(sp_perplot) <- c(1:nrow(sp_perplot))
##There are many plots with very few sp. We will need to decide which plots have enough sp to be able to run analysis. 
##To get community trait value, the sp for which we have traits will have to make up 80% of total cover. 
##So we will have to go back to the cover data (from the graphs and figures script) and see which plots fulfil this criterion. 
##But also.... That doesn't account for cover of the nurse
##Then we can try to gap fill plots that do not have enough sp. 

####EXPLORE ANALYSIS WITH SA DATA####
ft_sa <- trait_dat[which(trait_dat$Country == "South Africa") , ]

###Calculate FD for each plot###
#We will need to make a separate sp x trait matrix for each plot, since we have traits at the plot level.
#we will put each x matrix in a list, and then let dbFD run through each list element

IDlist <- c(unique(ft_sa$ID))
 #list to put the species x trait matrices in 

for (i in 1:length(IDlist)) {
  plot <- ft_sa[which(ft_sa$ID == IDlist[i]) , ]
  
  #make a species by trait matrix,x, for each plot
  x <- plot[ , c(which(colnames(plot) %in% c("sp_fullname", "MeanLL", 
                                               "MeanSLA", "MeanLDMC", "MeanLA", "MaxH", "MaxLS")))]
  row.names(x) <- x$sp_fullname #make the rownames the species names
  x <- x[, -which(colnames(x) == "sp_fullname")]
  
  xlist[[i]] <- x
}

###Now calculate FD for each plot
#Table to put results in 
column_headers <- c("ID", "nsp", "FRic", "qual.FRic", "FEve", "FDiv", "FDis", "RaoQ")
FD_results <- data.frame(matrix(nrow = length(IDlist), ncol = length(column_headers)))
colnames(FD_results) <- column_headers

for (k in 1:length(xlist)) {
  result <- dbFD(x = xlist[[k]], #for now, do not include weighting matrix
                 stand.x = T,    #standardise trait to mean 0 and unit variance
                 corr = "cailliez",  #if distance matrix cannot be represented in euclidean space, take the sqrt of the distance
                 calc.FRic = T, stand.FRic = F, #do not standardise FRic to the global FRic (which includes all sp)
                 m = 2,          #the number of PCoA axes to use when calculating FRic
                 #for 6 traits we could make m= 6??
                 calc.FDiv = T)
  
  FD_results[k , ]$ID <- IDlist[k]
  FD_results[k , ]$nsp <- result$nbsp
  FD_results[k , ]$FRic <- result$FRic
  FD_results[k , ]$qual.FRic <- result$qual.FRic
  FD_results[k , ]$FEve <- result$FEve
  FD_results[k , ]$FDiv <- result$FDiv
  FD_results[k , ]$FDis <- result$FDis
  FD_results[k , ]$RaoQ <- result$RaoQ
}
##Match Aridity and Graz with FD_results
FD_results <- merge(x = FD_results, y = filter(trait_dat[, c(1,15,16)], !duplicated(ID)), by = "ID", all.x = FALSE, all.y =FALSE)
###maxH and MaxLS is missing from Monsonia camdeboensis in ID = 248
###maybe we can borrow values from elsewhere or give it the avg?
###What does the function do with missing values?

###Does FD influence the plotlevel NIntc####
#Import NIntc results
nint_result <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation data\\results\\NIntc_results_allcountries_26Sep.csv", row.names = 1)
#Isolate SA plots
nint_sa <- nint_result[which(nint_result$country == "southafrica"),]
#Get the mean NIntc_richness per plot
nint_sa <- nint_sa[which(!is.na(nint_sa$NInta_richness)) , ] #remove NA values
nint_perplot <- as.data.frame(tapply(nint_sa$NIntc_richness, nint_sa$ID, FUN = mean))
colnames(nint_perplot) <- "avg_NIntc_richness"
nint_perplot$ID <- row.names(nint_perplot)
row.names(nint_perplot) <- c(1:nrow(nint_perplot))

#merge nint_perplot with FD_results
sa_dat <- merge(FD_results, nint_perplot, by = "ID")
sa_dat$GRAZ <- as.factor((sa_dat$GRAZ))
plot(sa_dat$avg_NIntc_richness ~ sa_dat$FRic)
plot(FRic ~ ARIDITY.v3, data = sa_dat)
plot(FRic ~ GRAZ, data = sa_dat)

###GLMM's
#make nint binomial
sa_dat$avg_NIntc_richness_binom <- (sa_dat$avg_NIntc_richness + 1)/2
mod1 <- glmmTMB(avg_NIntc_richness_binom ~ FRic, data = sa_dat, family = "binomial")
summary(mod1)


####How does the functional match change?####
##First step is to define functional groups
#the first argument will have to be the average trait for each species, here I just use a complete trait x species matrix from the list
clusters <- hclust(dist(xlist[[1]], method = "euclidean"), method = "ward.D")
plot(clusters)

#after looking at the plot, 4 clusters seem to be present
groups <- as.data.frame(cutree(clusters, k = 4))
colnames(groups) <- "fgr"
groups$spname <- rownames(groups)
rownames(groups) <- c(1:nrow(groups))

##now get the CWM of each fgr



####CODE THAT MAY BE USEFUL LATER####
#First explore the functions with the dummy dataset
traits <- dummy$trait
abun <- dummy$abun

test <- dbFD(x = traits, a = abun,
             #matrix of weights not included, all sp weighted equally
             stand.x = TRUE, #standardise all numeric traits to mean 0 and unit variance
             corr = "sqrt", #correction to use when spXsp distance matrix cannot be repsented in euclidean space. default is sqrt
             calc.FRic = TRUE, m = "max", #the number of PCoA axes to keep as ‘traits’ for calculating FRic. 
             #min is the min number of traits  number of traits that allows the species ≥ 2^traits condition to be met
             stand.FRic = FALSE,
             scale.RaoQ = FALSE, 
             calc.FGR = TRUE,
             clust.type = "ward",
             calc.CWM = TRUE, CWM.type = "all", 
             calc.FDiv = TRUE, dist.bin = 2,
             print.pco = TRUE, #return pcoa and eigenvalues
             messages = TRUE)

#Lets plot the PC axes
PC <- test$x.axes
plot(PC$A1, PC$A2, xlab = "PC1", ylab = "PC2")

##For now, let's work with plots that have >20 sp
p <- sp_perplot[which(sp_perplot$nsp == max(sp_perplot$nsp)) , ]$ID #get the plot with the highest number of species
iran <- trait_dat[which(trait_dat$ID == p) , ] #it is an iranian plot

#transform iran to a species by trait matrix
iran <- iran[ , c(which(colnames(iran) %in% c("sp_fullname", "MeanLL", 
                                                   "MeanSLA", "MeanLDMC", "MeanLA", "MaxH", "MaxLS")))]
row.names(iran) <- iran$sp_fullname
iran <- iran[, -which(colnames(iran) == "sp_fullname")]

test2 <- dbFD(iran, #matrix of abundances not included - all sp assumed to have equal abundance
            #matrix of weights also not included, all sp weighted equally
      stand.x = TRUE, #standardise all numeric traits to mean 0 and unit variance
     corr = "sqrt", #correction to use when spXsp distance matrix cannot be repsented in euclidean space. default is sqrt
     calc.FRic = TRUE, m = "min", #the number of PCoA axes to keep as ‘traits’ for calculating FRic. 
                                  #min is the min number of traits  number of traits that allows the species ≥ 2^traits condition to be met
     stand.FRic = FALSE,
     scale.RaoQ = FALSE, 
     calc.FGR = TRUE,
     clust.type = "ward",
     calc.CWM = TRUE, CWM.type = "all", 
     calc.FDiv = TRUE, dist.bin = 2,
     print.pco = TRUE, #return pcoa and eigenvalues
     messages = TRUE)

##Lets plot the species on the PC axes
#retreive principal components
PC_iran <- test2$x.axes
PC_iran$sp <- rownames(PC_iran)
rownames(PC_iran) <- c(1:nrow(PC_iran))

#add functional groups
groups <- as.data.frame(test2$spfgr)
colnames(groups) <- "functional_group"
groups$sp <- rownames(groups)
rownames(groups) <- c(1:nrow(groups))

merged <- merge(PC_iran, groups, by = "sp")

#now plot
plot(merged$A1, merged$A2, xlab = "PC1", ylab = "PC2", col = merged$functional_group)

plot(merged$A1, merged$A3, xlab = "PC1", ylab = "PC3", col = merged$functional_group)

