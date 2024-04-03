###Cleaning the Quadrat data
library(tidyverse)
library(tidylog)
library(readxl)
library(DescTools)

#import quadrat data for which we have trait data
#Get it into a usable format
#subset for the facilitation plots

FT_allsites <- read.csv("Functional trait data\\FT_all_sites.csv", row.names = 1) 
#countries in the FT data (spelling changed to match spelling in quadrat data)
FT_countries <- c("Algeria", "Argentina", "Australia", "Botswana","Brazil","Canada", "Chile", "China","Ecuador","Hungria", "Iran", 
                   "Israel","Kazajstan","Kenia","Mexico","Mongolia", "Namibia","Niger","Palestina","Peru","Portugal", "South Africa",
                   "Spain","Tunez","USA")
length(FT_countries)


#The quadrat data has a sheet for every country, lets extract only the countries for which we have FT data
# Specify the file path to your Excel file
excel_file <- "C:\\Users\\imke6\\Documents\\Msc Projek\\Functional trait analysis\\Functional trait data\\Raw data\\Cover quadrats_98sites_16_12_19.xlsx"

# Create a list to store the data frames for each selected sheet
selected_sheets_data <- list()

# Loop through each sheet and extract the data
for (sheet_name in FT_countries) {
  sheet_data <- read_excel(excel_file, sheet = sheet_name)
  selected_sheets_data[[sheet_name]] <- sheet_data
}# Now, selected_sheets_data contains data frames for the fac countries
length(selected_sheets_data)
names(selected_sheets_data)


##Each country's sheet contains multiple species x site matrices, one for each plot.
#each matrix is separated by a row of NA's. We need to cut the dataframe at each row of NA's 
# Function to cut dataframe at each row with NA in the first column
cut_dataframe <- function(df) {
  # Find the row indices where the first column is NA
  cut_indices <- which(is.na(df[ , 1]))
  
  # Split the dataframe at the identified indices
  split_dataframes <- split(df, cumsum(seq_along(df[ , 1]) %in% cut_indices))
  
  # Remove empty dataframes (resulting from consecutive NA rows)
  split_dataframes <- split_dataframes[sapply(split_dataframes, function(x) !all(is.na(x[ , 1])))]
  
  return(split_dataframes)
}

##Now cut the dataframes of each country
for (i in 1:length(selected_sheets_data)) {
  country <- as.data.frame(selected_sheets_data[[i]]) #select one country from the list
  country <- rbind(colnames(country), country) #change the column names to the first row
  colnames(country) <- c("Plotname", seq(1:(ncol(country)-1))) #give new column names
  if (i == 1) {
    country_cut <- cut_dataframe(country)
  } else {
    temp_country_cut <- cut_dataframe(country)
    country_cut <- c(country_cut, temp_country_cut)
  }
} #country_cut is a list of dataframes. each dataframe is a species x site matrix of a plot. It contains the matrices of all countries
length(country_cut) #326


##Get the names of all the plots in country_cut
cut_plotnames <- data.frame(plotnames = c(rep(NA, 326)))
for (t in 1:length(country_cut)) {
  plot <- country_cut[[t]]
  
  if(is.na(plot[1,1])) {#remove the first rowif it has NA's
    plot <- plot[-1 , ] 
  }
  
  plotname <- plot[1,1]
  cut_plotnames[t,1] <- plotname
} 




###Get all the quadrats in long format and rbind them to each other
for (t in 1:length(country_cut)) {
  focus <- as.data.frame(country_cut[[t]])
  colnames(focus)[1] <- "taxon" #change name of first column
  
  if(is.na(focus[1,1])) {#remove the first rowif it has NA's
    focus <- focus[-1 , ] 
  }
  
  focus$site_plot <- focus[1,1] #add a variable with the site name and plot
  focus <- focus[-1 , ] #remove the first row, it just contains the site name and quadrat names
  
  #remove the last two columns in each matrix
  if(focus$site_plot[1] %like% "%Verdelecho%" == FALSE) {
    focus <- focus[ , -which(colnames(focus) %in% c("101", "102"))] #remove rows called 101 and 102. 101 has only NA's, 102 has the species cover totals
  } else {
    focus <- focus[ , -which(colnames(focus) %in% c("121", "122"))] #Verdelecho has 120 quadrats surveyed, so different rows need to be removed
  }
  
  
  if(t == 1) {
    focus_long <- focus |>  #change to long format
      pivot_longer(cols = !c(taxon, site_plot), names_to = "quadrat", values_to = "percent_cover")
    
  } else {
    
    temp_focus_long <- focus |>  #change to long format
      pivot_longer(cols = !c(taxon, site_plot), names_to = "quadrat", values_to = "percent_cover")
    
    focus_long <- bind_rows(focus_long, temp_focus_long)
  }
}
length(unique(focus_long$site_plot)) #must be same length as country_cut

#Which plot is missing?
long_plotnames <- data.frame(plotname = c(unique(focus_long$site_plot)))
match(long_plotnames$plotname, cut_plotnames$plotname)
match(cut_plotnames$plotname, long_plotnames$plotname)
cut_plotnames$plotnames[202]
which(long_plotnames$plotname == "Natab1")
#Natab1 is missing from focus_long because it has no species


#Problems with NA values:
#For all sites except Verdelecho, 100 quadrats were surveyed.
#In these sites, all quadrats after 100 can be removed
#In Verdelecho, 120 quadrats were surveyed. All plots after 120 can be removed
#Site Lac du Bois has many NA's. Looks like they didnt do all the transects
#Contenda1 only did 80 transects
#Ciempozuelos, Monfrague, San MArtin, ZaragozaArido, ZaragozaSemiarido, Monsul, Puertodelascoberteras has NA's after transect 122, there are probably spaces in the excel file
#LondonFarm Syferkuil, Kruger and Mara only did 20 transects

#remove NA records
focus_long_working <- focus_long |> 
  filter(!is.na(percent_cover)) |> #remove NA records
  mutate(quadrat = as.numeric(quadrat))

#remove extra quadrats where they shouldn't be
for(t in 1:length(nrow(focus_long_working))) {  
if(focus_long_working[t, which(colnames(focus_long_working) == "site_plot")] %like% "%Verdelecho%" == FALSE) {
  focus_long_working <- focus_long_working |> 
    filter(quadrat < 101)
  }
}


###Now add plotinfo to focus_long####
#Import the information about the biodesert sites
siteinfo <- read.csv("C:\\Users\\imke6\\Documents\\Msc Projek\\Facilitation analysis\\Facilitation data\\BIODESERT_sites_information.csv") |> 
  select(ID,COU,SITE,SITE_ID,PLOT,GRAZ, ARIDITY.v3) |> 
  distinct() |> 
  mutate(SITE = str_to_lower(SITE), 
         SITE = str_replace_all(SITE, " ", "")) |> 
  mutate(site_plot = str_c(SITE, PLOT, sep = "")) |> 
  filter(!is.na(ID))


##Clean up the plot names in focus_long-working so that we can use them in joins
quadrat_survey <- focus_long_working |> 
  mutate(site_plot = str_to_lower(site_plot)) |> 
  #split the site_plot into 4 separate strings
  separate_wider_delim(site_plot, delim = " ", names = c("split1", "split2", "split3",  "split4", "split5"), 
                       too_few = "align_end") |> 
  #we need to eliminate cases where the string starts with a number
  #so set it to NA if split1 or split2 is a number
  select(!split1) |>  #split1 has only numbers so delete it
  #replace the numbers in the first position with NA
  mutate(split2 = replace(split2, split2 %in% c(NA, "83","84","85","86","87","88","93") , " "), 
         split3 = replace(split3, split3 %in% c(NA, "", "94", "95", "96", "97", "98", "99", "100", "101", "102", "103","104","105","106"), " ")) |> 
  
  mutate(split4 = replace(split4, split4 %like% c("%nyngan%"), "nyngam"),
         split4 = replace(split4, is.na(split4), " "),
         split4 = replace(split4, split4 %like% c("%wuxin%"), "wuxing")) |> 
  
  mutate(split5 = replace(split5, split5 == "kruger1", "krugerpark1"),
         split5 = replace(split5, split5 == "kruger2", "krugerpark2"),
         split5 = replace(split5, split5 ==  "kruger3", "krugerpark3"),
         
         split5 = replace(split5, split5 == "moztaza1", "mostaza1"),
         split5 = replace(split5, split5 == "moztaza2", "mostaza2"),
         split5 = replace(split5, split5 == "moztaza3", "mostaza3"),
         
         split5 = replace(split5, split5 == "mara1", "maraexperimentalfarm1"),
         split5 = replace(split5, split5 == "mara2", "maraexperimentalfarm2"),
         split5 = replace(split5, split5 == "mara3", "maraexperimentalfarm3"),
         split5 = replace(split5, split5 == "mara4", "maraexperimentalfarm4"),
         
         split5 = replace(split5, split5 %like% c("%sandvel3%"), "sandveld3")) |>  #change spelling so that it corresponds with fac_plotinfo
  
  mutate(site_plot = str_c(split3, split4, split5, sep = "")) |> #concatenate to make site_plot
  mutate(site_plot = str_trim(site_plot)) |>  #remove space in front of site_plot 
  select(!c(split2, split3, split4, split5)) |>  #delete these columns 
  inner_join(siteinfo, by = "site_plot")  #an inner join only keeps observations in x that have a match in y.
#dubois, detalca,  is only in the quadrat data, no mention of it in siteinfo
#these plots are thus removed from quadrat_survey


###Fix names in quadrat data####
clean_quadrat_survey <- quadrat_survey |> 
  #remove infraspecies ranks
  separate_wider_delim(taxon, delim = " ",   
                       names = c("split1", "split2", "split3",  "split4"), too_few = "align_start") |> 
  mutate(split2 = case_when(is.na(split2) ~ " ", .default = as.character(split2))) |> 
  mutate(taxon = str_c(split1, split2, sep = " ")) |> 
  select(!c(split1,split2,split3, split4)) |> 
  mutate(taxon = str_squish(taxon), #remove spaces before and after string
         taxon = str_to_sentence(taxon))  #make sentence case

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

data_harmony <- clean_quadrat_survey

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
head(change_tracker, 50)
change_tracker <- change_tracker |> 
  filter(!is.na(new_spec))

##Export the cleaned quadrat data for all plots
write.csv(data_harmony, "Functional trait data\\quadrat_survey_all_plots.csv")