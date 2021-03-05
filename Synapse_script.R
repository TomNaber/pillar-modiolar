#############################################################################################################
# This script was written for R version 4.0.0; last updated: 21.FEB.2021                                    #
#############################################################################################################

## runs the script from the designated folder (change the directory into your working directory) ##
  setwd("C:/Users/luscu/Desktop/Reijntjes et al., 2020 STAR Protocols")
# setwd('C:/Users/Daniël/Dropbox/daniel_paolo') another example of a working directory ##
# setwd('/home/paolot/Dropbox/MATLAB_Functions/daniel/') another example of a working directory ##

################################################ Set up ####################################################

rm(list=ls()) ## cleans the variables in the global environment ##
#getwd()       ## shows you the current working directory ##

#############################################################################################################
##Load required libraries (you need to install these packages before you can load them f.e. "install.pacakges("readxl") ##
library(readxl)
library(rgl)

#############################################################################################################
## Create a list of files that includes each file in the working directory ##
## with the .xlsx extension (so all excel files) ##
files <- dir(pattern= '*.xlsx')

## running "files" shows you the filenames in the list of files so you can doublecheck ##
## if they are included correctly ##
files

#############################################################################################################
## Set immunolabels. Here you must specify the classification of the immunolabels as used in the excel file ##
## Under the heading "Immunolabel" see example excel files. N.B. If this is incorrect the analysis will not run  ##
immunolabel1 <-"Ribbons"   
immunolabel2 <-"PSD"     

#############################################################################################################
figures = FALSE #Set figures to TRUE if you want to provide 3D figures of the plotted synapses with the plotted plane 

#############################################################################################################
##Creates required empty dataframes. These will later be used to form one dataframe called results ##
## Therefore, these dataframes are only used to run the script and do not provide data individually#

distances <- data.frame(lab1 = character(),
				lab2 = character(),
				dist = numeric(),
				lab1vol = numeric(),
				lab2vol = numeric()
)

resultsdf1 <- data.frame(PosX = numeric(),
			     	 PosY = numeric(),
			       PosZ = numeric(),
				 Immunolabel = character(),
				 Volume = numeric(),
				 PointID = character(),
				 Genotype = character(),
				 Expgroup = character(),
				 FileN = character()
)

resultsdf2 <- data.frame(PointID = character(),
				Pairedw = character(),
				Distance = numeric(),
				lab1vol = numeric(),
				lab2vol = numeric(),
				Ratio = numeric(),
				PairID = character(),
				fileN = character()
)

resultsdf3PCA <- data.frame(volTypePCA = character(),
				    fileN = character(),
				    PointID = character()
)

################################ Import functions ##########################################################
## Function #1: Calculate euclidean distances (for every ribbon/GluA pair (or other markers)) ##
source('get_couples.R')

## function 2: normalize the volume data for Ribbons and GluA to their respective medians ##

source('mediannorm.R')

## function 3: Plotting 3d representation and separate points based on PCAplane ##
source('separate_points.R')

## Run functions on files ##

for (file in files)
{	
	df <- readxl::read_xlsx(file, sheet =1, col_names =  TRUE) # Creates dataframe "df" by reading the excel file in "files" # 
	df<- mediannorm2(df) # runs the mediannorm script which normalizes all the volumesin "df" to the median volume #

	df$fileN <- rep(file, length(dim(df)[1])) # Creates an extra column in the dataframe df that ##
	## specifies the filename that the data is from  ##
	resultsdf1 <- rbind(resultsdf1,df) ## exports the normalized data to the temporary dataframe resultsdf1 ##
	##  and concatenates all data from all files  ##
	
	imlabel1  <- which(df$Immunolabel == immunolabel1) ## creates a length vector for the number of presynaptic marker counts in file ##
	## This vector is used to determine how many iterations of the "getcouples" function to run#
	resultsdf2 <- getCouples(df, distances, resultsdf2, dim(df[imlabel1,])[1]) # runs the getcouples function#
	# This function uses the euclidean distance between synaptic markers to form pairs between the closest pre and postsynaptic markers#
	
	print(file) ## check to see that the script is running and indicates progress ##
	resultsdf3PCA <- separate_pointsPCA(df, resultsdf3PCA, figures = figures, PCA = TRUE )
	## runse the separate points PCA script that calculates a plane that defines the pillar and modiolar ##
	## Classification of the synaptic markers##

}	

#############################################################################################################

## Remove data with distances over 1 micron ##

testlength<-length(which(resultsdf2$Distance > 1.0))

if(testlength > 1){ remove1 <- which(resultsdf2$Distance > 1)
resultsdf2 <- resultsdf2[-remove1,] }


## Merge temporary dataframes ##
results <- merge(resultsdf1, resultsdf2, c("fileN", "PointID"),all.x = TRUE, all.y = TRUE)
results <- merge(results, resultsdf3PCA, c("fileN", "PointID"),all.x = TRUE, all.y = TRUE)


## Remove unwanted data from results ##
rem <- c(which(results$Immunolabel == "Pillar"), which(results$Immunolabel == "Modiolar"), 
which(results$Immunolabel == "Nucleus"), which(results$Immunolabel == "Plane"),
which(results$Immunolabel == "HCbutt"),which(results$Immunolabel == "AltHCbutt"),which(results$Immunolabel == "Frame"))

results <- results[-rem,]
#############################################################################################################

## export the results as an excel file ##

write.csv(results,"Results.csv", row.names=F)
