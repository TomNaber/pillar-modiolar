getCouples <- function(df, distances, resultsdf2, Nimlabel1)  ## create the function "getcouples#
{ 
	label1  <- which(df$Immunolabel == immunolabel1) ## imports the immunolablel names from the main script ##
 	label2  <- which(df$Immunolabel == immunolabel2) ## imports the immunolablel names from the main script ##
	
 	## creates a matrix with the number of presynaptic markers as number of rows ##
 	## and the number of postsynaptic markers as the number columns ## 
 	## the matrix is assigned values of 0 to start with ##
 	
	nRows <-  dim(df[label1,])[1] 
	nCols <-  dim(df[label2,])[1] 
	tmp <- matrix(0, nRows, nCols)
	
	## Replace 0s in matrix with euclidean distances between each combination of pre and postsynaptic marker##
	## Subsequently pick the pair with the lowest euclid distance ##
	## This distance is filled in into a new dataframe "distances" ##
	## The pre and postsynaptic markers are then removed from the matrix and the next shortest distances is calculated etc.. ##
	
	for (irow in 1 : nRows)
	{
		tmp[irow, ] <- (((df$PosX[label1])[irow] - df$PosX[label2]) ^ 2 + 
			((df$PosY[label1])[irow] - df$PosY[label2]) ^ 2 + 
			((df$PosZ[label1])[irow] - df$PosZ[label2]) ^ 2) ^ (1/2)
	}

	lab1 <- df$PointID[label1][which(min(tmp) == tmp, arr.ind = TRUE)[1]]
	lab2 <- df$PointID[label2][which(min(tmp) == tmp, arr.ind = TRUE)[2]]
	lab1vol <- df$Volume[label1][which(min(tmp) == tmp, arr.ind = TRUE)[1]]
	lab2vol <- df$Volume[label2][which(min(tmp) == tmp, arr.ind = TRUE)[2]]

	distances <- rbind(distances, data.frame(
		lab1 = lab1,
		lab2 = lab2,
		dist = min(tmp, na.rm = TRUE),
		lab1vol = lab1vol,
		lab2vol = lab2vol)
	)
	newDat = (df[label1,])[df$PointID[label1] != lab1,]
	newDat = rbind(newDat, (df[label2,])[df$PointID[label2] != lab2, ])

	if ((dim(distances)[1] == Nimlabel1) | (dim(newDat)[1] < 1) | (length(newDat$PointID[label2]) <= 1))
	{	
		distances$ratio <- distances$lab1vol/distances$lab2vol
		distances$PairID <-1: (dim(distances)[1])		

		distances2 <- data.frame(distances$lab2, distances$lab1, distances$dist,
						 distances$lab1vol,distances$lab2vol,distances$ratio,distances$PairID) 
		colnames(distances2) <- colnames(distances)
		distances <- rbind(distances, distances2)
		distances$fileN <- rep(file, length(dim(distances)[1]))
		colnames(distances) <- c("PointID", "Pairedw", "Distance", "lab1vol", "lab2vol", "Ratio","PairID", "fileN")

		resultsdf2 <- rbind(resultsdf2, distances)

		return(resultsdf2)
	} else {
		getCouples(newDat, distances,resultsdf2, Nimlabel1)
	}
}

