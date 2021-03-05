separate_pointsPCA <- function(df,resultsdf3, figures = FALSE, PCA = FALSE ) ## creates function "separate pointsPCA" ##
{ 
	if (PCA == TRUE){
		dataList <- (na.omit(df$Volume))*7 ## increases the volume of the synaptic markers for better visual representation ##
	
		## Provides labels for all markers that need to be plotted and changes the colour of the plotted markers ##
		label1  <- which(df$Immunolabel == immunolabel1)
		Pillar   <- which(df$Immunolabel == "Pillar")
		Modiolar <- which(df$Immunolabel == "Modiolar")
		label2     <- which(df$Immunolabel == immunolabel2)
		Nucleus  <- which(df$Immunolabel == "Nucleus")

		goodcol <- as.numeric(df$Immunolabel)

		Alist <- c(label1,label2,Nucleus,Pillar,Modiolar)
		Blist <- c(rep(3,length(label1)),rep(2,length(label2)),rep(4,length(Nucleus)),
		rep(8,length(Pillar)),rep(8,length(Modiolar)))

		goodcol <- replace(goodcol, Alist, Blist)

	  ## Calculates the principal components through the XYZ coordinates through the nuclei ##
		## Because the data input is coordinates, the PCA components will indicate directional vectors ##
		## The vector for the first loading is used. This loading explains the variation along the IHC axis ##
		## Since we want to plot the plane along the IHC axis this is ideal ##
		xyz <- df[Nucleus,1:3]
		N <- nrow(xyz)
		mean_xyz <- apply(xyz, 2, mean)
		xyz_pca   <- princomp(xyz) 
		dirVector <- xyz_pca$loadings[, 1]   # PC1
		
		## Calucaltes a line to plot based on the directional vector of PC1 ##

		xyz_fit <- matrix(rep(mean_xyz, each = N), ncol=3) + xyz_pca$score[, 1] %*% t(dirVector) 

		t_ends <- c(min(xyz_pca$score[,1]) - 0.2, max(xyz_pca$score[,1]) + 0.2)  # for both ends of line
		endpts <- rbind(mean_xyz + t_ends[1]*dirVector, mean_xyz + t_ends[2]*dirVector)
		
		## Calculates the principal components through the XYZ coordinates through the synaptic markers ##
		## The rest same as above ##

		n2 <- c(label1 , label2)
		xyz2 <- df[n2,1:3]
		N2 <- nrow(xyz2)
		mean_xyz2 <- apply(xyz2, 2, mean)
		xyz_pca2   <- princomp(xyz2) 
		dirVector2 <- xyz_pca2$loadings[, 1]   # PC1

		xyz_fit2 <- matrix(rep(mean_xyz2, each = N2), ncol=3) + xyz_pca2$score[, 1] %*% t(dirVector2) 

		t_ends2 <- c(min(xyz_pca2$score[,1]) - 0.2, max(xyz_pca2$score[,1]) + 0.2)  # for both ends of line
		endpts2 <- rbind(mean_xyz2 + t_ends2[1]*dirVector2, mean_xyz2 + t_ends2[2]*dirVector2)
	
		## Identify three points (P,Q,R) to use to calculate the vector equation for a plane ##
		
		P <- xyz_fit[(length(xyz_fit[,1])/2),1:3] ## midpoint of the line through the nuclei ##
		Q <- xyz_fit2[1,1:3]                      ## Minimum (first point) of the line through the synaptic markers ##
		R <-xyz_fit2[length(xyz_fit2[,1]),1:3]    ## Maximum (last point) of the line through the synaptic markers ##
	
		## Vector equation of a plane in simplified steps ##
		
		PQ <- Q-P
		PR <- R-P
		i1 <- PQ[2]*PR[3]
		j1 <- PQ[3]*PR[1]
		k1 <- PQ[1]*PR[2]
		i2 <- PQ[3]*PR[2]
		j2 <- PQ[1]*PR[3]
		k2 <- PQ[2]*PR[1]
		a <- i1 - i2
		b <- j1 - j2
		c <- k1 - k2
		d = -(a*P[1] + b*P[2] + c*P[3]) # note the sign change 
		
		## Calculate the distance "-" or "+" distance of all the plotted synaptic markers to the plane ##
		## "dflab1" and "dflab2" identify on which side ("-" or "+" side) of the plane the pillar and modiolar facings are ##
		## synaptic markers are labelled accoringly as pillar or modiolar in a new dataframe ##

		dflab1 <- df[c(label1,Pillar,Modiolar),]
		dflab2 <- df[c(label2,Pillar,Modiolar),]
		
		Mxlab1 <- dflab1$PosX
		Mylab1 <- dflab1$PosY
		Mzlab1 <- dflab1$PosZ

		Mxlab2 <- dflab2$PosX
		Mylab2 <- dflab2$PosY
		Mzlab2 <- dflab2$PosZ

		pointslab1 <-(a*Mxlab1 + b*Mylab1 + c*Mzlab1 + d)/sqrt(a^2+b^2+c^2)
		pointslab2 <- (a*Mxlab2 + b*Mylab2 + c*Mzlab2 + d)/sqrt(a^2+b^2+c^2)

		neglab1 <- which(pointslab1 <0)
		poslab1 <- which(pointslab1 >0)
		neglab2 <- which(pointslab2 <0)
		poslab2 <- which(pointslab2 >0)
	
		neglistlab1ID <- dflab1$PointID[neglab1]
		poslistlab1ID <- dflab1$PointID[poslab1]
		neglistlab2ID <- dflab2$PointID[neglab2]
		poslistlab2ID <- dflab2$PointID[poslab2]

		if (any(dflab1$Immunolabel[neglab1] == "Pillar")) {
			pillarlab1ID <- as.character(neglistlab1ID)
			modiolarlab1ID <- as.character(poslistlab1ID)
		} else {
			modiolarlab1ID <- as.character(neglistlab1ID)
			pillarlab1ID <- as.character(poslistlab1ID)
		}

		if (any(dflab2$Immunolabel[neglab2] == "Pillar")) {
			pillarlab2ID <- as.character(neglistlab2ID)
			modiolarlab2ID <- as.character(poslistlab2ID)
		} else {
			modiolarlab2ID <- as.character(neglistlab2ID)
			pillarlab2ID <- as.character(poslistlab2ID)
		}

    ## creates a dataframe with an ID for each synaptic marker, and a specification for the side ##
		IDlist  <- c(pillarlab1ID, modiolarlab1ID, pillarlab2ID, modiolarlab2ID)
		voltype <- c(rep('Pillar', length(pillarlab1ID)), rep('Modiolar', length(modiolarlab1ID)),
			 rep('Pillar', length(pillarlab2ID)), rep('Modiolar', length(modiolarlab2ID)))

		resultsdf3PCA <- rbind(resultsdf3PCA, 
					data.frame(
					volTypePCA = voltype,
					fileN = rep(file, length(voltype)),
					PointID = IDlist
			 
					)
				)
	
		## If figure is set to TRUE in the master file, this code will create interactive 3D plots ##
		## These plots will indicate all synaptic markers and the plane ##
		## These plots can be used to inspect visually whether the plane has been positioned correctly ##
		
		if( figures == TRUE){
			open3d()
			plot3d(df$PosX, df$PosY, df$PosZ, col= as.numeric(df$Immunolabel),
			type="s", size=0, lit=FALSE)

			for(i in dataList) {
      			points3d(df$PosX[which(dataList == i)],
					df$PosY[which(dataList == i)], 
					df$PosZ[which(dataList == i)],
					col= goodcol[which(dataList == i)], size=i)
			}
	
		abclines3d(mean_xyz, a = dirVector, col="blue", lwd=2)     # mean + t * direction_vector
		for(i in 1:N) segments3d(rbind(xyz[i,], xyz_fit[i,]), col="green3")

		abclines3d(mean_xyz2, a = dirVector2, col="blue", lwd=2)     # mean + t * direction_vector
		for(i in 1:N2) segments3d(rbind(xyz2[i,], xyz_fit2[i,]), col="green3")
	
		planes3d(a,b,c,d)
		} 

		return(resultsdf3PCA)
	}
}
