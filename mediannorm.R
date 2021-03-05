mediannorm2 <- function(df){  ## creates function mediannorm2 ##
	
	label1  <- which(df$Immunolabel == immunolabel1)
	Pillar   <- which(df$Immunolabel == "Pillar")
	Modiolar <- which(df$Immunolabel == "Modiolar")
	label2     <- which(df$Immunolabel == immunolabel2)
	Nucleus  <- which(df$Immunolabel == "Nucleus")

	
	## calculate the median volume for the pre and postsynaptic label##
	## Then divides all values by their respective median ##
	meddf <-aggregate(Volume~Immunolabel, data = df, FUN=median, na.rm=TRUE)
	medlab1 <- meddf$Volume[which(meddf$Immunolabel == immunolabel1)]
 	medlab2 <- meddf$Volume[which(meddf$Immunolabel == immunolabel2)]

	df$Volume[label1] <- (df$Volume[label1]/medlab1)
	df$Volume[label2] <- (df$Volume[label2]/medlab2)
	df <- df	

	return(df)
}