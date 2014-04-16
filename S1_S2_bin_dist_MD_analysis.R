# bin dist analysis more dimensions




rm(list = ls())
setwd("~/Desktop/Repeat evolution/bin_synteny")


# Species workin on

spec1 <- "HUMAN"
spec2 <- "COW"


# if I could get this to work on everything i have saved then it will be much easier to run through everything
# the probelem is that I have to worry about all the pca stuff too


# load everything
# ensure there are the correct species in each variable

S1 <- read.table("./PCA/bins_hg19_pca_regions.txt")
S2 <- read.table("./PCA/bins_bt6_pca_order.txt")
s.bin <- read.table(paste("S_bins/" , spec1, "aligning", spec2, sep = ""), header = TRUE)
load("~/Documents/phd/Initial_PCA_results/HG19/HG19.pca")
S1.p <- pca
load("~/Documents/phd/Initial_PCA_results/BT6/pca.bt6.RObject")
S2.p <- pca


#####
#This works because the bin names are supposadly in the same order as their PCA values
# with future data this may not be the case

# merge bins adn coordinates into a single data frame 
S1.bin.coord <- data.frame(S1[,5], S1.p$x)
colnames(S1.bin.coord) <- c("S1.bin", paste("S1",colnames(S1.p$x), sep = ""))
S2.bin.coord <- data.frame(S2[,5], S2.p$x)
colnames(S2.bin.coord) <- c("S2.bin", paste("S2",colnames(S2.p$x), sep = ""))
bin.coord <- merge(merge(s.bin,S1.bin.coord), S2.bin.coord)




no <- 5
S1bin <- unique(bin.coord$S1.bin)
S1.N.dist <- NULL
S1ploters <- unique(bin.coord[,c(grep("S1.bin", colnames(bin.coord)), grep("S1PC", colnames(bin.coord)))])

for( i in seq(S1bin)){
	S1point <- S1ploters[S1ploters[,1] == S1bin[i],]
	
	# here is where we get the eucledian distance for all the points
	# the dataframe plottters contains all the information we are interested in 
	# neighbors <- (data.frame(ploters[,1],sqrt((ploters[,2] - point[,2])^2 + (ploters[,3] - point[,3])^2)))
		
	# the new way first gets x - y for each dimension and squeres it 
	# then it adds em up and finds the square root
	bigN <- NULL
	for( z in seq(length(S1point)-1)){
		N <- (S1point[,z + 1] - S1ploters[,z + 1])^2
		bigN <- cbind(bigN,N)
	}
	neighbors <- data.frame(S1ploters[,1] ,sqrt(rowSums(bigN)))

	# here we reorder based on the lowest number so as to identify the closest bins
	neighbors <- neighbors[order(neighbors[,2]),][2:(no+1),]
	neighbors[,3:4] <- c(rep(S1bin[i], no), rep(mean(neighbors[,2]), no))
	neighbors <- data.frame(S1.bin = neighbors[,3], neighbor = neighbors[,1],ave.dist = neighbors[,4], dist = neighbors[,2])
	S1.N.dist <- rbind(S1.N.dist, neighbors)
}

# we now have a table that has the distance to the 5 nearest neighbors of each point in human

# how the hell do we extend it to cow

#all the points i need for each calculation hvae been effectivly subsetted

# I got this wroking for multiple dimensions however it doesn't normailse correctly 
# I still don't know how to get it working for multi mappers

#((d(a1,b1)))

# this way we do each bin dist by itself rather than all at once and calculate the mean at the end


# here we can make a matrix after making the dist calculations


# calculating the distance for each point in S2
# fill in the distances then multiply by the proportions
#    b1 b2 b3  
# a1
# a2
#
# get the sum of the matrix and divide by the sum(pa) * sum(pb)






Dists <- NULL
for(i in seq(S1bin)){
	
	
	A_points <- bin.coord[bin.coord$S1.bin == S1bin[i],]
	# find the neigbors
	other_points <- bin.coord[bin.coord$S1.bin %in% S1.N.dist$neighbor[S1.N.dist$S1.bin == S1bin[i]],]
	
	#use the H.bin to get at each individual B point
	Bs <- unique(other_points$S1.bin)
	S2.distforall <- NULL
	#run a loop through each B point to get the distance of it
	for(bi in Bs){
		B_points <- other_points[other_points$S1.bin == bi,]
		Dist.matrix <- NULL
		dist_o1 <- NULL
		for(o in seq(dim(A_points)[1])){
			# get the distances for A1 to each B
			bigN <- NULL
			for( z in grep("S2PC",colnames(B_points))){
				N <- (B_points[,z] - A_points[o,z])^2
				bigN <- cbind(bigN,N)
				}
			DAtoB <- sqrt(rowSums(bigN))
			Dist.matrix <- rbind(Dist.matrix, DAtoB)
			}
		Dist.matrix <- matrix(data = Dist.matrix, nrow = dim(A_points)[1], ncol = dim(B_points)[1])
		A_p <- A_points$S2.Proportion
		B_p <- B_points$S2.Proportion
		S2.distfor1 <- sum(t(Dist.matrix*A_points$S2.Proportion)*B_points$S2.Proportion) / (sum(A_points$S1.Proportion)*sum(B_points$S1.Proportion))
		S2.distforall <- rbind(S2.distforall,S2.distfor1)
		}
	S2.mean.dist <- mean(S2.distforall)
	
	# get a DF with cols that show 
	row_final <- data.frame(S1.bin = S1bin[i],  S2.mean.dist , S1.mean.dist = S1.N.dist$ave.dist[S1.N.dist$S1.bin == S1bin[i]][1], S1PC1 = A_points$S1PC1[1], S1PC2 = A_points$S1PC2[1])
	Dists <- rbind(Dists, row_final)

	
	}
	

Dists$percentage.diff <- sqrt((((Dists$S1.mean.dist - Dists$S2.mean.dist)/ Dists$S1.mean.dist) * 100)^2)
Dists$diff <- sqrt((((Dists$S1.mean.dist - Dists$S2.mean.dist)))^2)
Dists$S1_to_S2_ratio <- Dists$S1.mean.dist / Dists$S2.mean.dist
Dists$S2_to_S1_ratio <- Dists$S2.mean.dist / Dists$S1.mean.dist
	
write.table(Dists, file = paste("Div_score_pc1_pc2/", spec1, "_to_", spec2, "_Mdim_dist", sep = ""), quote = FALSE, sep = "\t", row.names= FALSE )
	
	
	
	
library(ggplot2)
	
	qplot(Dists$S1PC1, Dists$S1PC2, colour = Dists$percentage.diff) + scale_colour_gradient(limits=c(min(Dists$percentage.diff), 200), high = "blue", low = "yellow")
	
	
	
	
	
	
	
	
	