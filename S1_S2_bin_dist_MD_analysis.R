# bin dist analysis more dimensions




rm(list = ls())
setwd("~/Desktop/Repeat evolution/bin_synteny")

# load everything

C <- read.table("./PCA/bins_hg19_pca_regions.txt")
#change
H <- read.table("./PCA/bins_bt6_pca_order.txt")
# change
s.bin <- read.table("S_bins/C_bin_aligning_H_bin_independant1", header = TRUE)
#change
load("~/Documents/phd/Initial_PCA_results/BT6/pca.bt6.RObject")
H.p <- pca
load("~/Documents/phd/Initial_PCA_results/HG19/HG19.pca")
pca$x[,1] <- pca$x[,1] * -1
B.p <- pca



# merge PC2 scores with appropriate bins 
Hb <- data.frame(H[,5], H.p$x)
Hname <- paste("H",colnames(H.p$x), sep = "")
colnames(Hb) <- c("H.bin", Hname)
Cb <- data.frame(C[,5], B.p$x)
Cname <- paste("C",colnames(B.p$x), sep = "")
colnames(Cb) <- c("C.bin", Cname)

Q <- merge(s.bin,Cb)
Q2 <- merge(Q, Hb)
Q2 <- data.frame(Q2[,1:4], Q2[,27:46], Q2[,5:26])
bin <- Q2


#####
#  #   In this part here i can definatly organsie something to measure it
#

# for each human bin find its x nearest neighbors
# and calculate the mean eucleidan distance 

# then calculate the mean eucledian distance for those same neighbors in cow

# findong neares neighbors 5 is our number for now
no <- 5
Hbin <- unique(bin[,1])
H.N.dist <- NULL
ploters <- unique(bin[c(1, grep("HPC", colnames(bin)))])

for( i in seq(Hbin)){
	point <- ploters[ploters[,1] == Hbin[i],]
	
	# here is where we get the eucledian distance for all the points
	# the dataframe plottters contains all the information we are interested in 
	# neighbors <- (data.frame(ploters[,1],sqrt((ploters[,2] - point[,2])^2 + (ploters[,3] - point[,3])^2)))
		
	# the new way first gets x - y for each dimension and squeres it 
	# then it adds em up and finds the square root
	bigN <- NULL
	for( z in seq(length(point)-1)){
		N <- (point[,z + 1] - ploters[,z + 1])^2
		bigN <- cbind(bigN,N)
	}
	neighbors <- data.frame(ploters[,1] ,sqrt(rowSums(bigN)))

	# here we reorder based on the lowest number so as to identify the closest bins
	neighbors <- neighbors[order(neighbors[,2]),][2:(no+1),]
	neighbors[,3:4] <- c(rep(Hbin[i], no), rep(mean(neighbors[,2]), no))
	neighbors <- data.frame(Hbin = neighbors[,3], neighbor = neighbors[,1],ave.dist = neighbors[,4], dist = neighbors[,2])
	H.N.dist <- rbind(H.N.dist, neighbors)
}

# we now have a table that has the distance to the 5 nearest neighbors of each point in human

# how the hell do we extend it to cow

#all the points i need for each calculation hvae been effectivly subsetted

# I got this wroking for multiple dimensions however it doesn't normailse correctly 
# I still don't know how to get it working for multi mappers

#((d(a1,b1)))

# this way we do each bin dist by itself rather than all at once and calculate the mean at the end




Dists <- NULL
for(i in seq(Hbin)){
	
	
	A_points <- bin[bin[,1] == Hbin[i],]
	other_points <- bin[bin[,1] %in% H.N.dist[H.N.dist[,1] == Hbin[i],2], ]
	
	#use the H.bin to get at each individual B point
	Bs <- unique(other_points$H.bin)
	cow.distforall <- NULL
	#run a loop through each B point to get the distance of it
	for(bi in Bs){
		B_points <- other_points[other_points$H.bin == bi,]
		
		dist_o1 <- NULL
		for(o in seq(dim(A_points)[1])){
			# get the distances for A1 to each B
			bigN <- NULL
			for( z in grep("CPC",colnames(B_points))){
				N <- (B_points[,z] - A_points[o,z])^2
				bigN <- cbind(bigN,N)
				}
			DAtoB <- sqrt(rowSums(bigN))
			# now have all the distances from A1 to B
			# next multuply them by there respective pB
			DAtoBpb <- DAtoB *  B_points$C.P
			# now we get all of the a1 to b distances and sum them, then we multiply by pA1
			allA1toB <- sum(DAtoBpb) * A_points[o,4]
		
			dist_o1 <- rbind(dist_o1, allA1toB)	
			}
		A_points$norm.dist <- dist_o1	
		cow.distfor1 = (sum(A_points$norm.dist)/ (sum(A_points$H.P) *  sum(B_points$H.P)))	
		cow.distforall <- rbind(cow.distforall,cow.distfor1)	
			
		}
	cow.mean.dist <- mean(cow.distforall)
	o1_final <- data.frame(H.bin = Hbin[i],  cow.mean.dist , human.mean.dist = H.N.dist[H.N.dist[,1] == Hbin[i],][1,3], HPC1 = A_points$HPC1[1], HPC2 = A_points$HPC2[1])
	Dists <- rbind(Dists, o1_final)

	
	}
	

Dists$percentage.diff <- sqrt((((Dists$human.mean.dist-Dists$cow.mean.dist)/ Dists$human.mean.dist) * 100)^2)
Dists$diff <- sqrt((((Dists$human.mean.dist-Dists$cow.mean.dist)))^2)

	
write.table(Dists, file = "Div_score_pc1_pc2/C_H_pc1_pc2_dists_MDim_new", quote = FALSE, sep = "\t", row.names= FALSE )
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	# for each A point get the distance to each B point
	dist_o1 <- NULL
	for(o in seq(dim(A_points)[1])){
		
		# distance calculation
		# calculating the mean distance in one step here
		#dist <- mean(sqrt((o2_points[,7] - o1_points[o,7])^2 +  (o2_points[,8] - o1_points[o,8])^2)*o2_points$C.P)
		
		# get the distances for A1 to each B
		bigN <- NULL
		for( z in grep("CPC",colnames(B_points))){
			N <- (B_points[,z] - A_points[o,z])^2
			bigN <- cbind(bigN,N)
			}
		DAtoB <- sqrt(rowSums(bigN))	
		# now have all the distances from A1 to B
		# next multuply them by there respective pB
		DAtoBpb <- DAtoB *  B_points$C.P
		# now we get all of the a1 to b distances and sum them, then we multiply by pA1
		allA1toB <- sum(DAtoBpb) * A_points[o,4]
		
		dist_o1 <- rbind(dist_o1, allA1toB)
		}
	A_points$norm.dist <- dist_o1
	
# calculate the last part getting all the a points and dividing by the total human total for cow A and B
# now the next part of it is to turn into to a mean rather than a sum so I think we divide by "no" to do that

	cow.mean.dist = (sum(A_points$norm.dist)/ (sum(A_points$H.P) *  sum(B_points$H.P)))/no
	
	# arrange the dists for cow, take all the primery dists and multiply by their length, then divide the mean by the sum of the human length	
	o1_final <- data.frame(H.bin = Hbin[i],  cow.mean.dist , human.mean.dist = H.N.dist[H.N.dist[,1] == Hbin[i],][1,3], HPC1 = A_points$HPC1[1], HPC2 = A_points$HPC2[1])
	
	
	
	Dists <- rbind(Dists, o1_final)
	}
# calculate the percentage difference petween the human and cow mean distances
# ((human - cow) / human )* 100 

Dists$percentage.diff <- sqrt((((Dists$human.mean.dist-Dists$cow.mean.dist)/ Dists$human.mean.dist) * 100)^2)
Dists$diff <- sqrt((((Dists$human.mean.dist-Dists$cow.mean.dist)))^2)




write.table(Dists, file = "Div_score_pc1_pc2/H_B_pc1_pc2_dists_MDim_new", quote = FALSE, sep = "\t", row.names= FALSE )
