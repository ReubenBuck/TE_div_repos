# multi species




cow.Dists <- read.table("Div_score_pc1_pc2/HUMAN_to_COW_Mdim_dist", header = TRUE)
sheep.Dists <- read.table("Div_score_pc1_pc2/HUMAN_to_SHEEP_Mdim_dist", header = TRUE)
horse.Dists <- read.table("Div_score_pc1_pc2/HUMAN_to_HORSE_Mdim_dist", header = TRUE)

library(ggplot2)
	
	
	
	
	
	
qplot(cow.Dists$S1PC1* -1 , cow.Dists$S1PC2, colour = cow.Dists$percentage.diff) + scale_colour_gradient(limits=c(min(cow.Dists$percentage.diff), 200), high = "red", low = "white")

qplot(sheep.Dists$S1PC1 * -1, sheep.Dists$S1PC2, colour = sheep.Dists$percentage.diff) + scale_colour_gradient(limits=c(min(sheep.Dists$percentage.diff), 200), high = "red", low = "white")	
qplot(horse.Dists$S1PC1* -1, horse.Dists$S1PC2, colour = horse.Dists$percentage.diff) + scale_colour_gradient(limits=c(min(horse.Dists$percentage.diff), 200), high = "red", low = "white")	

plot(horse.Dists$S1.bin,horse.Dists$percentage.diff)

plot(cow.Dists$S1.bin,cow.Dists$percentage.diff)

#cow.Dists$percentage.diff






