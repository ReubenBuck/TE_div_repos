# This will be our clean method for doing everything. I will try to use the smallest amount of variables as possible and use unambigous names for them

# plenty of comments

# the ability to upload any species genome
# score name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks
# awk 'int($4)>=200' cow2humanAlignment, thsi command is used in shell to get alignment pieces 200 bp long

# write down species here
spec1 <- "COW"
spec2 <- "HUMAN"

rm( list = ls())

setwd("~/Desktop/Repeat evolution/bin_synteny")
## write it in R to underastand the problem
a <- read.table("./usable_alignment/chr_align_C_H")
a[,1] <- rownames(a)
s1 <- read.table("./PCA/bins_bt6_pca_order.txt")
s2 <- read.table("./PCA/bins_hg19_pca_regions.txt")


require(GenomicRanges)
require(BSgenome.Hsapiens.UCSC.hg19)

# the conditional below is for the two different format chromosome information comes in for diffferent species

# if using LAST alignments there's a good chance I won't need this
s2.genome <- BSgenome.Hsapiens.UCSC.hg19

SEQ<- seqlengths(s2.genome)
if(length(grep("_", seqnames(s2.genome))) == 0){
	SEQ2 <- width(s2.genome$chrUn)
	names(SEQ2) <- names(s2.genome$chrUn)
	SEQ <- c(SEQ,SEQ2)
}
SEQ <- data.frame(names(SEQ), SEQ)
neg.ali <- a[a[,8] == "-",]
neg.ali <- merge(neg.ali, SEQ, by.x = 5, by.y = 1)
neg.ali[,7] <- neg.ali[,10] - neg.ali[,7]
neg.ali[,6] <- neg.ali[,10] - neg.ali[,6]
# remeber to reverse start and end so there are no negative widths
neg.ali <- data.frame(neg.ali$V1, neg.ali$V2, neg.ali$V3, neg.ali$V4, neg.ali$V5, neg.ali$V7, neg.ali$V6, neg.ali$V8, neg.ali$V9)
colnames(neg.ali) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
a <- a[a[,8] == "+",]
a <- rbind(a,neg.ali)
rm(neg.ali)
rm(SEQ)


#### remove one to many mappers from reference and query
a.S2.gr <- GRanges( seqnames = Rle(a[,5]),
					ranges = IRanges(start = a[,6], end = a[,7],names = a[,1]),
					alignmentID = a[,1])

a.S1.gr <- GRanges( seqnames = Rle(a[,2]),
					ranges = IRanges(start = a[,3], end = a[,4],names = a[,1]),
					alignmentID = a[,1])

# perform self overlaps to identify if regions overlap themselves
multiS2 <- as.matrix(findOverlaps(a.S2.gr, a.S2.gr))
multiS1 <- as.matrix(findOverlaps(a.S1.gr, a.S1.gr))
# basicly species 1 does not ovelap itself
# however species 2 does
# therefore regions in species 1 are mapping to the same place in species two
# non multis
n.multi <- multiS2[!(multiS2[,1] %in% unique(multiS2[duplicated(multiS2[,1]),1])),1]
a <- a[n.multi,]



##########
# Back to sorting
# data should now have no repeats
##########


# species 1 alignments
a.S1.gr <- GRanges( seqnames = Rle(a[,2]),
					ranges = IRanges(start = a[,3], end = a[,4], names = a[,1]),
					alignmentID = a[,1])
# species 2 alignments
a.S2.gr <- GRanges( seqnames = Rle(a[,5]),
					ranges = IRanges(start = a[,6], end = a[,7],names = a[,1]),
					alignmentID = a[,1])

s1.gr <- GRanges( seqnames = Rle(s1[,1]), 
					ranges = IRanges(start = s1[,2], end = s1[,3]),
					bin_ID = s1[,5])
					
s2.gr <- GRanges( seqnames = Rle(s2[,1]), 
					ranges = IRanges(start = s2[,2], end = s2[,3]),
					bin_ID = s2[,5])



# identify alignment segments that overlap S1 bins and S2 bins					
OL.s1 <- as.matrix(findOverlaps(a.S1.gr, s1.gr))
OL.s2 <- as.matrix(findOverlaps(a.S2.gr, s2.gr))

# merge so that it is possible to identify the alignment section between eahc s1 and s2 bin
colnames(OL.s1) = c("a","species1")
colnames(OL.s2) = c("a","species2")
Merge_OL.s1.s2 <- merge(OL.s1,OL.s2)

Merged.DF <- cbind(as.data.frame(a.S1.gr[Merge_OL.s1.s2$a], row.names = NULL, stringsAsFactors = FALSE)[,c(1:4,6)],
			 as.data.frame(a.S2.gr[Merge_OL.s1.s2$a], row.names = NULL, stringsAsFactors = FALSE)[,c(1:4,6)],
			 as.data.frame(s1.gr[Merge_OL.s1.s2$species1], row.names = NULL, stringsAsFactors = FALSE)[,c(1:4,6)],
			 as.data.frame(s2.gr[Merge_OL.s1.s2$species2], row.names = NULL, stringsAsFactors = FALSE)[,c(1:4,6)])
			 
colnames(Merged.DF) <- c(paste("aS1", colnames(Merged.DF)[1:5], sep = "_"),
			paste("aS2", colnames(Merged.DF)[6:10], sep = "_"),
			paste("binS1", colnames(Merged.DF)[11:15], sep = "_"),
			paste("binS2", colnames(Merged.DF)[16:20], sep = "_"))
			
for(i in seq(length(Merged.DF))){
	if(length(grep("seqname", colnames(Merged.DF)[i])) == 1){Merged.DF[,i] <- as.character(Merged.DF[,i])}
	else if(length(grep("alignmentID", colnames(Merged.DF)[i])) == 1){Merged.DF[,i] <- as.character(Merged.DF[,i])}
	else{Merged.DF[,i] <- as.numeric(Merged.DF[,i])}	
}			



# look for duplications as to locate alignment pieces on the edge of two bins
# some are duplicated more than once
De <- Merged.DF[Merged.DF$aS2_alignmentID   %in%  (Merged.DF$aS2_alignmentID[duplicated(Merged.DF$aS2_alignmentID)]),]

# identiying duplications creates a dataframe in which part of an alignment belongs in one bin and part belongs in the other bin
# The following loop will go through and change values so they sit inside bins
# ie if the end of an alignment is > then the end of the bin its aligned to then the end of the alignment is cut
# likiwise there will be a piece where that start of the same alignment is < the start of a bin, this will also be cut.
# this will be done for S1 and S2
# may not actullay need to be run through a for loop as im pretty sure each element does not interact with any other element
# therefore i could grab all of one type of element and make the changes all at once

# things on the boundries in both species form a negative byproduct which is discarded

i = NULL
# split border pieces accordingly using proportions in order to make it even between both species
for(i in seq(dim(De)[1])){
	
	#pieces where the start aS1 is less than the start of binS1
	# determine by what proportion this occurs and then remove it from the aignment piece in each species 
	if(De$aS1_start[i] < De$binS1_start[i]) {
		x <- (De$binS1_start[i] - De$aS1_start[i])/De$aS1_width[i]
		De$aS1_start[i] <- De$aS1_start[i] + (De$aS1_width[i]*x)
		De$aS1_width[i] <- De$aS1_width[i] - (De$aS1_width[i]*x)
		De$aS2_start[i] <- De$aS2_start[i] + (De$aS2_width[i]*x)
		De$aS2_width[i] <- De$aS2_width[i] - (De$aS2_width[i]*x)
    } 
 
	if(De$aS1_end[i] > De$binS1_end[i]){
		x <- (De$aS1_end[i] - De$binS1_end[i])/De$aS1_width[i]
 		De$aS1_end[i] <- De$aS1_end[i] - (De$aS1_width[i]*x)
		De$aS1_width[i] <- De$aS1_width[i] - (De$aS1_width[i]*x)
		De$aS2_end[i] <- De$aS2_end[i] - (De$aS2_width[i]*x)
		De$aS2_width[i] <- De$aS2_width[i] - (De$aS2_width[i]*x)
	} 

	if(De$aS2_start[i] < De$binS2_start[i]){
		x <- (De$binS2_start[i] - De$aS2_start[i])/De$aS2_width[i]
		De$aS2_start[i] <- De$aS2_start[i] + (De$aS2_width[i]*x)
		De$aS2_width[i] <- De$aS2_width[i] - (De$aS2_width[i]*x)
		De$aS1_start[i] <- De$aS1_start[i] + (De$aS1_width[i]*x)
		De$aS1_width[i] <- De$aS1_width[i] - (De$aS1_width[i]*x)
	} 
	
	if(De$aS2_end[i] > De$binS2_end[i]){ 
		x <- (De$aS2_end[i] - De$binS2_end[i])/De$aS2_width[i]
		De$aS2_end[i] <- De$aS2_end[i] - (De$aS2_width[i]*x)
		De$aS2_width[i] <- De$aS2_width[i] - (De$aS2_width[i]*x)
		De$aS1_end[i] <- De$aS1_end[i] - (De$aS1_width[i]*x)
		De$aS1_width[i] <- De$aS1_width[i] - (De$aS1_width[i]*x)
	} 
}
De <- De[!(De$aS1_width < 0),]
Merged.DF[as.numeric(rownames(De)),] <- De






# now work on seperating into larger groups

# work out the actual proportion of shared bins between two species
M.g <- Merged.DF
# turn widths into proportions 
M.g[,c(4,9)] <- M.g[,c(4,9)]/1500000

# M.a is supposed to be widths inside each bin

M.a <- M.g[,c(15,4,20,9)]
M.a <- as.matrix(M.a)

Hbins <- unique(M.a[,1])

s.bin <- NULL
for (i in seq(along=Hbins)){
	b <- M.a[(Hbins[i] == M.a[,1]),]
	if(length(b) == 4){b <- t(as.matrix(b)) }
	b1 <- unique(b[,3])
	c2 = NULL
	# need to add up all the portions of all the bins that align with S1
	for (o in seq(along=b1)){
		c <- b[(b1[o] == b[,3]),]
		if (length(c) == 4){
			c <- t(data.frame(c))}
		else if (length(c) != 4){
			c <- data.frame(c)}
		c1 <- data.frame(c[1,1], sum(c[,2]), c[1,3], sum(c[,4]))
		c2 <- rbind(c2,c1)
	}
	s.bin <- rbind(s.bin,c2)
}



colnames(s.bin) <- c("S1.bin", "S1.Proportion", "S2.bin", "S2.Proportion")

write.table(s.bin, file=paste("./S_bins/", spec1, "aligning", spec2, sep = ""), sep = "\t", quote = FALSE,row.names = FALSE)


# there is still errors in the overlaps becasue of absent ch names

