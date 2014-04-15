# This will be our clean method for doing everything. I will try to use the smallest amount of variables as possible and use unambigous names for them

# plenty of comments

# the ability to upload any species genome
# score name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks
# awk 'int($4)>=200' cow2humanAlignment, thsi command is used in shell to get alignment pieces 200 bp long




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



					
L.h <- as.matrix(findOverlaps(a.S1.gr, S1.gr))
L.c <- as.matrix(findOverlaps(a.S2.gr, S2.gr))

# lets try the merge function
colnames(L.h) = c("a","species1")
colnames(L.c) = c("a","species2")
M <- merge(L.h,L.c)
M.g <- cbind(as.data.frame(a.S1.gr[M[,1]], row.names = NULL),as.data.frame(a.S2.gr[M[,1]], row.names = NULL),as.data.frame(S1.gr[M[,2]], row.names = NULL),as.data.frame(S2.gr[M[,3]], row.names = NULL))
M.g <- M.g[,c(1,2,3,4,6,7,8,9,10,12,13,14,15,16,18,19,20,21,22,24)]
De <- M.g[M.g[,10] %in% M.g[duplicated(M.g[,10]),10],]


i = NULL
# human border pieces
for(i in seq(dim(De)[1])) 
     {
 if(De[i,2] < De[i,12]) {        
 x <- (De[i,12] - De[i,2])/De[i,4]
 De[i,2] <- as.integer(De[i,2] + (De[i,4]*x))
 De[i,4] <- as.integer(De[i,4] - (De[i,4]*x))
 De[i,7] <- as.integer(De[i,7] + (De[i,9]*x))
 De[i,9] <- as.integer(De[i,9] - (De[i,9]*x))
    } 
 
 if(De[i,3] > De[i,13]) 
 { 
 x <- (De[i,3] - De[i,13])/De[i,4]
 De[i,3] <- as.integer(De[i,3] - (De[i,4]*x))
 De[i,4] <- as.integer(De[i,4] - (De[i,4]*x))
 De[i,8] <- as.integer(De[i,8] - (De[i,9]*x))
 De[i,9] <- as.integer(De[i,9] - (De[i,9]*x))
  } 

if(De[i,7] < De[i,17]) 
 { 
 x <- (De[i,17] - De[i,7])/De[i,9]
 De[i,7] <- as.integer(De[i,7] + (De[i,9]*x))
 De[i,9] <- as.integer(De[i,9] - (De[i,9]*x))
 De[i,2] <- as.integer(De[i,2] + (De[i,4]*x))
 De[i,4] <- as.integer(De[i,4] - (De[i,4]*x))
 } 

 if(De[i,8] > De[i,18]) 
 { 
 x <- (De[i,8] - De[i,18])/De[i,9]
 De[i,8] <- as.integer(De[i,8] - (De[i,9]*x))
 De[i,9] <- as.integer(De[i,9] - (De[i,9]*x))
 De[i,3] <- as.integer(De[i,3] - (De[i,4]*x))
 De[i,4] <- as.integer(De[i,4] - (De[i,4]*x))
 } 
 }

M.g[as.numeric(rownames(De)),] <- De
M.g <- (M.g[!(M.g[,9] < 0), ])

# loop works well at spliting boarders 
# now work on seperating into larger groups

M.g[,c(4,9)] <- M.g[,c(4,9)]/1500000
M.a <- M.g[,c(15,4,20,9)]
M.a <- as.matrix(M.a)

Hbins <- unique(M.a[,1])

b2 <- NULL
for (i in seq(along=Hbins)){
	b <- M.a[(Hbins[i] == M.a[,1]),]
	if(length(b) == 4){b <- t(as.matrix(b)) }
	b1 <- unique(b[,3])
	c2 = NULL
	for (o in seq(along=b1)){
		c <- b[(b1[o] == b[,3]),]
		if (length(c) == 4){
			c <- t(data.frame(c))}
		else if (length(c) != 4){
			c <- data.frame(c)}
		c1 <- data.frame(c[1,1], sum(c[,2]), c[1,3], sum(c[,4]))
		c2 <- rbind(c2,c1)
	}
	b2 <- rbind(b2,c2)
}

s.bin <- b2

colnames(s.bin) <- c("H.bin", "H.P", "C.bin", "C.P")

write.table(s.bin, file="./S_bins/C_bin_aligning_H_bin_independant1", sep = "\t", quote = FALSE,row.names = FALSE)




