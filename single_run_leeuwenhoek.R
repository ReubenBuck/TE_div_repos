# Bin sort and MD analysis using data processed by reuben
#
#



# if it is possible get the neigborhood identification faster

# also install the table generation.


# then this whole thing can be run on leeuwenhoek and we cna use multiple processors and do multiple species at the same time



# sort bins up here and remove chromosomes later on
# Ican also run this on joys stuff but I'm going to have to run various scaling options
# and options on what variables I choose to keep in
# Then I can change the level at which we view genome maintenece 




# jobs to do on leuwenhoek
# ensure colnames are correct on joys files
# ensure we are actually able to pick and choose the variables we want
# make it so the file names all conform to each other
# download and process the appropriate alignmnet files on UCSC
# and also the chromosome information




rm( list = ls())


require("doParallel")
require(GenomicRanges)

#parallel options
cl <- makeCluster(8)
registerDoParallel(cl)


# option to remove unplaced chromosomes
rem.un <- "yes"

# option to decide the various neighborhood sizes to be used 
# write out a vector containg the sizes to analyse

sizes <- c(5,20,50,100,150,200,300,500,750)

# for spec names first letter upper, the rest lower
# also use common names except for cow which is Bovine

spec1 <- "Mouse"
spec2 <- "Human"
UCSCspec2 <- "hg19"

setwd("/scratch/reuben_counting/TE_div/")
## write it in R to underastand the problem
a <- read.table(paste("usable_alignemnts/chr_align_", spec1, "_", spec2, sep = ""))
a[,1] <- rownames(a)

# these files have bin information combined with dimensions
# conditional to make sure I'm using one of the two names used for this data


specfile <- list.files(paste("/scratch/ForReuben201404/",spec1, sep = ""))
if(any(specfile == "AllBinCounts.txt")){
		s1name <- paste("/scratch/ForReuben201404/",spec1,"/AllBinCounts.txt", sep = "")
		}else{
		s1name <- paste("/scratch/ForReuben201404/",spec1,"/binData.txt", sep = "")}

specfile <- list.files(paste("/scratch/ForReuben201404/",spec2, sep = ""))
if(any(specfile == "AllBinCounts.txt")){
                s2name <- paste("/scratch/ForReuben201404/",spec2,"/AllBinCounts.txt", sep = "")
                }else{
		s2name <- paste("/scratch/ForReuben201404/",spec2,"/binData.txt", sep = "")}

s1 <- read.table(s1name, header = TRUE)
if(rem.un == "yes"){
	if(length(grep("U", s1$chr)) > 0){s1 <- s1[-(grep( "U", s1$chr)),]}
	if(length(grep("_", s1$chr)) > 0){s1 <- s1[-(grep("_", s1$chr)),]}
	if(length(grep("M", s1$chr)) > 0){s1 <- s1[-(grep("M", s1$chr)),]}
	}
s2 <- read.table(s2name, header = TRUE)
if(rem.un == "yes"){
	if(length(grep("U", s2$chr)) > 0){s2 <- s2[-(grep( "U", s2$chr)),]}
	if(length(grep("_", s2$chr)) > 0){s2 <- s2[-(grep("_", s2$chr)),]}
	if(length(grep("M", s2$chr)) > 0){s2 <- s2[-(grep("M", s2$chr)),]}
	}

# s1 and s2 require editing to remove unnecsary features for our purposes
slist <- list(s1,s2)
mb <- 1000000
bin.size = 500000

#which columns to keep 
keep.NGenes = "no"
keep.NG4s = "no"
keep.NCpGI = "no"
keep.CpGBP = "no"
keep.GC = "no"

for(i in seq(along=slist)){
      count <- slist[[i]]
      count <- count[count$Known >= bin.size,]
      count[,5:(length(count)-1)] <- count[,5:(length(count)-1)]/mb   
      if(keep.NGenes == "no"){count <- count[,!(colnames(count) == "NGenes")]}
      if(keep.NG4s ==	"no"){count <- count[,!(colnames(count) == "NG4s")]}
      if(keep.NCpGI ==	"no"){count <- count[,!(colnames(count) == "NCpGI")]}
      if(keep.CpGBP  ==	"no"){count <- count[,!(colnames(count) == "CpGBP")]}
      if(keep.GC ==	"no"){count <- count[,!(colnames(count) == "GC")]}
      count <- count[,!(colnames(count) == "Known")]
      colnames(count)[1:4] <- c("chr", "binID", "start", "end")
      count$binID <- 1:dim(count)[1]
      slist[[i]] <- count
}

s1 <- slist[[1]]
s2 <- slist[[2]]


# handling negative alignments
con <- gzcon(url(paste("http://hgdownload.soe.ucsc.edu/goldenPath/",UCSCspec2,"/database/chromInfo.txt.gz", sep="")))
txt <- readLines(con)
dat <- read.table(textConnection(txt))
SEQ <- dat[,1:2]
colnames(SEQ) <- c("chr", "size")
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
# however species 2 does
# therefore regions in species 1 are mapping to the same place in species two
# non multis
n.multi <- multiS2[!(multiS2[,1] %in% unique(multiS2[duplicated(multiS2[,1]),1])),1]
a <- a[n.multi,]
a[,1] <- rownames(a) <- 1:dim(a)[1]

#####
#
# Take out the repeat sequences

# "A" is a back up variable
A <- a

# Readin Repeat files


S1rep <- read.table(paste("Rep_files/all_chr_rep_", spec1, sep=""))
S2rep <- read.table(paste("Rep_files/all_chr_rep_", spec2, sep=""))

S1rep.gr <- reduce(GRanges(seqnames = Rle(S1rep[,1]), ranges = IRanges(start = S1rep[,2], end = S1rep[,3])))
S2rep.gr <- reduce(GRanges(seqnames = Rle(S2rep[,1]), ranges = IRanges(start = S2rep[,2], end = S2rep[,3])))


# reestablish alignment files after duplicate have been taken out
a.S2.gr <- GRanges( seqnames = Rle(a[,5]),
					ranges = IRanges(start = a[,6], end = a[,7],names = a[,1]),
					alignmentID = a[,1])

a.S1.gr <- GRanges( seqnames = Rle(a[,2]),
					ranges = IRanges(start = a[,3], end = a[,4],names = a[,1]),
					alignmentID = a[,1])


## seems like some of the cow things are multimapping
#  they do have scores so anything that overlaps the wrong way can be removed

##################
#
#     PART 1: removing repeats from humam
#
##################

# find the area from S1 that does not contain repeats
set.gr <-  setdiff(a.S1.gr,S1rep.gr)

# identify the nonrepeat sections that have been joined due to setdiff function 
OV <- as.matrix(findOverlaps(set.gr,a.S1.gr))
set.df <- as.data.frame(set.gr)
ali <- data.frame(OV[,1],set.df[OV[,1], 1:3], a[OV[,2], 1:4])

# the ali dataframe is the correct size just the numbers are wrong, all it needs is correct splitting and it will be fine
n.dup <- ali[!(ali[,1] %in% unique(OV[duplicated(OV[,1]),1])),]
# this covers all of the setdiffds that merged two or more a
dup.a <- ali[ali[,1] %in% unique(OV[duplicated(OV[,1]),1]),]
colnames(dup.a) <- c("p.id", "chr", "p.s", "p.e", "a.id", "chr", "a.s","a.e")


#start
# as <= ps & ae < pe
# with start we change the end
start <- dup.a[(dup.a$a.s - dup.a$p.s <= 0) & (dup.a$a.e - dup.a$p.e < 0),]
start[,4] <- start[,8]

# mid 
# as > ps & ae < pe
# with mid we change the end and start
mid <- dup.a[(dup.a$a.s - dup.a$p.s > 0) & (dup.a$a.e - dup.a$p.e < 0),]
mid[,3] <- mid[,7]
mid[,4] <- mid[,8]

# end
#as > ps & ae >= pe
# with end we change the start
end <- dup.a[(dup.a$a.s - dup.a$p.s > 0) & (dup.a$a.e - dup.a$p.e >= 0),]
end[,3] <- end[,7]

dup.a2 <- rbind(start,mid,end)
colnames(n.dup) <- colnames(dup.a2)

Dali <- rbind(n.dup,dup.a2)
Dali <- Dali[,c(5,2,3,4)]
dali.gr <- GRanges(seqnames = Rle(Dali[,2]), ranges = IRanges(start = Dali[,3], end = Dali[,4]))

##################
#
#     PART 2: removing the missing part from species 2
#
##################

# this step is the next challange to get something to actually work
# the above step works but i will go back to making it faster
# it does also need to be reversed to cow and so will this section later on
# after that we have independance

# seems preety easy 

# CsN = (Cl/Hl * dHs) + CsO
# the same for end

# take human repeats from cow
# take cow repeats from cow


Q <- as.matrix(findOverlaps(dali.gr,a.S1.gr))
ali.2 <- data.frame(Dali[Q[,1],], a[Q[,2],])

colnames(ali.2) <- c("U.a.id", "U.S1.chr", "U.S1.p.s", "U.S1.p.e", "I.a.id", "I.S1.chr", "I.S1.start", "I.S1.end", "I.S2.chr", "I.S2.start", "I.S2.end")


ali.2$S2.length <- ali.2$I.S2.end - ali.2$I.S2.start
ali.2$S1.length <- ali.2$I.S1.end - ali.2$I.S1.start
# remove those with width 0 as they no longer exist as alignment pieces
ali.2 <- ali.2[!(ali.2$S2.length < 1 | ali.2$S1.length < 1 ),]
ali.2$N.S2.chr <- ali.2$I.S2.chr
ali.2$N.S2.start <- as.integer(((ali.2$S2.length/ali.2$S1.length)*(ali.2$U.S1.p.s - ali.2$I.S1.start)) + ali.2$I.S2.start)
ali.2$N.S2.end <- as.integer(((ali.2$S2.length/ali.2$S1.length)*(ali.2$U.S1.p.e - ali.2$I.S1.start)) + ali.2$I.S2.start)

a <- ali.2[,c("U.a.id", "U.S1.chr" , "U.S1.p.s", "U.S1.p.e", "N.S2.chr", "N.S2.start", "N.S2.end")]
a[,1] <- rownames(a) <- 1:dim(a)[1]



##################
#
#     PART 3: removing repeats from cow
#
##################

a.S2.gr <- GRanges( seqnames = Rle(a[,5]),
					ranges = IRanges(start = a[,6], end = a[,7],names = a[,1]),
					alignmentID = a[,1])

a.S1.gr <- GRanges( seqnames = Rle(a[,2]),
					ranges = IRanges(start = a[,3], end = a[,4],names = a[,1]),
					alignmentID = a[,1])

set.gr <-  setdiff(a.S2.gr,S2rep.gr)

OV <- as.matrix(findOverlaps(set.gr,a.S2.gr))
set.df <- as.data.frame(set.gr)
ali <- data.frame(OV[,1],set.df[OV[,1], 1:3], a[OV[,2], c(1,5:7)])

# the ali dataframe is the correct size just the numbers are wrong, all it needs is correct splitting and it will be fine
n.dup <- ali[!(ali[,1] %in% unique(OV[duplicated(OV[,1]),1])),]
# this covers all of the setdiffds that merged two or more a
dup.a <- ali[ali[,1] %in% unique(OV[duplicated(OV[,1]),1]),]
colnames(dup.a) <- c("p.id", "chr", "p.s", "p.e", "a.id", "chr", "a.s","a.e")

#start
# as <= ps & ae < pe
# with start we change the end
start <- dup.a[(dup.a$a.s - dup.a$p.s <= 0) & (dup.a$a.e - dup.a$p.e < 0),]
start[,4] <- start[,8]

# mid 
# as > ps & ae < pe
# with mid we change the end and start
mid <- dup.a[(dup.a$a.s - dup.a$p.s > 0) & (dup.a$a.e - dup.a$p.e < 0),]
mid[,3] <- mid[,7]
mid[,4] <- mid[,8]

# end
#as > ps & ae >= pe
# with end we change the start
end <- dup.a[(dup.a$a.s - dup.a$p.s > 0) & (dup.a$a.e - dup.a$p.e >= 0),]
end[,3] <- end[,7]

dup.a2 <- rbind(start,mid,end)
colnames(n.dup) <- colnames(dup.a2)

Dali <- rbind(n.dup,dup.a2)
Dali <- Dali[,c(5,2,3,4)]
dali.gr <- GRanges(seqnames = Rle(Dali[,2]), ranges = IRanges(start = Dali[,3], end = Dali[,4]))

##################
#
#     PART 4: removing the missing part from human
#
##################



Q <- as.matrix(findOverlaps(dali.gr,a.S2.gr))

# there may be a problem with the missing chromaosmes in cow
# reverse columns so that cow comes first
ali.2 <- data.frame(Dali[Q[,1],], a[Q[,2],])


colnames(ali.2) <- c("U.a.id", "U.S2.chr", "U.S2.p.s", "U.S2.p.e", "I.a.id", "I.S1.chr", "I.S1.start", "I.S1.end", "I.S2.chr", "I.S2.start", "I.S2.end")

ali.2$S2.length <- ali.2$I.S2.end - ali.2$I.S2.start
ali.2$S1.length <- ali.2$I.S1.end - ali.2$I.S1.start
ali.2 <- ali.2[!(ali.2$S2.length < 1 | ali.2$S1.length < 1 ),]
ali.2$N.S1.chr <- ali.2$I.S1.chr
ali.2$N.S1.start <- as.integer(((ali.2$S1.length/ali.2$S2.length)*(ali.2$U.S2.p.s - ali.2$I.S2.start)) + ali.2$I.S1.start)
ali.2$N.S1.end <- as.integer(((ali.2$S1.length/ali.2$S2.length)*(ali.2$U.S2.p.e - ali.2$I.S2.start)) + ali.2$I.S1.start)
# it seems now that i have new alignments that no longer contain human repeats 


a <- ali.2[,c("U.a.id", "N.S1.chr" , "N.S1.start", "N.S1.end", "U.S2.chr", "U.S2.p.s", "U.S2.p.e")]
a[,1] <- rownames(a) <- 1:dim(a)[1]

if(rem.un == "yes"){
	if(length(grep("U", a$N.S1.chr)) > 0){a <- a[-(grep( "U", a$N.S1.chr)),]}
	if(length(grep("_", a$N.S1.chr)) > 0){a <- a[-(grep("_", a$N.S1.chr)),]}
	if(length(grep("M", a$N.S1.chr)) > 0){a <- a[-(grep("M", a$N.S1.chr)),]}
	if(length(grep("U", a$U.S2.chr)) > 0){a <- a[-(grep( "U", a$U.S2.chr)),]}
	if(length(grep("_", a$U.S2.chr)) > 0){a <- a[-(grep("_", a$U.S2.chr)),]}
	if(length(grep("M", a$U.S2.chr)) > 0){a <- a[-(grep("M", a$U.S2.chr)),]}
	}


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

# check to see if repeats and alignments are overlapping
# program will quit if there is an error					
XS1 <- as.matrix(findOverlaps(a.S1.gr, S1rep.gr))
if(length(XS1)>0){print("repeats and alignments are overlapping")
		quit(save = "no")}
XS2 <- as.matrix(findOverlaps(a.S1.gr, S1rep.gr))
if(length(XS2)>0){print("repeats and alignments	are overlapping")
                quit(save = "no")}


s1.gr <- GRanges( seqnames = Rle(s1$chr), 
					ranges = IRanges(start = s1$start, end = s1$end),
					bin_ID = s1$binID)
					
s2.gr <- GRanges( seqnames = Rle(s2$chr), 
					ranges = IRanges(start = s2$start, end = s2$end),
					bin_ID = s2$bin)



# identify alignment segments that overlap S1 bins and S2 bins					
OL.s1 <- as.matrix(findOverlaps(a.S1.gr, s1.gr))
OL.s2 <- as.matrix(findOverlaps(a.S2.gr, s2.gr))

# merge so that it is possible to identify the alignment section between eahc s1 and s2 bin
colnames(OL.s1) = c("a","species1")
colnames(OL.s2) = c("a","species2")
Merge_OL.s1.s2 <- merge(OL.s1,OL.s2)

Merged.DF <- cbind(as.data.frame(a.S1.gr[Merge_OL.s1.s2$a], row.names = NULL)[,c(1:4,6)],
			 as.data.frame(a.S2.gr[Merge_OL.s1.s2$a], row.names = NULL)[,c(1:4,6)],
			 as.data.frame(s1.gr[Merge_OL.s1.s2$species1], row.names = NULL)[,c(1:4,6)],
			 as.data.frame(s2.gr[Merge_OL.s1.s2$species2], row.names = NULL)[,c(1:4,6)])
			 
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





# if I could get this to work on everything i have saved then it will be much easier to run through everything
# the probelem is that I have to worry about all the pca stuff too


# load everything
# ensure there are the correct species in each variable



#####
#This works because the bin names are supposadly in the same order as their PCA values
# with future data this may not be the case

# merge bins adn coordinates into a single data frame 
S1.bin.coord <- data.frame(s1$binID, s1[,5:length(s1)])
colnames(S1.bin.coord) <- c("S1.bin", paste("S1TE",colnames(s1[,5:length(s1)]), sep = "_"))
S2.bin.coord <- data.frame(s2$binID, s2[,5:length(s2)])
colnames(S2.bin.coord) <- c("S2.bin", paste("S2TE",colnames(s2[,5:length(s2)]), sep = "_"))
bin.coord <- merge(merge(s.bin,S1.bin.coord), S2.bin.coord)


###################
###
#
###    Begin to run parallel loop here that can handle the process at various sizes of "no"
#
###
###################







dist.all <- foreach(i = sizes) %dopar% { no <- i
	
no <- i
S1bin <- unique(bin.coord$S1.bin)
S1.N.dist <- NULL
S1ploters <- unique(bin.coord[,c(grep("S1.bin", colnames(bin.coord)), grep("S1TE", colnames(bin.coord)))])

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

# make this faste over multiple dimensions
# if possible find a way to vectorise it rater than loop it


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
			for( z in grep("S2TE",colnames(B_points))){
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
	row_final <- data.frame(S1.bin = S1bin[i],  S2.mean.dist , S1.mean.dist = S1.N.dist$ave.dist[S1.N.dist$S1.bin == S1bin[i]][1])
	Dists <- rbind(Dists, row_final)

	
	}
	

Dists$percentage.diff <- sqrt((((Dists$S1.mean.dist - Dists$S2.mean.dist)/ Dists$S1.mean.dist) * 100)^2)
Dists$diff <- sqrt((((Dists$S1.mean.dist - Dists$S2.mean.dist)))^2)
Dists$S1_to_S2_ratio <- Dists$S1.mean.dist / Dists$S2.mean.dist
Dists$S2_to_S1_ratio <- Dists$S2.mean.dist / Dists$S1.mean.dist
	
write.table(Dists, file = paste("Div_score/", spec1, "_to_", spec2, "_Mdim_dist_scaled_N.size_", no, sep = ""), quote = FALSE, sep = "\t", row.names= FALSE )
	
	
}
	
	



