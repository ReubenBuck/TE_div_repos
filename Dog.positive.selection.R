





#get positive selection regions in dog
library(ggplot2)
require(gdata)
require(RCurl)
require(ggbio)
require(GenomicRanges)
url <- "http://www.plosgenetics.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pgen.1002316.s013"
pos.sel <- read.xls(url)
pos.sel$chr <- paste("chr", pos.sel$chr, sep = "")



spec1 <- "DOG"
spec2 <- "HUMAN"

setwd("~/Desktop/Repeat evolution/bin_synteny")

Dists <- read.table(file = paste("Div_score_pc1_pc2/", spec1, "_to_", spec2, "_Mdim_dist_scaled", sep = ""), header = TRUE)


s1 <- read.table(paste("~/Documents/phd/Desktop analyses/new_PCA/sort results/",spec1,"_Scaled", sep = ""), header = TRUE)
s1.gr <- GRanges( seqnames = Rle(s1$chr), 
					ranges = IRanges(start = s1$start, end = s1$end),
					bin_ID = s1$binID)
s2 <- read.table(paste("~/Documents/phd/Desktop analyses/new_PCA/sort results/",spec2,"_Scaled", sep = ""), header = TRUE)





# get PCA coordinate to do some sort of MDS
pca <- prcomp(s1[,5:length(s1)], scale.=FALSE)
biplot(pca , xlabs = rep(".", dim(pca$x)[1]))
plot(pca$x, pch = 16)
PC <- data.frame(S1.bin = s1[,1], S1PC1 = pca$x[,1], S1PC2 = pca$x[,2])
merg <- merge(Dists, PC)
Interesting <- pca$rotation[,1]

qplot(merg$S1PC1, merg$S1PC2, colour = merg$S2_to_S1_ratio) + scale_colour_gradient(limits=c(min(merg$S2_to_S1_ratio),max(merg$S2_to_S1_ratio)), high = "blue", low = "yellow")
qplot(merg$S1PC1, merg$S1PC2, colour = merg$percentage.diff) + scale_colour_gradient(limits=c(min(merg$percentage.diff),max(merg$percentage.diff)), high = "blue", low = "yellow")


# get the count values of each TE
s1a <- s1
colnames(s1a)[1] <- "S1.bin"
merg2 <- merge(Dists,s1a)
plot(merg2$percentage.diff, merg2$CR1__Eutheria)


# plot things on chromosomes

# need to get seq lengths
UCSCspec1 <- "hg19"
con <- gzcon(url(paste("http://hgdownload.soe.ucsc.edu/goldenPath/",UCSCspec1,"/database/chromInfo.txt.gz", sep="")))
txt <- readLines(con)
dat <- read.table(textConnection(txt))
slen <- dat[,2]
names(slen) <- dat[,1]
A <- as.data.frame(seqinfo(s1.gr))
slen <- slen[names(slen) %in% rownames(A)]
seqlengths(s1.gr) <- slen[order(names(slen))]
q <- autoplot(sort(s1.gr),layout = "karyogram",fill = "orange")
q <- q + layout_karyogram(sort(reduce(bt)), cgeom = "rect", ylim = c(11,21), color = "blue")
print(q)

# colour in areas of various TE dists
# see if there is anything noticable 
quant <- quantile(merg2$percentage.diff)
      
Q1 <- GRanges(seqnames= Rle(merg2$chr[merg2$percentage.diff > quant[1] & merg2$percentage.diff < quant[2]]), ranges = IRanges(start = merg2$start[merg2$percentage.diff > quant[1] & merg2$percentage.diff < quant[2]], end = merg2$end[merg2$percentage.diff > quant[1] & merg2$percentage.diff < quant[2]]))

Q2 <- GRanges(seqnames= Rle(merg2$chr[merg2$percentage.diff > quant[2] & merg2$percentage.diff < quant[3]]), ranges = IRanges(start = merg2$start[merg2$percentage.diff > quant[2] & merg2$percentage.diff < quant[3]], end = merg2$end[merg2$percentage.diff > quant[2] & merg2$percentage.diff < quant[3]]))

Q3 <- GRanges(seqnames= Rle(merg2$chr[merg2$percentage.diff > quant[3] & merg2$percentage.diff < quant[4]]), ranges = IRanges(start = merg2$start[merg2$percentage.diff > quant[3] & merg2$percentage.diff < quant[4]], end = merg2$end[merg2$percentage.diff > quant[3] & merg2$percentage.diff < quant[4]]))

Q4 <- GRanges(seqnames= Rle(merg2$chr[merg2$percentage.diff > quant[4] & merg2$percentage.diff < quant[5]]), ranges = IRanges(start = merg2$start[merg2$percentage.diff > quant[4] & merg2$percentage.diff < quant[5]], end = merg2$end[merg2$percentage.diff > quant[4] & merg2$percentage.diff < quant[5]]))

pos.sel.gr <- GRanges(seqnames = Rle(pos.sel$chr), ranges = IRanges(start = pos.sel$start, end = pos.sel$end))

q <- autoplot(sort(s1.gr),layout = "karyogram" ,fill = NA, col = NA)
q <- q + layout_karyogram(sort(reduce(Q1)), cgeom = "rect", ylim = c(1,10), fill = "blue")
q <- q + layout_karyogram(sort(reduce(Q2)), cgeom = "rect", ylim = c(1,10), fill = "purple")
q <- q + layout_karyogram(sort(reduce(Q3)), cgeom = "rect", ylim = c(1,10), fill = "red")
q <- q + layout_karyogram(sort(reduce(Q4)), cgeom = "rect", ylim = c(1,10), fill = "orange")
q <- q + layout_karyogram(sort(reduce(pos.sel.gr)), cgeom = "rect", ylim = c(11,21), fill = "pink")
print(q)


# So a PCA without unplaced chromosomes as they seem to have a very high density for dog specific repeats
s1new <- read.table(paste("~/Documents/phd/Desktop analyses/new_PCA/sort results/",spec1,"_COUNT", sep = ""), header = TRUE)
merg4 <- merg2[-(grep("Un", merg2[,"chr"])),]

# do verious plots and look for corelations
qplot(merg4$S1PC1, merg4$S1PC2, colour = merg4$percentage.diff) + scale_colour_gradient(limits=c(min(merg4$percentage.diff),max(merg4$percentage.diff)), high = "blue", low = "yellow")
pp <- prcomp(merg4[,11:length(merg4)])
biplot(pp)
merg4$S1PC1 <- pp$x[,1]
merg4$S1PC2 <- pp$x[,2]
C <- cor(merg4[,c(4,11:length(merg4))])
qplot(merg4$S1PC1, merg4$S1PC2, colour = merg4$S2_to_S1_ratio) + scale_colour_gradient(limits=c(min(merg4$S2_to_S1_ratio),max(merg4$S2_to_S1_ratio)), high = "blue", low = "yellow")

plot(merg4$S2_to_S1_ratio,merg4$SINE2.tRNA__Canidae )

# Found regions of positive selection
# lets just colour in bins that overlap with positive selection 

# we could potentialy do a regression analysis to see if there is some sort of relationship with areas identified with being positivly selected and the TE divergence

pos.sel.fdr5 <- pos.sel[pos.sel[,"breed"] == "EBD",]

pos.sel.fdr5.gr <- GRanges(seqnames = Rle(pos.sel.fdr5$chr), ranges = IRanges(start = pos.sel.fdr5$start, end = pos.sel.fdr5$end))
pos.sel.gr <- GRanges(seqnames = Rle(pos.sel$chr), ranges = IRanges(start = pos.sel$start, end = pos.sel$end))

merg4.gr <- GRanges(seqnames = Rle(merg4$chr), ranges = IRanges(start = merg4$start, end = merg4$end))

OV.pos <- as.matrix(findOverlaps(merg4.gr,pos.sel.fdr5.gr, minoverlap = 500000))


merg4$pos <- 1
merg4[unique(OV.pos[,1]),"pos"] <- 2

plot(merg4$S1PC1,merg4$percentage.diff , col = merg4$pos, pch = 16)
plot(merg4$S1PC1[merg4$pos == 2],merg4$percentage.diff[merg4$pos == 2] , pch = 16)












