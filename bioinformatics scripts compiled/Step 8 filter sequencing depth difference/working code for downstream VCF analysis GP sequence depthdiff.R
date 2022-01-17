#######prepping the data

??read.vcfR

library(vcfR)
setwd("D:/Dropbox/RA/Kat/GP bioinformatics/Filtering_steps/1 Basic filters for sequencing depth")
vcf <- read.vcfR("vcf_filtered_MAFmeanDP10minDP5.vcf")

#chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)

#lets fix up the naming scheme first

colnames(vcf@gt)

names <- read.csv("names.csv", header=FALSE)
names$V1
colnames(vcf@gt) <- names$V1

colnames(vcf@gt)
vcf@meta

gt <- extract.gt(vcf, IDtoRowNames = F)
fixed <- getFIX(vcf)
snps <- cbind(fixed[,1:5], gt)
head(snps)[,1:10]


snps.1 <- as.data.frame(as.matrix(snps))
snps.1$CHROM <- as.character(as.factor(snps.1$CHROM))
snps.1$POS <- as.character(as.factor(snps.1$POS))
snps.1$ID <- as.character(as.factor(snps.1$ID))
snps.1$identifier <- with(snps.1, paste0(CHROM, POS, ID))
snps.3 <- snps.1[,1:(ncol(snps.1)-1)]


##############

######## read depth

read.depth <- extract.gt(vcf, element="AD")
length(unique(rownames(read.depth)))
nrow(read.depth)
read.depth.ref <- masplit(read.depth, record = 1, sort=0)
read.depth.snp <- masplit(read.depth, record = 2, sort=0)
# this worked. Can generate estimates for combined avg read depth plus min read depth per ind! (do as combined. the relative read depth was 
# ridiculous)
# So the code calculates the read depth value off only sections which have at least 1 alt/ref allele count (ignores 0 and the MAF issue)
# So realistically the read count ref ans SNP are done. Just need to do combined\
# So my thinking for the read.depth.ref and read.depth.snp data frames is I want to filer out any genotype sequences that are less than 2 for either 
# but only if it is a heterozygote. I don't want to work in averages here or remove the actual loci. I want to case anything as less than 
# 2 as missing and recode to NA as needed for individuals per loci.

# So to do this I need a simple loop to read the numbers in either read.depth.ref or read.depth.snp but many say 0 so I need to have a secondary command that reads in zero 
# and ignores it. 

#it might just be easier to do it from combined and do a min 

read.depth.combined <- (read.depth.ref+read.depth.snp)

head(read.depth.ref)[1:6,1:6]
read.depth.ref.count<- rowSums(read.depth.ref, na.rm=T)
head(read.depth.ref.count)

read.depth.ref <- as.data.frame(read.depth.ref)
read.depth.ref$length <- rep(NA)

n <- ncol(read.depth.ref)-1
for (r in 1:nrow(read.depth.ref)) {
  read.depth.ref$length[r] <- length(which(read.depth.ref[r,1:n] !=0))
}

read.depth.ref.avg <- read.depth.ref.count/read.depth.ref$length
head(read.depth.ref.avg)
summary(read.depth.ref.avg)

read.depth.snp.count<- rowSums(read.depth.snp, na.rm=T)
head(read.depth.snp.count)

read.depth.snp <- as.data.frame(read.depth.snp)
read.depth.snp$length <- rep(NA)

n2 <- ncol(read.depth.snp)-1
for (r in 1:nrow(read.depth.snp)) {
  read.depth.snp$length[r] <- length(which(read.depth.snp[r,1:n2] !=0))
}

read.depth.snp.avg <- read.depth.snp.count/read.depth.snp$length


head(read.depth.combined)[1:6,1:6]
read.depth.combined.count<- rowSums(read.depth.combined, na.rm=T)
head(read.depth.combined.count)

read.depth.combined <- as.data.frame(read.depth.combined)
read.depth.combined$length <- rep(NA)

n <- ncol(read.depth.combined)-1
for (r in 1:nrow(read.depth.combined)) {
  read.depth.combined$length[r] <- length(which(read.depth.combined[r,1:n] !=0))
}

read.depth.combined.avg <- read.depth.combined.count/read.depth.combined$length
head(read.depth.combined.avg)

summary(read.depth.combined.avg)
summary(read.depth.snp.avg)
summary(read.depth.ref.avg)




#Figure out what you want the read depth average cut off to be
length(read.depth.ref.avg)
length(read.depth.snp.avg)
length(which(read.depth.ref.avg >= 1))
length(which(read.depth.snp.avg >= 1))
length(which(read.depth.ref.avg < 1))
length(which(read.depth.snp.avg < 1))
length(which(read.depth.combined.avg >= 15))

hist(read.depth.ref.avg)
hist(read.depth.snp.avg)

coverage.rd <- cbind(read.depth.ref.avg, read.depth.snp.avg)
coverage.rd <- as.data.frame(coverage.rd)
coverage.rd$snp.index <- 1:nrow(coverage.rd)
coverage.rd1 <- coverage.rd[which(coverage.rd$read.depth.ref.avg >= 1.0),]
nrow(coverage.rd)
nrow(coverage.rd1)
coverage.rd2 <- coverage.rd1[which(coverage.rd1$read.depth.snp.avg >= 1.0),]
nrow(coverage.rd2)

par(mfrow=c(2,2))
hist(read.depth.ref.avg, main="Read depth of ref allele", xlab="Read depth")
hist(read.depth.snp.avg, main="Read depth of snp allele", xlab="Read depth")
hist(coverage.rd2$read.depth.ref.avg, main="Read depth of ref allele > 2.0", xlab="Read depth")
hist(coverage.rd2$read.depth.snp.avg, main="Read depth of snp allele > 2.0", xlab="Read depth")

index <- 1:nrow(snps.3)
snps.index <- cbind(index, snps.3)
snps.rd <- snps.index[which(snps.index[,1] %in% coverage.rd2$snp.index),]
nrow(snps.rd)

# Filtered for a read depth minimum of 4
snps.rd

#Filter by coverage between the SNP and the ref allele to make sure there isn't calling bias
coverage.rd2$max <- pmax(coverage.rd2$read.depth.ref.avg, coverage.rd2$read.depth.snp.avg)
coverage.rd2$diff <- ((abs(coverage.rd2$read.depth.ref.avg - coverage.rd2$read.depth.snp.avg))/(coverage.rd2$max))*100
hist(coverage.rd2$diff, main="Coverage difference", xlab="% diff in coverage")
###haven't actually run a filter yet but all fell below threshold difference except 1
length(which(coverage.rd2$diff <= 80))
length(which(coverage.rd2$diff >= 80))

coverage.rd3 <- coverage.rd2[which(coverage.rd2$diff <= 80),]
coverage.rd3 <- coverage.rd3[which(coverage.rd3$diff >= 20),]

snps.rd1_2 <- snps.index[which(snps.index[,1] %in% coverage.rd3$snp.index),]
# we want to use the Dart based filter of 5X difference. this is all whch are above 80% 

snps.rd2 <- snps.rd1_2
#### Call rate
snps.rd[1:10,1:10]
callrate <- apply(snps.rd1_2, 1, function(x) 100-(sum(is.na(x))/(ncol(snps.rd1_2)-5))*100)
callrateind <- apply(snps.rd1_2,MARGIN=2, function(x) 100-(sum(is.na(x))/(nrow(snps.rd1_2)-6))*100)


sum(callrateind>90)

sum(callrate>0)

snps.rd2$callrate <- callrate

# this one is fine. 
snps.rd2 <- subset(snps.rd2, snps.rd2$callrate>=0)
collen <- ncol(snps.rd2)
snps.rd2 <- snps.rd2[,-collen]



callrateind
hist(callrateind)
callrateind <- callrateind[-(1:6)]

snps.rd3 <-rbind(snps.rd2,callrateind)

snps.rd3 <- t(snps.rd3)
snps.rd3 <- as.data.frame(snps.rd3)
collen <-ncol(snps.rd3)
snps.rd3[1:6,collen] <- 99

snps.rd3 <- subset(snps.rd3, snps.rd3[,collen] >= 10)

snps.rd3 <- snps.rd3[,-collen]

snps.rd3 <- t(snps.rd3)
snps.rd3 <- as.data.frame(snps.rd3)

snps.rd3

nrow(snps.rd3)
ncol(snps.rd3)
# Now that they are filtered by call rate and ind call rate
# We can do MAF


snps.rd3




##hist(callrate, main="Call Rate", xlab="Call Rate")
#```
#Number of SNPs retained with call rate threshold:
#  ```{r}
length(which(callrate >= 30))
###haven't actually run a filter yet but all fell below threshold difference so don't need to

#snps.rd$callrate <- callrate
callrate[1:10]

# we have 
coverage.rd2
callrate
callrateind

#snps.rd <- snps.index[which(snps.index[,1] %in% coverage.rd2$snp.index),]


#######MAF filter
hist(callrate)
snps.rd4 <- snps.rd3
snps.rd4 <- as.data.frame(snps.rd4)
snps.rd4$refcount <- rep(NA)
n3 <- ncol(snps.rd4)
for (r in 1:nrow(snps.rd4)) {
  snps.rd4$refcount[r] <- 2*(length(which(snps.rd4[r,7:n3] == "0/0"))) + 
    length(which(snps.rd4[r,7:n3] == "0/1"))
}

snps.rd4$altcount <- rep(NA)
for (r in 1:nrow(snps.rd4)) {
  snps.rd4$altcount[r] <- 2*(length(which(snps.rd4[r,7:n3] == "1/1"))) + 
    length(which(snps.rd4[r,7:n3] == "0/1"))
}


snps.rd4$minor <- pmin(snps.rd4$refcount, snps.rd4$altcount)
snps.rd4$total <- snps.rd4$refcount + snps.rd4$altcount
snps.rd4$maf <- snps.rd4$minor/snps.rd4$total
hist(snps.rd4$maf, main="Minor Allele Frequency", xlab="MAF")
#```
#Number of SNPs retained with MAF threshold:
#  ```{r}
length(which(snps.rd4$maf > 0.05))
#```








snps.rd5 <-subset(snps.rd4,snps.rd4$maf > 0.0000)
# find final 5 columns added in
collen<- ncol(snps.rd5)

snps.rd5[1:10,(collen-4):collen]
snps.rd5[1:10,559:571]
# remove final 5
snps.rd5 <- snps.rd5[,-c((collen-4):collen)]
snps.rd5[1:10,559:565]

snps.rd5


# So there are a few steps required to filter the vcf file back to good use. First I need to filter the vcf@fix component to remove 
# all loci already removed (this can be done with a dplyr filter command based of POS)

# That is the easy step. The hard step is the vcf@gt section. This needs to be filtered by loci and by individual to remove all out of them
# Except there is no identification for each locus 

# i need to index the gt markers. I also need an index to remain in the SNPs, and filter through that index. 

vcf2 <- vcf

# Figure out what the class of index is (note it's integer)
class(snps.rd5$index)
snps.rd4$index
snps.rd3$index
class(snps.rd2$index)


# generate an index column for future filtering
vcf2gt <- vcf2@gt
vcf2gt <- as.data.frame(vcf2gt)
vcf2gt$index <-coverage.rd$snp.index

# convert final snp.rd to correct class for integer column
snps.rd6 <- snps.rd5
snps.rd6$index <- as.integer(snps.rd6$index)

library(dplyr)

# Filter out index through dplyr function
library(dplyr)

vcf2gtfiltered <- vcf2gt%>%
  filter((index %in% snps.rd6$index))


# success, now repeat for individuals as inds were removed to 

# so after reading, it is not practical to filter by column in dpylyr... again. So transpose time

colnames(snps.rd6)
colnames(vcf2gtfiltered)
snps.rd6$FORMAT <- NA
colnames(snps.rd6)
snps.rd7 <- t(snps.rd6)
snps.rd7 <- as.data.frame(snps.rd7)
snps.rd7$inds <- rownames(snps.rd7)

vcf2gtfiltered2 <- t(vcf2gtfiltered)
vcf2gtfiltered2 <- as.data.frame(vcf2gtfiltered2)
vcf2gtfiltered2$inds <- rownames(vcf2gtfiltered2)


vcf2gtfiltered2 <- vcf2gtfiltered2%>%
  filter((inds %in% snps.rd7$inds))

vcf2gtfiltered3 <- t(vcf2gtfiltered2)
vcf2gtfiltered3 <- as.data.frame(vcf2gtfiltered3)


# This worked and inds have now been removed, Note still need to remove index column and ind row

# what column is index and what row is inds 
# note they are the final col and row in the dataframe

vcf2gtfiltered3[1:5,550:555]
vcf2gtfiltered3[1:5,555]

vcf2gtfiltered3[2764:2764,550:555]
collen <- ncol(vcf2gtfiltered3)
rowlen <-nrow(vcf2gtfiltered3)

vcf2gtfiltered3 <- vcf2gtfiltered3[,-collen]
vcf2gtfiltered3 <- vcf2gtfiltered3[-rowlen,]

vcf2gtfiltered3[1:5,550:554]
vcf2gtfiltered3[2762:2763,550:554]



# Next is vcf2@fix

vcf2fix <- vcf2@fix
vcf2fix <- as.data.frame(vcf2fix)

vcf2fix2 <- vcf2fix%>%
  filter((POS %in% snps.rd6$POS))

vcf2fix2

# weird there is a double here
# Need to remove the two duplicate positions!

vcf2fix2

vcf2fix2 = vcf2fix2[!duplicated(vcf2fix2$POS),]

#
#  vcf2gtfiltered3 is the filtered vcf@gt
# vcf2fix2 is the filtered vcf@fix
vcf2gtfiltered3
vcf2fix2

vcf2gtfiltered3 <- as.matrix(vcf2gtfiltered3)
vcf2fix2 <- as.matrix(vcf2fix2)

vcf2@fix <- vcf2fix2
vcf2@gt <- vcf2gtfiltered3

# Need to remove names from rows and columns still

rownames(vcf@fix)
rownames(vcf2@fix) <- NULL
rownames(vcf@gt)
rownames(vcf2@gt) <- NULL

colnames(vcf@fix)
colnames(vcf@gt)
colnames(vcf2@gt)
setwd("D:/Dropbox/RA/Kat/GP bioinformatics/Filtering_steps/2 filter sequencing depth difference")

# This appeared to have worked!
# And the writing of the vcf worked to!
?write.vcf()
write.vcf(vcf2, file = "GP_filtered_step2_DPdiff5x.vcf.gz")
#write.vcf(vcf, file = "GP_raw_vcf.vcf.gz")

nrow(vcf2)
length(colnames(vcf2@gt))









# heterozygosity + techrep and reproducibility (when calculated later)





#Heterozygosity
snps.rd2 <- snps.rd
snps.rd2$na <- rep(NA)
snps.rd2$na <- apply(snps.rd2, 1, function(x) sum(is.na(x)))
snps.rd2$seq <- rep(NA)
snps.rd2$seq <- ncol(snps.rd2) - snps.rd2$na - 8
snps.rd2$hets <- apply(snps.rd2, 1, function(x) length(which(x == "0/1")))
het_count.rd <- 100*(snps.rd2$hets/snps.rd2$seq)
hist(het_count.rd, main="Heterozygosity", xlab="Proportion of heterozygotes at SNP")



# Reproducibility
#create a csv file of 1 column with the names of all sampeles which have or are technical replicates

tech.reps <- read.csv('technicalreplicatesordered.csv', header=F)
snp.tech.reps.rd <- snps.rd[,which(colnames(snps.rd) %in% tech.reps$V1)]
head(colnames(snp.tech.reps.rd))
snp.tech.reps.rd <- snp.tech.reps.rd[,order(names(snp.tech.reps.rd))]
head(colnames(snp.tech.reps.rd))
snp.tech.match.rd <- matrix(nrow=nrow(snp.tech.reps.rd), ncol=ncol(snp.tech.reps.rd))


for(r in 1:nrow(snp.tech.reps.rd)){
  for(c in seq(1,ncol(snp.tech.reps.rd), 2)) {
    if(!is.na(snp.tech.reps.rd[r,c]) && !is.na(snp.tech.reps.rd[r,(c+1)]) && snp.tech.reps.rd[r,c] != snp.tech.reps.rd[r,c+1]){
      snp.tech.match.rd[r,c] <- "ERROR"
    } else {
      snp.tech.match.rd[r,c] <- NA
    }
  }
}

snp.tech.match.rd <- as.data.frame(as.matrix(snp.tech.match.rd))
snp.tech.match.rd$error <- rep(NA)
n4 <- ncol(snp.tech.match.rd) - 1
for (r in 1:nrow(snp.tech.match.rd)) {
  snp.tech.match.rd$error[r] <- length(which(snp.tech.match.rd[r,1:n4] == "ERROR"))
}

snp.tech.match.rd$reproducibility <- rep(NA)
n5 <- (ncol(snp.tech.match.rd) - 2)/2
snp.tech.match.rd$reproducibility <- 100-((snp.tech.match.rd$error/n5)*100)
summary(snp.tech.match.rd$reproducibility)
summary(snp.tech.match.rd$error)
error.rate <- 100-snp.tech.match.rd$reproducibility
summary(error.rate) #mean error rate pre-filtering (aside from read depth filtering)
hist(snp.tech.match.rd$reproducibility, main="Reproducibility between technical replicates (pairs)", cex.main=0.8, xlab="Reproducibility %")


#coverage.rd2 <- coverage.rd1[which(coverage.rd1$read.depth.snp.avg > 2.5),]
#filters by Rep avg

snp.tech.reps.rd2$reproducibility <- snp.tech.match.rd$reproducibility

snps.rd3.index <- subset(snp.tech.reps.rd2, snp.tech.reps.rd2$reproducibility >= 99)

snps.rd3 <- snps.rd2[which(snps.rd2[,1] %in% snps.rd3.index$index),]

