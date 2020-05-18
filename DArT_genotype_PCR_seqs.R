# The purpose of this code is to identify SNP sequences that distinguish between
# the four E.diaphana genotypes maintained at the University of Melbourne so that
# they might be used in a PCR assay to identify the morphologically identical
# animals.

library(dplyr)
library(tibble)
setwd("~/repositories/DArT_SNP")


# Import raw DArT data as a reference dataframe
# Hyphens have been recoded as NA to prevent SNP codes converting to factor levels
snpSeqs <- read.csv(file="Report_DExa16-2399_SNP_singlerowNA.csv", header = T, skip = 6, row.names = 1)

# Subset the DArT report to remove the AIMS SeaSim samples and sample 15
# Also omit columns with irrelevant DArT data
snpSeqsSub <- snpSeqs[,c(2, 9, 20:32, 34, 36, 38, 40:48)]

# Convert dataframe to tibble for filtering
snpSeqsSubTbl <- as_tibble(snpSeqsSub)

# Convert rownames to values of variable called SnpDesc to keep dplyr happy
snpSeqsSubTbl <- rownames_to_column(snpSeqsSubTbl, var = "SnpDesc")

# Remove loci (rows) wih reproducibility <100%
snpSeqsSubTbl <- filter(snpSeqsSubTbl, RepAvg == 1)

# Remove loci (rows) with missing values (call rate < 1)
snpSeqsSubTbl <- filter(snpSeqsSubTbl, CallRate == 1)

# Subset samples by genotype
ep1 <- select(snpSeqsSubTbl, c(1,2,9,21,22,23,15,17,26,13))
ep2 <- select(snpSeqsSubTbl, c(1,2,8,6,24,10,16))
ep3 <- select(snpSeqsSubTbl, c(1,2,12,18,25,19,27,20,7))
ep4 <- select(snpSeqsSubTbl, c(1,2,5,11,14))

# Identify homozygous non-SNP loci in each genotype
# Strict threshold = 0, but final result is no SNPs!
ep1nonsnp <- rowSums(ep1[,3:10]) <= 2
ep2nonsnp <- rowSums(ep2[,3:7]) <= 2
ep3nonsnp <- rowSums(ep3[,3:9]) <= 2
ep4nonsnp <- rowSums(ep4[,3:5]) <= 2

# Retain only SNP loci that discriminate from the others unambiguously
ep1unamb <- ep1[ep2nonsnp == T & ep3nonsnp == T & ep4nonsnp == T,]
ep2unamb <- ep2[ep1nonsnp == T & ep3nonsnp == T & ep4nonsnp == T,]
ep3unamb <- ep3[ep1nonsnp == T & ep2nonsnp == T & ep4nonsnp == T,]
ep4unamb <- ep4[ep1nonsnp == T & ep2nonsnp == T & ep3nonsnp == T,]

# Identify SNP loci with maximal SNP designation
# No loci possess 100% SNP designation, so relax criteria to (no. of samples x 2)-1
ep1snpMax <- rowSums(ep1unamb[,3:10]) >= 15
ep2snpMax <- rowSums(ep2unamb[,3:7]) >= 9
ep3snpMax <- rowSums(ep3unamb[,3:9]) >= 13
ep4snpMax <- rowSums(ep4unamb[,3:5]) >= 5

# Retain only SNP loci with maximal SNP designation
ep1.snp <- ep1unamb[ep1snpMax,]
ep2.snp <- ep2unamb[ep2snpMax,]
ep3.snp <- ep3unamb[ep3snpMax,]
ep4.snp <- ep4unamb[ep4snpMax,]

# Write ouput for each genotype to csv file
write.table(x = ep1.snp, file = "ep1snpData.csv", sep = ",", col.names = NA)
write.table(x = ep2.snp, file = "ep2snpData.csv", sep = ",", col.names = NA)
write.table(x = ep3.snp, file = "ep3snpData.csv", sep = ",", col.names = NA)
write.table(x = ep4.snp, file = "ep4snpData.csv", sep = ",", col.names = NA)
