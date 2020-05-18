library(ade4)
library(adegenet)
library(dartR)
setwd("~/repositories/DArT_SNP")


# Read files into a genlight object
geno1 <- gl.read.dart(filename = "Report_DExa16-2399_SNP_singlerow.csv", ind.metafile = "Epallida_Genotype_Metadata2.csv", topskip = 6)


# Remove A2, A4, A24, A25 and E.p._A_15. Remove monomorphic loci. Recalc stats.
geno2 <- gl.drop.ind(geno1, ind.list = c("A2", "A4", "A24", "A25", "E.p._A_15"), mono.rm = T, recalc = T)
rm(geno1)


# Remove loci with average reproducibility of <100%
geno3 <- gl.filter.repavg(geno2, t = 1)
rm(geno2)


# Retain only loci without missing values
geno4 <- gl.filter.callrate(geno3, method = "loc", threshold = 1)
rm(geno3)


# Generate and plot PCoA to assess delineation of samples
geno4pcoa <- gl.pcoa(geno4)
gl.pcoa.plot(geno4pcoa, geno4, labels = "pop")
rm(geno4)
rm(geno4pcoa)


# Based on clustering in the PCoA, there are 4 genotypes:
# Ep1 == E.p._A01, E.p._A02, E.p._A03, E.p._A07, E.p._A08, E.p._A09, E.p._A12, E.p._A21
# Ep2 == SII, E.p._S01, E.p._A10, E.p._A20, E.p._A22
# Ep3 == E.p._A06, E.p._A11a, E.p._A11b, E.p._A14, E.p._A17, E.p._A18, E.p._A19
# Ep4 == SI, SA, SB
