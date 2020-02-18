setwd("~/R/eDNA/")

library(vegan)
#library(DESeq2) #for rlog()
#library(ggplot2)
browseVignettes("vegan")
help("rarefy")
data("BCI")

specnumber(BCI, MARGIN = 2) # MARGIN = 2, it finds frequencies of species
specnumber(BCI, MARGIN = 1) # observed number of species per row

S <- specnumber(BCI, MARGIN = 1) # returns the number of species observed in a row! 
raremax <- min(specnumber(BCI, MARGIN = 1))
Srare <- rarefy(BCI, raremax)

plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species") #77
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)

####################################################
quantile(rowSums(BCI))

spa <- specaccum(BCI)
plot(spa)

mod <- decorana(BCI)
plot(mod)

#To express richness for the same number of individ-uals, we can use:
Srar <- rarefy(BCI, min(rowSums(BCI)))
S100 <- rarefy(BCI, 100)
S2 <- rarefy(BCI,2)

par(mfrow = c(2,2))
plot(S, Srar, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species") #340
plot(S, S100, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species") #100
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species") #77
plot(S, S2, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species") #2



rarecurve(BCI, step = 1, sample = raremax, col = "blue", cex = 0.6)
rarecurve(BCI, step = 10, sample = raremax, col = "blue", cex = 0.6)
rarecurve(BCI, step = 100, sample = raremax, col = "blue", cex = 0.6)
# horizontal lines for the rarefied species richnesses!
rarecurve(BCI[1,], step = 10, sample = raremax, col = "blue", cex = 0.6) #for line 1
rarecurve(BCI[c(1:5),], step = 10, sample = raremax, col = "blue", cex = 0.6) #for line 1 -5


###############################
df <- read.csv2("eDNA_counts.csv", row.names = 1)
df[is.na(df)] <- 0 # replace NAs to 0
df.t <- as.data.frame(t(df)) #transpose

#To express richness for the same number of individ-uals, we can use:
raremax <- min(rowSums(df.t))
S <- specnumber(df.t, MARGIN = 1) # returns the number of species observed in a row! 

Srar <- rarefy(df.t, min(rowSums(df.t)))
S100 <- rarefy(df.t, 100)
S2 <- rarefy(df.t,2)

plot(S, Srar, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
plot(S, S100, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
plot(S, S2, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

rarecurve(df.t, step = 10, sample = raremax, col = "blue", cex = 0.6)
rarecurve(df.t[5,], step = 10, sample = raremax, col = "blue", cex = 0.6)


# rlog transform
rlog(df, blind=F)
# normTransform
normTransform(df)


