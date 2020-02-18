
setwd("Y:/R/eDNA/2017Weil/")
# Load Packages
library(vegan)
library(dplyr)
#browseVignettes("vegan")

#Exp <- as.character("Prossen 2017")
Exp <- as.character("Weil 2017") #Experiment = Exp

### Load Data #####
coldata <- read.csv2("coldata.csv", row.names = 1)
df <- read.csv2("eDNA_counts.csv", row.names = 1)
df[is.na(df)] <- 0 # replace NAs to 0
###################

pdf(file = paste0(Exp,"_eDNA_analysis.pdf"), 
    #paper = "a4",
    onefile = T,
    bg = "transparent", #Background color
    fg ="black",        #Foreground color
    width = 15.6,
    height = 9.2) 

par(mfrow=c(2,2))

### Plot Reads (seqzise) ####
seqsize <- apply(df, 2, sum)
xx <- barplot(seqsize,
              ylim = c(0,1.2*max(seqsize)),
              ylab = "Number of read counts",
              #xlab = "",
              main = paste0("Read counts - ",Exp)
        )
text(x=xx, y=seqsize, labels = seqsize, pos = 3, cex = 0.8, col = "red")
min(seqsize) # 18446
############################

### DESeq2 Normalization  ##########
library(DESeq2)

all(colnames(df) == rownames(coldata))
match(colnames(df),row.names(coldata))

# ~ Specimen
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = coldata,
                              design = ~ Specimen)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
norm.counts <- counts(dds, normalized=TRUE)

seqsize <- apply(norm.counts, 2, sum)
xx <- barplot(seqsize,
              ylim = c(0,1.2*max(seqsize)),
              ylab = "DESeq2 norm read counts",
              #xlab = "",
              main = paste0("DESeq2 norm counts - ",Exp," [~ SampleType + Specimen]")
)
seqsize <- round(seqsize, digits = 1)
text(x=xx, y=seqsize, labels = seqsize, pos = 3, cex = 0.8, col = "red")
normCounts <- as.data.frame(norm.counts)
normCounts.t <- as.data.frame(t(norm.counts))
normCounts.t.no0 <- normCounts.t[,colSums(normCounts.t)>0]

min(seqsize) # 31332.5
X <- c(18446,min(seqsize))
xx <- barplot(X,
        ylab = "Minimum Read count size",
        main = "Comparison min. read count size Unnormalized vs DESeq2 Norm.",
        ylim = c(0,1.2*max(X))
        )
text(x=xx, y=X, 
     labels = c(paste0("Unnorm - ",X[1]),paste0("DESeq2 norm - ",X[2])), 
     pos = 3, cex = 0.8, col = "red")
#############################

par(mfrow=c(1,2))
### Continue with Vegan #####
counts <- as.data.frame(t(df)) #transpose
counts.no0 <- counts[,colSums(counts)>0] # Filter data (removing 0 count species)

# Remove contaminations: 
#sort(colnames(counts)) #let's see what we have here ...

CrapCheck <- c("Danio rerio", "Oryzias latipes","Pimephales promelas")
for(CRAP in CrapCheck){
    if (CRAP %in% colnames(counts) == T){
        message(paste(" Oh no!",CRAP,"contamination! :/ => Removing this species from the list."))
        x <- match(CRAP, colnames(counts)) # returns colNbr of Drerio in counts
        counts <- counts[,-x]
        counts.no0 <- counts.no0[,-x]
        rm(x)
    } else {
        message(paste(" Juhu! No",CRAP,"contamination! :)"))
    }
}

### Plot species per sample ####
specNbr <- specnumber(counts.no0)
X <- specNbr
xx <- barplot(X,
              ylim = c(0,1.4*max(X)),
              ylab = "N of observed Species",
              xlab = "Sample Specimen",
              main = paste0("Number of observed Species - ",Exp)
)
text(x=xx, y = X, labels = X, pos = 3, cex = 1, col = "red")
###############################

### Venn #####
X <- counts.no0
S1_1 <- as.data.frame(which(specnumber(X[1,], MARGIN = 2) > 0))
S1_2 <- as.data.frame(which(specnumber(X[2,], MARGIN = 2) > 0))
S1_3 <- as.data.frame(which(specnumber(X[3,], MARGIN = 2) > 0))

x <- list(S1_1 = rownames(S1_1),
          S1_2 = rownames(S1_2),
          S1_3 = rownames(S1_3))
venn <- gplots::venn(x, show.plot = T)
attr(venn,"intersections")$`S1_1:S1_2:S1_3` #List of intersecting species 
###########

### RAREFY ####  ... NOT DONE  YET ... 

### SPECACCUM #####
X <- counts.no0[c(4:nrow(counts.no0)),]
perm <- 1000
sp1 <- specaccum(X, "exact", permutations = perm)
sp2 <- specaccum(X, "random", permutations = perm)

sp2
summary(sp2)

par(mfrow=(c(1,2)))
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     ylab = "Number of observed Species",
     xlab = "Number of eDNA sample replicates",
     main = paste0("Species Accumulation Curve - ",Exp))
boxplot(sp2, col="yellow", add=TRUE, pch="+")
#sp3 <- specaccum(X, "rarefaction") this looks weird ... 
#plot(sp3, col="yellow", add=TRUE, pch="+")

## Fit Arrhenius models to all random accumulations
# WARNING! For many permutations this takes super long! 
mods <- fitspecaccum(sp2, "arrh")
plot(mods, col="hotpink",
     ylab = "Number of observed Species",
     xlab = "Number of eDNA sample replicates",
     main = "Arrhenius models for random accum")
boxplot(sp2, col = "yellow", border = "blue", lty=1, cex=0.3, add= TRUE)####

## Fit Lomolino model to the exact accumulation
mod1 <- fitspecaccum(sp1, "lomolino")
coef(mod1)
fitted(mod1)

## Add Lomolino model using argument 'add'
plot(mod1, add = TRUE, col=2, lwd=2)

## Use nls() methods to the list of models
sapply(mods$models, AIC)
#################

## Accumulation model POOLACCUM #####
par(mfrow=c(1,1))
pool <- poolaccum(X, 
                  permutations = perm,
                  minsize = 2) # estimates the extrapolated species richness in a species pool
summary(pool, display = "chao")
plot(pool,
     main = paste0("Extrapolated species richness models - ",Exp)
)
## Quantitative model
estimateR(X[,])
###################################
dev.off()

save.image(paste0(getwd(),"/.RData"))

# Create relative abundance table
#x <- counts.no0
#apply(x, 2, colSums)  #MARGIN (1=rows, 2=cols)

