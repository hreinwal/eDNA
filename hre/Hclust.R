# Bray Curtis

# Pearson Correlation 
df.rld <- read.csv2(file = "df.csv",
                    header = TRUE, 
                    row.names = 1)

dis <- vegdist(1-cor(df.rld, method = "pearson"))

plot(hclust(dis, method = "complete")) # clustering: "single", "complete", "average", "mcquitty", "median" or "centroid"

plot(hclust(dis, method = "complete"), # clustering: "single", "complete", "average", "mcquitty", "median" or "centroid"
     labels = paste(coldata$Condition),
     ylab = "Pearson correlation")

plot(hclust(dis, method = "complete"), # clustering: "single", "complete", "average", "mcquitty", "median" or "centroid"
     labels = paste(coldata$Condition, coldata$Tank, sep="-"),
     ylab = "Pearson correlation")
     

# Spearman Correlation 