# Cumulative miRNA repression curve

This script provides the code to generate a cumulative miRNA repression curve. To investigate whether miR-9-alt regulates predicted miR-9-alt-specific targets in LGG, we took advantage of the RNA-Seq and miRNA-Seq data in The Cancer Genome Atlas (TCGA). Specifically, we compared mRNA expression profiles between patients with high levels of miR-9-alt (top 5%, 25 samples) and relatively low levels of miR-9-alt (bottom 5%, 25 samples).

##### Diagram of the analysis pipeline used to obtain the Cumulative curves
<img src="https://github.com/Gu-Lab-RBL-NCI/Cumulative-miRNA-repression-curve/blob/master/scheme.png">

## The current scripts for the generation of the cumulative miRNA repression curve execute the follow:

**1. Import miRNA targets predicted from [TargetScan](http://www.targetscan.org//vert_50/seedmatch.html) [(Lewis et al., 2005)](https://www.ncbi.nlm.nih.gov/pubmed/15652477)**

```
targets <- read.xls ("TargetScan_miR-9-5p.predicted_targets.xlsx", sheet = "TargetScan5", header = TRUE) #all targets 5.2 version
alternative <- read.xls ("TargetScanCustom5.2_miR9-5p-alternative.xlsx", sheet = "Sheet1", header = TRUE)
```
 - Filter-out targets poorly conserved genes
 - Establish unique miR-9-can and miR-9-alt target genes

```
targets <- targets[which(targets$Conserved.sites.total>0),] #filter-out targets poorly conserved only
intersect_grouped_alternative <- intersect(grouped[,c(7)],alternative[,c(1)])
intersect_grouped_canonical <- intersect(grouped[,c(7)], targets[,c(1)])
intersect_alternative_canonical <- intersect(alternative[,c(1)], targets[,c(1)])
newgrouped <- grouped[ !(grouped$gene %in% intersect_grouped_canonical), ]
newgrouped <- newgrouped[ !(newgrouped$gene %in% intersect_grouped_alternative), ]
purecanonical <- targets[ !(targets$gene %in% intersect_alternative_canonical), ]
purealternative <- alternative[ !(alternative$gene %in% intersect_alternative_canonical), ]
```

**2. Import the normalized mRNA-seq and miR-9-can/alt expression from LGG patients**

```
expression <- read.table("LGG-GBM.gene_expression.normalized.txt", header = TRUE)
miRNA <- read.table("summary2-mir-9-cleavage-LGG.txt", header = TRUE)
```

*Note that the file `LGG-GBM.gene_expression.normalized.txt` needs to be downloaded separately from [The Cancer Genome Atlas datasets](https://tcga-data.nci.nih.gov/docs/publications/lgggbm_2016/).*

**3. Import LGG patients data manifest**

**4. Set up high/low miRNA level expressing groups and minimum expression of the genes**
```
#Parameters
n <- 25 #number of samples
threshold <- 200 #threshold, minimum gene expression
```

**5. Calculate for each gene the fold-change between the high/low miRNA level expressing groups**

```
total <- total[order(total$RPMalternative),] 
lowalternative <- total[c(2:(n+1)), ]
highalternative <- total[(nrow(total)-n+1):nrow(total), ]
diff_A5P <- log2(mean(highalternative$RPMalternative)) - log2(mean(lowalternative$RPMalternative))

lowcanonical <- total[(1:n), ]
highcanonical <- total[(nrow(total)-n+1):nrow(total), ]
diff_C5P <- log2(mean(highcanonical$RPMcanonical)) - log2(mean(lowcanonical$RPMcanonical))
```

**6. Calculate the cumulative fold-change for the base line genes (no targets) and genes with targets for miR-9-alt**

```
##Alternative targets
#Baseline genes 
datLog2FC <-  as.numeric(grouped$Log2FC_A5P) #alternative FC
cumfreq0 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(grouped))) 
plot(breaks, cumfreq0, type="n", xlab = "mRNA fold-change Log2", ylab = "Cumulative fraction miR-9 alterantive targets")
lines(breaks, cumfreq0)

datLog2FC <-  as.numeric(newgrouped$Log2FC_A5P) #alternative FC without targets
cumfreq1 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(newgrouped))) 
lines(breaks, cumfreq1, col="pink")

alternativemachedtargets <- merge(grouped, purealternative, by="gene")
#alternativemachedtargets <- merge(grouped, alternative, by="gene")
write.table(alternativemachedtargets, "alternativemachedtargetsv2.txt", sep="\t", append = FALSE, row.names = F)

#one target site
alternative1 <- alternativemachedtargets[ which(alternativemachedtargets$total==1),]
alternative1_line <- c(0, cumsum(table(cut(alternative1$Log2FC_A5P, breaks, right=FALSE))/nrow(alternative1)))
lines(breaks, alternative1_line, col="blue")

#two target site
alternative2 <- alternativemachedtargets[ which(alternativemachedtargets$total>1),]
alternative2_line <-c(0, cumsum(table(cut(alternative2$Log2FC_A5P, breaks, right=FALSE))/nrow(alternative2)))
lines(breaks, alternative2_line, col="red")

#3 or more target site --> too few targets
alternativemore3 <- alternativemachedtargets[ which(alternativemachedtargets$X8mer==1 & alternativemachedtargets$total==1),]
alternative3_line <-c(0, cumsum(table(cut(alternativemore3$Log2FC_A5P, breaks, right=FALSE))/nrow(alternativemore3)))
lines(breaks, alternative3_line, col="green")
```


### References:

[Lewis, B.P., Burge, C.B., and Bartel, D.P. (2005). Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microRNA targets. Cell 120, 15–20.](https://www.ncbi.nlm.nih.gov/pubmed/15652477)
