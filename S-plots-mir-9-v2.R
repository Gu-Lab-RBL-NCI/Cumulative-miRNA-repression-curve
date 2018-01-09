setwd("~/Desktop")
#install.packages("dplyr")
library(dplyr)

require(gdata)
#targets <- read.xls ("TargetScan7.1__miR-9-5p.predicted_targets.xlsx", sheet = "TargetScan7", header = TRUE) #all targets
targets <- read.xls ("TargetScan_miR-9-5p.predicted_targets.xlsx", sheet = "TargetScan5", header = TRUE) #all targets 5.2 version
#targets <- read.xls ("conserved-canonical-targets-mir-9.xlsx", sheet = "Sheet1", header = TRUE) #conserved targets custom 5.2v
alternative <- read.xls ("TargetScanCustom5.2_miR9-5p-alternative.xlsx", sheet = "Sheet1", header = TRUE)
colnames(targets)[1] <- "gene"

#expression <- read.csv("GBMLGG_EB_RmDiffFullGenesRanRmDup.csv", header = TRUE)
expression <- read.table("LGG-GBM.gene_expression.normalized.txt", header = TRUE)
miRNA <- read.table("summary2-mir-9-cleavage-LGG.txt", header = TRUE)

survival <- read.xls ("1486946718491-manifest.xlsx", sheet = "Sheet1", header = TRUE)
survival <- survival[,c(2,26)]
colnames(survival)[1] <- "ID"
survival$ID <- substr(survival$ID, 1, 28)
miRNA <- merge(miRNA, survival, by="ID")

miRNA$samples <- substr(miRNA$ID, 1, 12)
index <- which(duplicated(miRNA$samples))
miRNA <- miRNA[-index, ]

expression$samples<-rownames(expression)
row.names(expression) <- NULL

total <- merge(miRNA, expression, by = "samples")

genes_survival <- total
genes_survival <- genes_survival[complete.cases(genes_survival),]
genes_survival_miRNA <- genes_survival[,c(1:20)]
#write.table(genes_survival_miRNA, "genes_survival_miRNA.txt", sep="\t", append = FALSE)

genes_survival_gene <-as.data.frame(cbind(genes_survival$days_to_death, genes_survival$MAML1))
genes_survival_gene <- genes_survival_gene[order(genes_survival_gene$V2),] 
#write.table(genes_survival_gene, "genes_survival_gene.txt", sep="\t", append = FALSE)

total <- total[,-c(1,2)]

hist(total$RPMalternative, breaks = 30)
hist(total$RPMcanonical, breaks = 60)
summary(total[,c(16,17)])

#Parameters
n <- 25 #number of samples
threshold <- 200
#threshold <- 200

#total <- total[order(total$RPMcanonical),] 
total <- total[order(total$RPMalternative),] 
lowalternative <- total[c(2:(n+1)), ]
highalternative <- total[(nrow(total)-n+1):nrow(total), ]
diff_A5P <- log2(mean(highalternative$RPMalternative)) - log2(mean(lowalternative$RPMalternative))

lowcanonical <- total[(1:n), ]
highcanonical <- total[(nrow(total)-n+1):nrow(total), ]
diff_C5P <- log2(mean(highcanonical$RPMcanonical)) - log2(mean(lowcanonical$RPMcanonical))

#n=124 25% sample
#lowalternative2 <- total[ which(total$RPMalternative < 16928.6),]
#highalternative2 <- total[ which(total$RPMalternative > 26179.0),]
#lowcanonical <- total[ which(total$RPMcanonical < 535240),]
#highcanonical <- total[ which(total$RPMcanonical > 724121),]

numcol <- ncol(total)

### Restart
grouped <- data.frame(matrix(ncol = numcol, nrow =0))
colnames(grouped) <- colnames(total)

grouped <- rbind(sapply(lowalternative, mean), sapply(highalternative, mean), sapply(lowcanonical, mean), sapply(highcanonical, mean))
grouped <-t(grouped)
grouped <-as.data.frame(grouped)
colnames(grouped)[1] <- "Low_A5P"
colnames(grouped)[2] <- "High_A5P"
colnames(grouped)[3] <- "Low_C5P"
colnames(grouped)[4] <- "High_C5P"

grouped <- grouped[-c(1:17),]
grouped <-grouped[ which(grouped$Low_A5P >threshold & grouped$Low_C5P >threshold & grouped$High_A5P >threshold & grouped$High_C5P >threshold ),]

grouped$Log2FC_A5P <- log2(grouped$High_A5P) - log2(grouped$Low_A5P)
grouped$Log2FC_C5P <- log2(grouped$High_C5P) - log2(grouped$Low_C5P)

grouped$gene<-rownames(grouped)
rownames(grouped) <- NULL
targets <- targets[which(targets$Conserved.sites.total>0),] #filter-out targets poorly conserved only
intersect_grouped_alternative <- intersect(grouped[,c(7)],alternative[,c(1)])
intersect_grouped_canonical <- intersect(grouped[,c(7)], targets[,c(1)])
intersect_alternative_canonical <- intersect(alternative[,c(1)], targets[,c(1)])
newgrouped <- grouped[ !(grouped$gene %in% intersect_grouped_canonical), ]
newgrouped <- newgrouped[ !(newgrouped$gene %in% intersect_grouped_alternative), ]
purecanonical <- targets[ !(targets$gene %in% intersect_alternative_canonical), ]
purealternative <- alternative[ !(alternative$gene %in% intersect_alternative_canonical), ]

both <- targets[ (targets$gene %in% intersect_alternative_canonical), ]

breaks <-  seq(-1.5, 1.5, by=0.05)

##Canonical targets
#Baseline genes Canonical
datLog2FC <-  as.numeric(grouped$Log2FC_C5P) #canonical FC
cumfreq0 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(grouped))) 

plot(breaks, cumfreq0, type="n", xlab = "mRNA fold-change Log2", ylab = "Cumulative fraction miR-9 canonical targets")
lines(breaks, cumfreq0)

datLog2FC <-  as.numeric(newgrouped$Log2FC_C5P) #canonical FC without targets
cumfreq1 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(newgrouped))) 
lines(breaks, cumfreq1, col="pink")


canonicalmachedtargets <- merge(grouped, purecanonical, by="gene")
#write.table(canonicalmachedtargets, "canonicalmachedtargets.txt", sep="\t", append = FALSE, row.names = F)
canonicalmachedtargets <- merge(grouped, targets, by="gene")

#one target site
canonical1 <- canonicalmachedtargets[ which(canonicalmachedtargets$Conserved.sites.total==1),]
canonical1_line <- c(0, cumsum(table(cut(canonical1$Log2FC_C5P, breaks, right=FALSE))/nrow(canonical1)))
lines(breaks, canonical1_line, col="blue")

#two target site
canonical2 <- canonicalmachedtargets[ which(canonicalmachedtargets$Conserved.sites.total==2),]
canonical2_line <- c(0, cumsum(table(cut(canonical2$Log2FC_C5P, breaks, right=FALSE))/nrow(canonical2)))
lines(breaks, canonical2_line, col="red")

#three or more target site
canonical3 <- canonicalmachedtargets[ which(canonicalmachedtargets$Conserved.sites.total>=3),]
canonical3_line <- c(0, cumsum(table(cut(canonical3$Log2FC_C5P, breaks, right=FALSE))/nrow(canonical3)))
lines(breaks, canonical3_line, col="green")

#write.table(rbind(breaks, cumfreq0, cumfreq1, canonical1_line, canonical2_line, canonical3_line), "mir-9-5p-canonical-1-2-3.txt", sep="\t", append = FALSE)

diff_C5P
diff_A5P
##Alternative targets

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

write.table(rbind(breaks, cumfreq0, cumfreq1, alternative1_line, alternative2_line, alternative3_line), "mir-9-5p-alternative-1-2-3.txt", sep="\t", append = FALSE)
write.table(rbind(breaks, cumfreq0, cumfreq1, alternative1_line, alternative2_line, alternative3_line), "mir-9-5p-alternative-1-2.txt", sep="\t", append = FALSE, row.names = FALSE)

#both canonical and alternative
datLog2FC <-  as.numeric(grouped$Log2FC_A5P) #alternative FC
cumfreq0 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(grouped))) 
plot(breaks, cumfreq0, type="n", xlab = "mRNA fold-change Log2", ylab = "Cumulative fraction miR-9 canonical + alterantive targets")
lines(breaks, cumfreq0)

datLog2FC <-  as.numeric(newgrouped$Log2FC_A5P) #alternative FC without targets
cumfreq1 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(newgrouped))) 
lines(breaks, cumfreq1, col="pink")

bothmachedtargets <- merge(grouped, both, by="gene")

#one target site
both1 <- bothmachedtargets[ which(bothmachedtargets$Conserved.sites.total==1),]
both1_line <- c(0, cumsum(table(cut(both1$Log2FC_A5P, breaks, right=FALSE))/nrow(both1)))
lines(breaks, both1_line, col="blue")
alternative1 <- alternativemachedtargets[ which(alternativemachedtargets$total==1),]
alternative1_line <- c(0, cumsum(table(cut(alternative1$Log2FC_A5P, breaks, right=FALSE))/nrow(alternative1)))
lines(breaks, alternative1_line, col="red")
canonical1 <- canonicalmachedtargets[ which(canonicalmachedtargets$Conserved.sites.total==1),]
canonical1_line <- c(0, cumsum(table(cut(canonical1$Log2FC_A5P, breaks, right=FALSE))/nrow(canonical1)))
lines(breaks, canonical1_line, col="green")
#write.table(rbind(breaks, cumfreq0, cumfreq1, canonical1_line, alternative1_line, both1_line), "mir-9-5p-both.txt", sep="\t", append = FALSE)

###############

#one target site 8mer
#alt8mer <- alternativemachedtargets[ which(alternativemachedtargets$X8mer==1),]
#alt8mer <- alt8mer[order(alt8mer$Log2FC_A5P),] 
#alt8mer_line <- c(0, cumsum(table(cut(alt8mer$Log2FC_A5P, breaks, right=FALSE))/nrow(alt8mer)))
#lines(breaks, alt8mer_line, col="pink")

#one target site 7mer-m8
#alt7merm8 <- alternativemachedtargets[ which(alternativemachedtargets$X7mer.m8==1),]
#alt7merm8 <- alt7merm8[order(alt7merm8$Log2FC_A5P),] 
#alt7merm8_line <- c(0, cumsum(table(cut(alt7merm8$Log2FC_A5P, breaks, right=FALSE))/nrow(alt7merm8)))
#lines(breaks, alt7merm8_line, col="green")

#one target site 7mer-A1
#alt7mera1 <- alternativemachedtargets[ which(alternativemachedtargets$X7mer.1A==1),]
#alt7mera1 <- alt7mera1[order(alt7mera1$Log2FC_A5P),] 
#alt7mera1_line <- c(0, cumsum(table(cut(alt7mera1$Log2FC_A5P, breaks, right=FALSE))/nrow(alt7mera1)))
#lines(breaks, alt7mera1_line, col="purple")

#candiategenes <- rbind(alt7mera1[c(1:15),], alt7merm8[c(1:15),], alt8mer[c(1:15),])
#write.table(candiategenes, "mir-9-5p-alternative-candidate-genes2.txt", sep="\t", append = FALSE)
#grouped2 <- data.frame(matrix(ncol = numcol, nrow =0))
#colnames(grouped2) <- colnames(total)
#grouped2 <- rbind(sapply(lowalternative2, mean), sapply(lowalternative2, sd), sapply(highalternative2, mean), sapply(highalternative2, sd))
#grouped2 <-t(grouped2)
#grouped2 <-as.data.frame(grouped2)
#grouped2$gene<-rownames(grouped2)
#rownames(grouped2) <- NULL
#colnames(grouped2)[1] <- "Low_A5P_mean"
#colnames(grouped2)[2] <- "Low_A5P_sd"
#colnames(grouped2)[3] <- "High_A5P_mean"
#colnames(grouped2)[4] <- "High_A5P_sd"
#grouped2$Log2FC_A5P <- log2(grouped2$High_A5P_mean) - log2(grouped2$Low_A5P_mean)

#grouped3 <- merge(grouped2, candiategenes, by="gene")
#write.table(grouped3, "mir-9-5p-alternative-candidate-genes2-1-4-quartile.txt", sep="\t", append = FALSE)


#write.table(rbind(breaks, cumfreq0, alt8mer_line, alt7merm8_line, alt7mera1_line), "mir-9-5p-alternative-8-7m8-7a1.txt", sep="\t", append = FALSE)
