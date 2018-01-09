# Cumulative miRNA repression curve

This script provides the code to generate a cumulative miRNA repression curve. To investigate whether miR-9-alt regulates predicted miR-9-alt-specific targets in LGG, we took advantage of the RNA-Seq and miRNA-Seq data in The Cancer Genome Atlas (TCGA). Specifically, we compared mRNA expression profiles between patients with high levels of miR-9-alt (top 5%, 25 samples) and relatively low levels of miR-9-alt (bottom 5%, 25 samples).

##### Diagram of the analysis pipeline used to obtain the Cumulative curves
<img src="https://github.com/Gu-Lab-RBL-NCI/Cumulative-miRNA-repression-curve/blob/master/scheme.png">

## The current scripts for the generation of the cumulative miRNA repression curve execute the follow:

**1. Import miRNA targets predicted from [TargetScan](http://www.targetscan.org//vert_50/seedmatch.html) [(Lewis et al., 2005)](https://www.ncbi.nlm.nih.gov/pubmed/15652477)**
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

**2. Import the normalized mRNA-seq expression from LGG patients**

*Note that the file `LGG-GBM.gene_expression.normalized.txt` needs to be downloaded separately from [The Cancer Genome Atlas datasets](https://tcga-data.nci.nih.gov/docs/publications/lgggbm_2016/).*

**3. Import LGG patients data manifest**

**4. Set up high/low miRNA level expressing groups and minimum expression of the genes**
```
#Parameters
n <- 25 #number of samples
threshold <- 200 #threshold, minimum gene expression
```

**5. Calculate for each gene the fold-change between the high/low miRNA level expressing groups**




### References:

[Lewis, B.P., Burge, C.B., and Bartel, D.P. (2005). Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microRNA targets. Cell 120, 15–20.](https://www.ncbi.nlm.nih.gov/pubmed/15652477)
