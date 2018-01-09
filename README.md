# Cumulative miRNA repression curve

This script provides the code to generate a cumulative miRNA repression curve. To investigate whether miR-9-alt regulates predicted miR-9-alt-specific targets in LGG, we took advantage of the RNA-Seq and miRNA-Seq data in The Cancer Genome Atlas (TCGA). Specifically, we compared mRNA expression profiles between patients with high levels of miR-9-alt (top 5%, 25 samples) and relatively low levels of miR-9-alt (bottom 5%, 25 samples).

##### Diagram of the analysis pipeline used to obtain the Cumulative curves
<img src="https://github.com/Gu-Lab-RBL-NCI/Cumulative-miRNA-repression-curve/blob/master/scheme.png">

## The current scripts for the generation of the cumulative miRNA repression curve execute the follow:

**1. Import miRNA targets predicted from [TargetScan](http://www.targetscan.org//vert_50/seedmatch.html) [(Lewis et al., 2005)](https://www.ncbi.nlm.nih.gov/pubmed/15652477)**

**2. pri-miRNA secondary structure**


*Note that the file `LGG-GBM.gene_expression.normalized.txt` needs to be downloaded separately from [The Cancer Genome Atlas datasets](https://tcga-data.nci.nih.gov/docs/publications/lgggbm_2016/).


Cumulative fraction plot of fold-change in expression of mRNAs between the top and low levels of the miR-9-alt in patients from LGG
