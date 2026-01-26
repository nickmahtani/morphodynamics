
This one: ![[Pasted image 20260119092216.png]]

But before I get ahead of myself, I need to pre-process the data. Basically do what this paragraph says:

We used Cell Ranger (v3.0.2) with the default parameters to obtain transcript count matrices by aligning the sequencing reads to the human genome and transcriptome (hg38, provided by 10x Genomics, v3.0.0). Count matrices were further preprocessed using the Seurat R package (v4.3.0) and R (version 4.4.0)[63](https://www.nature.com/articles/s41586-025-09151-3#ref-CR63 "Stuart, T. et al. Comprehensive integration of single-cell data. Cell 177, 1888–1902.e21 (2019)."). First, cells were filtered on the basis of the number of detected genes and the fraction of mitochondrial genes. As sequencing depth varied between time points, the threshold of the number of detected genes was set individually for each sample. For the scRNA-seq time course, the number of detected genes was between 1,000 and 7,500 and the mitochondrial genes threshold was less than 10%. Three thousand variable features were used and the number of PCA was set to 50. For Fig. [1b–d](https://www.nature.com/articles/s41586-025-09151-3#Fig1), the total number of cells per day after preprocessing were day 5 = 5,481, day 7 = 8,183, day 7 = 4,912, day 16 = 6,571, day 21 = 7,962 and day 30 = 7,950.

## Step 1

- Use cell ranger: but is this in R or in bash? I will come back to this later

## Step 2

- Assuming you already got the cell ranger files, generated using 10x genomics and preprossesed using cell ranger
- I downloaded the seurat files [here](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15057) to do some analysis assuming that they contain the files aready being run through cell ranger. 

```

```