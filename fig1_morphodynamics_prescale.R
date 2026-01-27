setwd("~/Documents/Courses/QuadBio/Morphodynamics")
library(Seurat)
library(dplyr)
library(simspec)

# Preprocessing of the data

# Loading the data
d5 <- readRDS("data/Seurat_object_raw_HB4_D5.rds")
d7 <- readRDS("data/Seurat_object_raw_HB4_D7.rds")
d11 <- readRDS("data/Seurat_object_raw_HB4_D11.rds")
d16 <- readRDS("data/Seurat_object_raw_HB4_D16.rds")
d21 <- readRDS("data/Seurat_object_raw_HB4_D21.rds")
d30 <- readRDS("data/Seurat_object_raw_HB4_D30.rds")

# Quality control
objs <- list(d5, d7, d11, d16, d21, d30)
# to match the QC performed in the paper, I will also nFeature_RNA to be between 1000 to 7,500 and the mitochondiral gnes threshold to be less than 10%
objs <- lapply(objs, function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT[-\\.]")
  subset(x, subset = nFeature_RNA >1000 & nFeature_RNA <7500 & percent.mt <10)
})


#VlnPlot(objs[[1]], features = c("nFeature_RNA", "percent.mt"), ncol = 2)



# to make sure im on the right track, i also count the number of cells
#lapply(objs, ncol)


# Normalization

objs <- lapply(objs, NormalizeData)


#Merging
ntimecourse <- merge(
  objs[[1]],
  y = objs[2:6]
)


# Finding variable features
ntimecourse <- FindVariableFeatures(ntimecourse, nfeatures = 3000)

# making the genes_use file
var_genes <- VariableFeatures(ntimecourse)
tfs <-readLines("transcription_factors.txt")
cc_genes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
genes_use_raw <- union(var_genes,tfs)
genes_use <- setdiff(genes_use_raw,cc_genes)

# computing cell cyle score
ntimecourse<-CellCycleScoring(ntimecourse,
                              s.features = cc.genes.updated.2019$s.genes,
                              g2m.features = cc.genes.updated.2019$g2m.genes,
                             set.ident = FALSE )

# Scaling data
ntimecourse<-ScaleData(ntimecourse, 
                       features = genes_use,
                       vars.to.regress = c("S.Score", "G2M.Score"))

# PCA time! 
ntimecourse <- RunPCA(ntimecourse, features = genes_use, npcs = 50)


# integration with CSS time!

ntimecourse <- cluster_sim_spectrum(ntimecourse, label_tag = "orig.ident", 
                                    cluster_resolution = 0.3, dims_use = 1:10)


#umap
ntimecourse <- RunUMAP(ntimecourse, 
                       reduction = "css", 
                       dims = 1:ncol(Embeddings(ntimecourse, "css")),
                       spread = 1.0,
                       min.dist = 0.4)
#clustering
ntimecourse <- FindNeighbors(ntimecourse, reduction = "css", 
                             dims = 1:ncol(Embeddings(ntimecourse, "css"))) %>%
  FindClusters(resolution = 0.6)

plot1 <- UMAPPlot(ntimecourse, group.by="orig.ident")
plot2 <- UMAPPlot(ntimecourse, label = T)


saveRDS(ntimecourse, file="My RDS Files/ntimecourse.rds")











# Regress out the cell cyle scores from the CSS embeddings
css_embeddings <- Embeddings(ntimecourse, "css")

# Creating a matrix to store the residuals
css_regressed <- matrix(nrow = nrow(css_embeddings),
                        ncol = ncol(css_embeddings))

# Regresssing out S.Score and G2M.Score from each CSS dimension
for (i in 1:ncol(css_embeddings)){
  model <- lm(css_embeddings[,i] ~ ntimecourse$S.Score + ntimecourse$G2M.Score)
  css_regressed[,i] <- residuals(model)
}

# Preserve row and column names
rownames(css_regressed) <- rownames(css_embeddings)
colnames(css_regressed) <- colnames(css_embeddings)

# Replace the CSS reduction with the regressed version
ntimecourse[["css"]] <- CreateDimReducObject(
  embeddings = css_regressed,
  key = "CSS_",
  assay = DefaultAssay(ntimecourse)
)











# Finding variable features

objs <- lapply(objs, FindVariableFeatures, nfeatures = 3000)


# getting the variable genes first
var_genes_list <- lapply(objs, VariableFeatures)
# getting the union of all variable genes
var_genes <- unique(unlist(var_genes_list))
# loading the list of tf from the text file i provided
tfs <-readLines("transcription_factors.txt")
# getting cell cycle genes from Seurat
cc_genes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
#combining and deleting
genes_use_raw <- union(var_genes,tfs)
genes_use <- setdiff(genes_use_raw,cc_genes)

# computing cell cycle scores

objs <- lapply(objs, function(x) {
  CellCycleScoring(
    object = x,
    s.features = cc.genes.updated.2019$s.genes,
    g2m.features = cc.genes.updated.2019$g2m.genes,
    set.ident = FALSE
  )
})

# Z-scaling and regressing out cell cycle scores
objs <- lapply(objs, function(x){
  ScaleData(
    x,
    features = genes_use,
    vars.to.regress = c("S.Score", "G2M.Score")
  )
})

# PCA time!
objs <- lapply(objs, function(x){
  RunPCA(
    x,
    features = genes_use,
    npcs = 50,
    verbose = FALSE
  )
})


## Integration of different time points


# meta cells stuff, i will come back to this at the end (i think this
#happens before the cell scorring thing)


#Merging
ntimecourse <- merge(
  objs[[1]],
  y = objs[2:6]
)

# trying something, perhaps i need to run pCA after the merging
ntimecourse <- RunPCA(ntimecourse, features = genes_use, npcs = 50)

# Data integration using CSS

library(simspec)

ntimecourse <- cluster_sim_spectrum(ntimecourse, label_tag = "orig.ident", 
                                    cluster_resolution = 0.3, dims_use = 1:10)

# i skip the regression for now

#umap
ntimecourse <- RunUMAP(ntimecourse, reduction = "css", 
                       dims = 1:ncol(Embeddings(ntimecourse, "css")))
#clustering
ntimecourse <- FindNeighbors(ntimecourse, reduction = "css", 
                             dims = 1:ncol(Embeddings(ntimecourse, "css"))) %>%
  FindClusters(resolution = 0.6)












# getting the list of TF from AnimalTFDB

tfdb <- read.delim("Homo_sapiens_TF.txt", header = TRUE)
tfs <- tfdb$Symbol
# Remove blanks
tfs <- tfs[tfs!=""]

#saving into a new file
writeLines(tfs,"transcription_factors.txt")


# displaying the top 20 variable features
top_features <- head(VariableFeatures(d5),20)
plot1 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot2


# Linear dimensionality reduction using PCA
# number of PCA set to 50 as in the paper

d5 <- RunPCA(d5, npcs = 50)
d7 <- RunPCA(d7, npcs = 50)
d11 <- RunPCA(d11, npcs = 50)
d16 <- RunPCA(d16, npcs = 50)
d21 <- RunPCA(d21, npcs = 50)
d30 <- RunPCA(d30, npcs = 50)

ElbowPlot(d5, ndims = ncol(Embeddings(d5, "pca")))

