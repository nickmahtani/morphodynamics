# ==============================================================================
# Figure 1b Reproduction: Timnecourse Integration (Day 5-30)
# ==============================================================================
  
setwd("~/Documents/Courses/QuadBio/Morphodynamics")

library(Seurat)
library(dplyr)
library(simspec)
library(MetBrewer)
library(ggplot2)



# ==============================================================================
# 1. Loading individual timepoint RDS files
# ==============================================================================
  

d5 <- readRDS("data/Seurat_object_raw_HB4_D5.rds")
d7 <- readRDS("data/Seurat_object_raw_HB4_D7.rds")
d11 <- readRDS("data/Seurat_object_raw_HB4_D11.rds")
d16 <- readRDS("data/Seurat_object_raw_HB4_D16.rds")
d21 <- readRDS("data/Seurat_object_raw_HB4_D21.rds")
d30 <- readRDS("data/Seurat_object_raw_HB4_D30.rds")

objs <- list(d5, d7, d11, d16, d21, d30)

# ==============================================================================
# 2. Quality control: filtering samples
# ==============================================================================

# to match the QC performed in the paper, I will also nFeature_RNA to be between 
# 1000 to 7,500 and the mitochondrial genes threshold to be less than 10%
  
objs <- lapply(objs, function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT[-\\.]")
  subset(x, subset = nFeature_RNA >1000 & nFeature_RNA <7500 & percent.mt <10)
})

# Verifying that I am on the right track with the paper:
#print(sapply(objs, ncol))
# Expected: D5=5481, D7=8183, D11=4912, D16=6571, D21=7962, D30=7950 

# Visualization with violin plot
#VlnPlot(objs[[1]], features = c("nFeature_RNA", "percent.mt"), ncol = 2)

# ==============================================================================
# 3. Log-normalize each sample
# ==============================================================================

objs <- lapply(objs, NormalizeData)

# ==============================================================================
# 4. Merging all timepoints
# ==============================================================================

timecourse <- merge(objs[[1]], y = objs[2:6])

# ==============================================================================
# 5. Finding variable features
# ==============================================================================

timecourse <- FindVariableFeatures(timecourse, nfeatures = 3000)

# ==============================================================================
# 6. Cell cycle scoring: must happen before scaling since the paper uses cell
# cycle scores for regression
# ==============================================================================

timecourse<-CellCycleScoring(
  timecourse,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident = FALSE
  )

# ==============================================================================
# 7. Scale all genes with cell cycle regression
# ==============================================================================


# Creating genes file with tfs and variable genes as describe in the paper

var_genes <- VariableFeatures(timecourse)
tfs <-readLines("transcription_factors.txt")
cc_genes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
genes_use_raw <- union(var_genes,tfs)
genes_use <- setdiff(genes_use_raw,cc_genes)

timecourse<-ScaleData(
  timecourse, 
  features = genes_use,
  vars.to.regress = c("S.Score", "G2M.Score")
)

# ==============================================================================
# 8. PCA time!
# ==============================================================================

# in my previous code, i had used 50 not 20

timecourse <- RunPCA(timecourse, features = genes_use, npcs = 20)

# ==============================================================================
# 9. Integration with CSS
# ==============================================================================
# i previously used cluster resolution of 0.8

timecourse <- cluster_sim_spectrum(
  timecourse, 
  label_tag = "orig.ident", 
  cluster_resolution = 0.6, 
  dims_use = 1:10)

# ==============================================================================
# 10. Regress cell cycle effects from CSS embeddings
# ==============================================================================

timecourse <- regress_out_from_embeddings(
  timecourse,
  reduction = "css",
  vars_to_regress = c("S.Score", "G2M.Score"),
  reduction.name = "CSSccremoved",
  reduction.key = "CSSccremoved_"
)


# ==============================================================================
# 11. Clustering 
# ==============================================================================
# previously 7.5
timecourse <- FindNeighbors(
  timecourse, 
  reduction = "CSSccremoved", 
  dims = 1:ncol(Embeddings(timecourse, "CSSccremoved"))
  ) %>%
  FindClusters(resolution = 0.6)

# ==============================================================================
# 12. UMAP
# ==============================================================================

timecourse <- RunUMAP(
  timecourse, 
  reduction = "CSSccremoved", 
  dims = 1:ncol(Embeddings(timecourse, "CSSccremoved")),
  spread = 1.0,
  min.dist = 0.4)




# ==============================================================================
# Fancy plotting the UMAP with cluster annotations
# ==============================================================================


Idents(timecourse) <- "seurat_clusters"
umap_clust <- DimPlot(
  timecourse, 
  reduction = "umap",
  label = TRUE, 
  pt.size = 2, 
  label.size = 5, 
  cols = met.brewer("Austria", n =16)
  ) & NoAxes()


# ==============================================================================
# Dotplot analysis to help with cluster labelling
# ==============================================================================


Idents(timecourse) <- "seurat_clusters"

genes_set <- c( "POU5F1","KRT8","MKI67","SFRP2","SOX2","TCF7L2","PAX6","NES",
                "VIM","HES4","HES5","STMN2","GLI3","CDH2","SOX21","OTX2","LHX2",
                "FEZF1","FEZF2","SOX1","FOXG1","SIX3","RAX","SIX6","PAX2",
                "ZIC5","FOXO1","FOXO3","CDKN1A","HESX1","ZIC1","PRTG","BMP7",
                "FGF8","NRG1","NR2F1","NR2F2","PTX3","WLS","EMX1","EMX2","DLX2",
                "NKX2-1","TAL2","RSPO3","LHX5","LHX9","WNT4","WNT2B","WNT8B",
                "SST","PROM1","ITGA6","ITGB1","BTG2","OCLN",
                "CYP26A1","SOX10","FOXD3")

dotplot <- DotPlot(
  object = timecourse,
  features = genes_set,
  cluster.idents = T,
  dot.scale = 6,
  dot.min = 0,
  scale.by = "radius"
)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_colour_gradient2(low="darkgreen", mid="grey90", high="purple4")+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2)
  
dot_data <- dotplot$data


# to add the labels, I filter for high pct.exp and avg.exp, grouped by cluster
dot_data <- dot_data %>%
  group_by(id) %>%  # id is the cluster column in DotPlot data
  arrange(id, desc(pct.exp), desc(avg.exp.scaled)) %>%
  filter(pct.exp > 20)%>%  
  slice_head(n=5)


# ==============================================================================
# Average expression of manually selected marker genes to also help with labeling
# ==============================================================================


Idents(timecourse) <- "seurat_clusters"

# Manually curated list of marker genes based on the paper 

marker_sets <- list(
  Neurectoderm           = c("POU5F1", "KRT8", "SFRP2", "PROM1"),
  Late_Neurectoderm      = c("SOX2", "NES", "VIM", "HES4", "HES5", "PAX6"),
  Prosencephalic_prog    = c("OTX2", "SOX21", "HESX1", "ZIC1", "CDH2"),
  Late_Prosencephalic    = c("TCF7L2", "GLI3", "PRTG", "BMP7", "FGF8"),
  Telencephalic_prog     = c("FOXG1", "FEZF1", "FEZF2", "EMX1", "EMX2", "SIX3", "LHX2"),
  Diencephalic_prog      = c("RAX", "SIX6", "TAL2", "RSPO3", "WNT4", "NRG1", "NR2F1"),
  TelDie_neurons         = c("STMN2", "DLX2", "SST"),
  Unknown_prolif         = c("MKI67", "POU5F1")
)

# making a clean gene list from the named list above
all_markers <-unique(unlist(marker_sets))

# getting normalized expression for all marker genes
expr_data <- FetchData(timecourse, vars = all_markers, layer = "data")
expr_data$cluster <- timecourse$seurat_clusters

# Computing the mean per cluster
avg_expr <- expr_data %>%
  group_by(cluster) %>%
  summarise(across(all_of(all_markers), mean)) %>%
  tibble::column_to_rownames("cluster")%>%
  t()

# printing scores for each cell type
for (ct in names(marker_sets)) {
  genes <- marker_sets[[ct]]
  scores <- colMeans(avg_expr[genes, ])
  cat(paste0(ct, ":\n"))
  cat("  Top clusters: ")
  top <- sort(scores, decreasing = TRUE)[1:3]
  cat(paste(paste0(names(top), " (", round(top, 2), ")"), collapse = ", "), "\n")
}

# ==============================================================================
# Timepoint composition per cluster to also help with labeling
# ==============================================================================

# Printing timepoint composition per cluster to help
comp <- table(timecourse$seurat_clusters, timecourse$orig.ident)
comp_pct <- round(prop.table(comp, margin = 1) * 100, 1)
print(comp_pct)


# ==============================================================================
# Assignment of new cluster ids based on my judgment from the  above data
# ==============================================================================

# saving the current cluster ids as backup
timecourse[["cluster.ids"]] <- Idents(object = timecourse)

# enumerate current cluster IDs and the labels for them
cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,12,13)

class <-  c("Neurectoderm", #  0
             "Telencephalic progenitors",#  1
             "Neurectoderm", #  2
             "Telencephalic progenitors", #3
             "Late Prosencephalic progenitors",#  4
             "Diencephalic progenitors", #  5
             "Neurectoderm", #  6
             "Telencephalic progenitors",# 7
             "Late Neurectoderm",# 8
             "Late Neurectoderm",# 9
             "Late Neurectoderm",#   10
             "Tel/Die neurons",#  11
             "Prosencephalic progenitors", # 12
             "Prosencephalic progenitors")#     13

timecourse@meta.data[,'class'] <- plyr::mapvalues(x = timecourse@meta.data$cluster.ids, from = cluster.ids, to = class)

timecourse$class <- factor(
  timecourse$class, 
  levels = c("Neurectoderm","Late Neurectoderm","Prosencephalic progenitors",
             "Telencephalic progenitors","Late Prosencephalic progenitors",
             "Diencephalic progenitors","Tel/Die neurons"))

color_anno = c("#ff6e1e","#e6c262","#7fb59e","#049983","#016260","#264653","#930707","#14231d","#0c3447")

classdimplot <- DimPlot(
  timecourse, 
  reduction = "umap", 
  pt.size = 0.8, label = F,
  order = T, 
  group.by = 'class', 
  cols = color_anno) & NoAxes()

# ==============================================================================
# UMAP grouping by timepoint - Figure 2c
# ==============================================================================

timecourse$timepoint <- factor(
  timecourse$orig.ident,
  levels = c("HB4_D5", "HB4_D7", "HB4_D11", "HB4_D16", "HB4_D21", "HB4_D30"),
  labels = c("5", "7", "11", "16", "21", "30")
)

timepoint_colors <- c(
  "5"  = "#3d0c02",  
  "7"  = "#b94a1a",  
  "11" = "#e8a830", 
  "16" = "#b8cce0",  
  "21" = "#2d6b4a",  
  "30" = "#1a2f4a"   
)

timepointdimplot <- DimPlot(
  timecourse, 
  reduction = "umap", 
  pt.size = 0.3, label = F,
  order = T, 
  group.by = 'timepoint', 
  cols = timepoint_colors) +
  labs(title = "Timepoint", color = "Time (days)")+
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 14),  # left align title
    legend.key.size = unit(1.2, "cm")  # bigger legend dots
  ) &
  NoAxes()



# ==============================================================================
# Bar plot
# ==============================================================================


# what is the cell identity distribution in each of the days




# calculate cluster proportions
bar_data <- timecourse@meta.data %>%
  group_by(timepoint, class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(proportion = n / sum(n))

# plot
barplot <- ggplot(bar_data, aes(x = timepoint, y = proportion, fill = class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_anno) +
  labs(
    x = "Time (day)",
    y = "Cluster proportion",
    fill = NULL  # removes legend title
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )














# Archive




# to see the columns in the metadata
#colnames(Timecourse[[]])


#plot1 <- UMAPPlot(timecourse, group.by="orig.ident")
plot2 <- UMAPPlot(timecourse, label = T)

#if needed
# Flip horizontal (mirror left-right)
#timecourse@reductions$umap@cell.embeddings[,1] <- -timecourse@reductions$umap@cell.embeddings[,1]                                                                                                                                                                                                                                                                                                                                                                     

# Flip vertical (mirror up-down)                                                                                                                                                                                                                                                                                                                                                                                                                                      
#timecourse@reductions$umap@cell.embeddings[,2] <- -timecourse@reductions$umap@cell.embeddings[,2]



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




# getting the list of TF from AnimalTFDB

tfdb <- read.delim("Homo_sapiens_TF.txt", header = TRUE)
tfs <- tfdb$Symbol
# Remove blanks
tfs <- tfs[tfs!=""]

#saving into a new file
writeLines(tfs,"transcription_factors.txt")



