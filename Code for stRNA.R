library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)
library(harmony)
library(DoubletFinder)
library(glmGamPoi)
library(Rcpp)
library(openxlsx)
library(tibble)
library(here)
library(phylogram)
library(ggrastr)
library(ggbeeswarm)
library(beeswarm)
library(vipor)
library(Cairo)
library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(gt)
library(RColorBrewer)
library(gridExtra)

Brain_ST = Load10X_Spatial(data.dir = "./cHC-1L",
                           filename = "filtered_feature_bc_matrix.h5",
                           assay = "Spatial")

Brain_ST <- SCTransform(Brain_ST, assay = "Spatial", verbose = FALSE)
Brain_ST <- RunPCA(Brain_ST, assay = "SCT", verbose = FALSE)
Brain_ST <- FindNeighbors(Brain_ST, reduction = "pca", dims = 1:30)
Brain_ST <- FindClusters(Brain_ST, verbose = FALSE,resolution = 0.2)
Brain_ST <- RunUMAP(Brain_ST, reduction = "pca", dims = 1:30)

p1 <- DimPlot(Brain_ST, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(Brain_ST, label = TRUE, label.size = 3)
p1 + p2

#load scRNA
Brain_scRNA <- CreateSeuratObject(counts = Read10X("./P01T/filtered_feature_bc_matrix/"), project = "scRNA" )
Brain_scRNA[["percent.mt"]] <- PercentageFeatureSet(Brain_scRNA, pattern = "^MT-")
Brain_scRNA <- subset(Brain_scRNA, subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 6000 & 
                        percent.mt < 30)

Brain_scRNA <- NormalizeData(Brain_scRNA)
Brain_scRNA <- FindVariableFeatures(Brain_scRNA, selection.method = "vst", nfeatures = 2000)
Brain_scRNA <- ScaleData(Brain_scRNA)
Brain_scRNA <- RunPCA(Brain_scRNA, npcs = 50)
Brain_scRNA <- FindNeighbors(Brain_scRNA, dims = 1:30)
Brain_scRNA <- FindClusters(Brain_scRNA, resolution = 0.5)
Brain_scRNA <- RunUMAP(Brain_scRNA, dims = 1:30)

Idents(Brain_scRNA) <- "seurat_clusters"
p11 <- DimPlot(Brain_scRNA, reduction = "umap", pt.size=0.5, label = F,repel = TRUE)

Brain_scRNA.markers <- FindAllMarkers(Brain_scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Brain_scRNA.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(data.frame(ID=rownames(Brain_scRNA.markers),Brain_scRNA.markers), file = "cluster_marker_genes.txt", sep = "\t", col.names = T, row.names = F, quote = F)

#Annotation of cell subpopulations based on differentially expressed genes
Brain_scRNA <- RenameIdents(Brain_scRNA,"0"="CD4",
                            "1"="CD8", 
                            "2"="Neo", 
                            "3"= "CD8_T", 
                            "4"= "CD14+_Monocytes", 
                            "5"= "NK_Cells",
                            "6"= "Ductal_Cells",
                            "7"= "cDC2",
                            "8"= "CD56+ NK_Cells",
                            "9"= "Treg",
                            "10"= "B",
                            "11"= "Endothelial_Cells",
                            "12"= "Pan_CD8",
                            "13"= "Monocytes",
                            "14"= "EpithelialCells",
                            "15"= "Smooth_Muscle_cell",
                            "16"= "Mast_Cells",
                            "17"= "cDC1",
                            "18"= "Plasmacytoid_Dendritic_Cells",
                            "19"= "SmoothMuscle"
)

Brain_scRNA@meta.data$celltype <- Idents(Brain_scRNA)
p22 <- DimPlot(Brain_scRNA, reduction = 'umap', 
               label = TRUE, pt.size = 0.5) + NoLegend()
p11 + p22

Brain_scRNA_reference <- SCTransform(Brain_scRNA, ncells = 7181, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
DimPlot(Brain_scRNA_reference, group.by = "Subtype", label = TRUE)
Brain_ST_2 <- Brain_ST

###Annotate stRNA data based on scRNA cellular subtype information
anchors <- FindTransferAnchors(reference = Brain_scRNA_reference, 
                               query = Brain_ST_2, 
                               normalization.method = "SCT")

predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = Brain_scRNA_reference$Subtype, 
                                  prediction.assay = TRUE,
                                  weight.reduction = Brain_ST_2[["pca"]], 
                                  k.weight = 40,
                                  dims = 1:30)

class(predictions.assay)
predictions.assay@data[,1:4]

predictions.res <- predictions.assay@data
rownames(predictions.res) <- make.names(rownames(predictions.res))
cell_types <- rownames(predictions.res)
for(cell_type in cell_types){
  Brain_ST_2 <- AddMetaData(object= Brain_ST_2,metadata = predictions.res[cell_type,],col.name = cell_type)
}

Brain_ST_2[["predictions"]] <- predictions.assay

DefaultAssay(Brain_ST_2) <- "predictions"
SpatialFeaturePlot(Brain_ST_2, features = c("Ductal.Cells","Pan.CD8"),
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE)
























