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
library(data.table)

cdsAB.anchors <- FindIntegrationAnchors(object.list = list(P01T,P02T,P03T,P04T,P05T,P06T,P07T,T01,T02), dims = 1:20)
cdsAB.integrated <- IntegrateData(anchorset = cdsAB.anchors, dims = 1:20)
DefaultAssay(cdsAB.integrated) <- "integrated"
cdsAB.integrated <- ScaleData(cdsAB.integrated, features = rownames(cdsAB.integrated))
cdsAB.integrated <- RunPCA(cdsAB.integrated, npcs = 50, verbose = FALSE)
ElbowPlot(cdsAB.integrated, ndims = 50)
cdsAB.integrated <- FindNeighbors(cdsAB.integrated, dims = 1:20)
saveRDS(cdsAB.integrated,"./CAA.integrated")
cdsAB.integrated <- RunUMAP(cdsAB.integrated, dims = 1:20)

for (i in c(0.3,0.4,0.5,1)) {
  cdsAB.integrated <- FindClusters(cdsAB.integrated, resolution = i)
  print(DimPlot(cdsAB.integrated, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}
cdsAB.integrated <- FindClusters(cdsAB.integrated, resolution = 0.5)

p1 <- DimPlot(cdsAB.integrated, reduction = "umap", label = TRUE)
p2 = DimPlot(cdsAB.integrated, reduction = "umap", group.by='orig.ident')
plotc <- p1+p2
plotc

seu_obj <- cdsAB.integrated
seu_obj.markers <- FindAllMarkers(seu_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seu_obj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(data.frame(ID=rownames(seu_obj.markers),seu_obj.markers), file = "cluster_marker_genes.txt", sep = "\t", col.names = T, row.names = F, quote = F)

top10 = seu_obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seu_obj,top10$gene,size = 3)+scale_fill_gradientn(colors = c("darkblue", "white", "red")) #??ÿ????Ⱥtop5???????Ƶ???ͼ
p = DotPlot(seu_obj,features = unique(top10$gene),assay = "RNA")
p

saveRDS(seu_obj,"./CAA_umap.rds")

###Annotation with ScType algorithm
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

#Cell type assignment
source("./annotation/gene_sets_prepare.R")
# load cell type annotation function
source("./annotation/sctype_score_.R")

# DB file
db_ = "./annotation/ScTypeDB_full.xlsx"
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

#Finally, let's assign cell types to each cluster
es.max = sctype_score(scRNAseqData = seu_obj[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(seu_obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seu_obj@meta.data[seu_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu_obj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown" ????��ע????Ϣ????Ϊunknown
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

seu_obj@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seu_obj@meta.data$customclassif[seu_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

View(seu_obj@meta.data)
p1 = DimPlot(seu_obj, reduction = "umap", group.by='customclassif')
p2 = DimPlot(seu_obj, reduction = "umap", group.by='orig.ident')
p3 = DimPlot(seu_obj, reduction = "umap", group.by='seurat_clusters')
plotc <- p1+p2+p3

###Add neo information
metatable <- read_excel("./metadata/Neo_metadata.xlsx")
metadata <- FetchData(seu_obj, "cell.name")
metadata$cell_id <- rownames(metadata)
metadata <- left_join(x = metadata, y = metatable, by = "cell.name")
rownames(metadata) <- metadata$cell_id
seu_obj <- AddMetaData(seu_obj, metadata = metadata)

Idents(seu_obj) <- seu_obj@meta.data$Specific
DimPlot(seu_obj, cells.highlight=WhichCells(seu_obj, idents="Neo"),
        cols.highlight = "red" ) 

DimPlot(seu_obj, cells.highlight=WhichCells(seu_obj, idents="Pre-Neo"),
        cols.highlight = "blue",
        sizes.highlight = 1) 

DimPlot(seu_obj,
        cells.highlight=list(
          c1=WhichCells(seu_obj, idents="Pre-Neo"),
          c2=WhichCells(seu_obj, idents="Neo")
        ),
        cols.highlight=c("red", "blue") 
        
)

seu_obj <- readRDS("./Neo_annotation.rds")

#monocle
library(tidyverse)
library(ggpubr)
library(cowplot)
library(dplyr)
library(ggplot2)
library(tidyr)
library(piano)
library(msigdbr)
library(monocle)
library(RColorBrewer)

seu_obj$customclassif <- as.character(seu_obj$customclassif) 
seu_obj$customclassif_2 <- as.character(seu_obj$customclassif) 

seu_obj$customclassif[seu_obj$customclassif_2 == "Trm"] <- "Trm" 
seu_obj$customclassif[seu_obj$customclassif_2 == "Tem"] <- "Tem"
seu_obj$customclassif[seu_obj$customclassif_2 == "Exhausted"] <- "Exhausted"
seu_obj$customclassif[seu_obj$customclassif_2 == "NME1+"] <- "NME1+"
seu_obj$customclassif[seu_obj$customclassif_2 == "Tinn"] <- "Tinn"
seu_obj$customclassif[seu_obj$customclassif_2 == "HSP+"] <- "HSP+"
seu_obj$customclassif[seu_obj$customclassif_2 == "Temra"] <- "Temra"

SO.traj1 <- seu_obj[,seu_obj$customclassif %in% c("Trm", "Tem", "Exhausted", "NME1+", "Tinn", "HSP+", "Temra")] 
View(SO.traj1@meta.data)

pd <- new('AnnotatedDataFrame', data = SO.traj1@meta.data)  
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(SO.traj1), row.names = row.names(SO.traj1)))
cds <- newCellDataSet(as(SO.traj1@assays$integrated@data, "matrix"), phenoData = pd, featureData = fd, expressionFamily = negbinomial.size()) 
cds$clusters <- SO.traj1$customclassif

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds) 

markers <- FindAllMarkers(object = SO.traj1, min.pct = 0.25, thresh.use = 0.25)
markers <- subset(markers, p_val_adj < 0.05)
order.genes <- unique(as.character(markers$gene))

cds <- setOrderingFilter(cds, order.genes)
cds <- reduceDimension(cds = cds, max_components = 3,method = 'DDRTree')  
cds <- orderCells(cds)

saveRDS(cds,"./cds.rds")
cds <- readRDS("./cds.rds")

plot_cell_trajectory(cds,color_by = "customclassif")
ggsave("celltype.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
plot_cell_trajectory(cds,color_by = "State")
ggsave("State.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
plot_cell_trajectory(cds,color_by = "Pseudotime")
ggsave("Pseudotime.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
plot_cell_trajectory(cds,color_by = "customclassif")+facet_wrap(~customclassif,nrow=1)
ggsave("celltypeb.pdf",device = "pdf",width = 21,height = 9,units = c("cm"))

###Plot
seu_obj <- readRDS("./CAA_Neo_annotation.rds")
library(ClusterGVis)
library(org.Hs.eg.db)
library(Seurat)
library(dplyr)

levels(seu_obj)
Idents(seu_obj) <- seu_obj@meta.data$customclassif

seu_obj.markers.all <- Seurat::FindAllMarkers(seu_obj,
                                              only.pos = TRUE,
                                              min.pct = 0.25,
                                              logfc.threshold = 0.25)
seu_obj.markers <- seu_obj.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

head(seu_obj.markers)

st.data <- prepareDataFromscRNA(object = seu_obj,
                                diffData = seu_obj.markers,
                                showAverage = TRUE)

str(st.data)

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.25,
                        topn = 5,
                        seed = 5201314)

head(enrich)

markGenes = unique(seu_obj.markers$gene)[sample(1:length(unique(seu_obj.markers$gene)),60,
                                                replace = F)]

visCluster(object = st.data,
           plot.type = "line")

pdf('sc1.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:9)) #cluster.order
dev.off()

pdf('sc2.pdf',height = 10,width = 14,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 40,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:9),
           go.col = rep(jjAnno::useMyCol("stallion",n = 8),each = 5),
           add.bar = T)
dev.off()

st.data <- prepareDataFromseu_obj(object = seu_obj,
                                  diffData = seu_obj.markers,
                                  showAverage = TRUE,
                                  keep.uniqGene = FALSE,
                                  sep = "_")

df <- st.data$wide.res

visCluster(object = st.data,
           plot.type = "line")

pdf('sc3.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:9))
dev.off()

seu_obj.markers1 <- seu_obj.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 6, wt = avg_log2FC)

st.data <- prepareDataFromseu_obj(object = seu_obj,
                                  diffData = seu_obj.markers1,
                                  showAverage = FALSE)
pdf('sc4.pdf',height = 10,width = 8,onefile = F)

visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           markGenes = unique(seu_obj.markers1$gene),
           markGenes.side = "left",
           annoTerm.data = enrich,
           #show_column_names = F,
           line.side = "left",
           cluster.order = c(1:9),
           add.bar = T)
dev.off()



















