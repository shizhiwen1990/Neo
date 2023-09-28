library(Seurat)
library(dplyr)
library(sctransform)
library(DoubletFinder)

###1. Data loading and QC
#sample list
sample_id <- read_excel("./patients_metadata.xlsx", range = cell_cols("A:A")) %>% .$sample_id

#import cellranger files from different data sets
for (i in seq_along(sample_id)){
  assign(paste0("scs_data", i), Read10X(data.dir = paste0("./cellranger/", sample_id[i], "/filtered_feature_bc_matrix")))
}

#create seurat objects from cellranger files 
for (i in seq_along(sample_id)){
  assign(paste0("seu_obj", i), CreateSeuratObject(counts = eval(parse(text = paste0("scs_data", i))), project = sample_id[i], min.cells = 3))
}

### merge data sets
seu_obj <- merge(seu_obj1, y = c(seu_obj2, seu_obj3, seu_obj4, seu_obj5, seu_obj6, seu_obj7, seu_obj8, seu_obj9, seu_obj10), add.cell.ids = sample_id, project = "PACA")

### calculate mitochondrial, hemoglobin and ribosomal gene counts
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^MT-", col.name = "pMT")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^HBA|^HBB", col.name = "pHB")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^RPS|^RPL", col.name = "pRP")
nFeature_lower <- 200
pMT_upper <- 25
qcparams <- c("nFeature_RNA", "pMT")
seu_obj <- subset(seu_obj, subset = nFeature_RNA > nFeature_lower & pMT < pMT_upper)

###Doublet removal with the assumption that doublets represent 6% of cells
FindDoublets <- function(sample_id, seurat_aggregate) {
  seu_obj <- seurat_aggregate
  seurat_obj <- subset(seu_obj, idents = sample_id) #metadata$sample_id
  View(seurat_obj@meta.data)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- RunPCA(seurat_obj)
  # ElbowPlot(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:9)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:9)
  DimPlot(seurat_obj)
  #ggsave2("SuppFig2B.png", path = "./results", width = 30, height = 30, units = "cm")
  
  ## pK Identification (no ground-truth) Ѱ??????PKֵ---------------------------------------------------------------------------------------
  sweep.res.list_lung <- paramSweep_v3(seurat_obj, PCs = 1:9, sct = F)
  sweep.stats_lung <- summarizeSweep(sweep.res.list_lung, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_lung)
  pK <- bcmvn_kidney %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  seurat_doublets <- doubletFinder_v3(seurat_obj, PCs = 1:9, pN = 0.25, pK = pK,
                                      nExp = round(0.05*length(seurat_obj@active.ident)), 
                                      reuse.pANN = FALSE, sct = F)
  
  #write.xlsx(seurat_doublets@meta.data, "./doublets_metadata.xlsx", sheetName="Sheet1",col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE) #????metadata??Ϣ??????ÿ??ϸ????????��Դ???ٴ???Ϣ??ϸ??ע????Ϣ
  
  #View(seurat_doublets@meta.data)
  
  # create doublet groupings and visualize results
  DF.class <- names(seurat_doublets@meta.data) %>% str_subset("DF.classifications")
  pANN <- names(seurat_doublets@meta.data) %>% str_subset("pANN")
  p1 <- ggplot(bcmvn_kidney, aes(x=pK, y=BCmetric)) +
    geom_bar(stat = "identity") + 
    ggtitle(paste0("pKmax=",pK)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) #PKֵ??ͼ
  p2 <- DimPlot(seurat_doublets,reduction = "umap", group.by = DF.class) #Doubletչʾ
  p3 <- FeaturePlot(seurat_doublets, features = pANN)
  outFile <- paste0("doublets.",sample_id,".pdf")
  dir.create("plots", showWarnings = FALSE)
  pdf(here("plots",outFile))
  print(p1) # need to use print() when drawing pdf in a function call
  print(p2)
  print(p3)
  dev.off()
  # create a df of barcodes and doublet designations
  df_doublet_barcodes <- as.data.frame(cbind(rownames(seurat_doublets@meta.data), seurat_doublets@meta.data[[DF.class]]))
  return(df_doublet_barcodes)
}
# filter out doublets prior to snRNA preprocessing
Idents(seu_obj) <- "orig.ident"
list.doublet.bc <- lapply(sample_id, function(x) {FindDoublets(x, seurat_aggregate = seu_obj)}) #?ò?????ÿ?????????��???doublets????
#??ȡ????ϸ??doublets??????????????ϸ??or˫ϸ??
doublet_id <- list.doublet.bc %>%
  bind_rows() %>%
  dplyr::rename("doublet_id" = "V2") %>%
  tibble::column_to_rownames(var = "V1") # this is the barcode column
table(doublet_id) # quantify total doublet vs. singlet calls (expect ~5% doublets)
seu_obj <- AddMetaData(seu_obj,doublet_id)
Idents(seu_obj) <- "doublet_id"
seu_obj <- subset(seu_obj,idents = "Singlet") 

#Export input data for machine learning model
seu_obj <- SCTransform(seu_obj, vars.to.regress = c("pMT","nCount_RNA"), verbose = TRUE) 
SCT_matrix <- seu_obj@assays$SCT@scale.data
input_data <- scale(SCT_matrix)
write.table(data.frame(ID=rownames(input_data),input_data), file = "input_data.txt", sep = "\t", col.names = T, row.names = F, quote = F)


















