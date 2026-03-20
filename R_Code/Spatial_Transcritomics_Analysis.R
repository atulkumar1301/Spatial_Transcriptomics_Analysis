library(Seurat)
library (future)
plan ("multisession", workers = 10)
library(ggplot2)
library(arrow)
library(sctransform)
library(spacexr)


path <- "~/OneDrive - University of Eastern Finland/Projects/Minna_Spatial_Transcriptomics_Data/Control_Aortas_01__20250416__172150/Xenium_Ranger/baysor-demo/outs/"

# Load the Xenium data
xenium.obj <- LoadXenium(path, fov = "fov", segmentations = "cell")

# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

xenium.obj$nCount_Xenium

xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")

xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))

xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)

xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)

xenium.obj <- FindClusters(xenium.obj, resolution = 0.3)

DimPlot(xenium.obj)

FeaturePlot(xenium.obj, features = c("Myh11", "Acta2", "Tagln", "Col1a1", "Dcn", "Lum", "Cdh5", "Pecam1", "Vwf", "Cd68", "Adgre1", "Gpnmb"))

ImageDimPlot(xenium.obj, cols = "polychrome", size = 0.75)


query.counts <- GetAssayData(xenium.obj, assay = "Xenium", layer = "counts")[, Cells(xenium.obj)]

coords <- GetTissueCoordinates(xenium.obj, which = "centroids")

rownames(coords) <- coords$cell

coords$cell <- NULL

query <- SpatialRNA(coords, query.counts, colSums(query.counts))


mouse.aortas.ref <- readRDS("~/OneDrive - University of Eastern Finland/Projects/Minna_Spatial_Transcriptomics_Data/Mouse_Aortas.rds")
mouse.aortas.ref <- UpdateSeuratObject(mouse.aortas.ref)

Idents(mouse.aortas.ref) <- "cell_type"
# remove CR cells because there aren't enough of them for annotation

#mouse.aortas.ref <- subset(mouse.aortas.ref, subset = cell_type != "epithelial cell")

counts <- GetAssayData(mouse.aortas.ref, assay = "RNA", layer = "counts")

cluster <- as.factor(mouse.aortas.ref$cell_type)

names(cluster) <- colnames(mouse.aortas.ref)

nUMI <- mouse.aortas.ref$nCount_RNA

names(nUMI) <- colnames(mouse.aortas.ref)

nUMI <- colSums(round(counts))

levels(cluster) <- gsub("/", "-", levels(cluster))

reference <- Reference(round(counts), cluster, nUMI)

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8, CELL_MIN_INSTANCE = 15)

RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

annotations.df <- RCTD@results$results_df

annotations <- annotations.df$first_type

names(annotations) <- rownames(annotations.df)

xenium.obj$predicted.celltype <- annotations

keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]

xenium.obj <- subset(xenium.obj, cells = keep.cells)

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "predicted.celltype",
                              niches.k = 5, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype", size = 1.5, cols = "polychrome",
                              dark.background = F) + ggtitle("Cell type")

niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") +
  scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))

celltype.plot | niche.plot
