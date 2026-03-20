library(Seurat)
library(reticulate)
library(anndata)

data <- read_h5ad("~/OneDrive - University of Eastern Finland/Projects/Minna_Spatial_Transcriptomics_Data/Updated_mouse_Aortas_data.h5ad")
data <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs,min.features = 500, min.cells = 30)
saveRDS(data,"~/OneDrive - University of Eastern Finland/Projects/Minna_Spatial_Transcriptomics_Data/Mouse_Aortas.rds")
