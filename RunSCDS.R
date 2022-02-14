library(Seurat)
library(scds)
library(dplyr)
library(tictoc)

tic("scds doublets")
(seurat %>% 
 as.SingleCellExperiment %>% 
 cxds_bcds_hybrid)@colData[,c('cxds_score',
                              'bcds_score',
                              'hybrid_score')
                          ]%>% 
as.data.frame -> scds.doublet.profiles
toc()


meta             <- merge(seurat@meta.data, scds.doublet.profiles, by.x=0, by.y=0)
rownames(meta)   <- meta$Row.names
meta$Row.names   <- NULL
seurat@meta.data <- meta
