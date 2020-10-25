seu <- ASI.ILC1 

tic("Rescaling and Dimensionality reduction")
DefaultAssay(seu) <- "RNA"
seu  <- FindVariableFeatures(seu,  selection.method = "vst", nfeatures = 2000)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 3*1000 * 1024^2)
seu <- ScaleData(seu,  features=rownames(seu), block.size = 500, min.cells.to.block = 1000)
plan("sequential")

seu <- RunPCA(seu,  npcs=50, features=VariableFeatures(object=seu),  verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- RunUMAP(object=seu,  dims=1:20, verbose=FALSE, min.dist=0.5, n.neighbors = 30L) # uwot in seurat3.1 is problematic
#seu <- RunUMAP(object=seu,  dims=1:30, verbose=FALSE, min.dist=0.5, n.neighbors = 30L, umap.method = "umap-learn", metric = "correlation")

seu <- FindClusters(seu, algorithm=4, resolution=0.5, group.singletons=TRUE, verbose=FALSE)
toc()

seu -> ASI.ILC1
