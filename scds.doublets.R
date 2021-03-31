sce <- SingleCellExperiment(list(counts=as.matrix(seu@assays$RNA@counts)))
colData(sce)<-DataFrame(seu@meta.data)
logcounts(sce) = log1p(counts(sce))
vrs            = apply(logcounts(sce),1,var)
pc             = rsvd::rpca(t(logcounts(sce)[order(vrs,decreasing=TRUE)[1:100],]))
reducedDim(sce, "umap")<-seu@reductions$umap@cell.embeddings
rm(vrs,pc)

sce = cxds(sce,retRes = TRUE)
sce = bcds(sce,retRes = TRUE,verb=TRUE)
sce = cxds_bcds_hybrid(sce)

plotReducedDim(sce, "umap",colour_by = "cxds_score")
plotReducedDim(sce, "umap",colour_by = "bcds_score")
plotReducedDim(sce, "umap",colour_by = "hybrid_score")


scds.meta <- data.frame(colData(sce)[,c("cxds_score","bcds_score","hybrid_score")])
new.seu.meta <- cbind(seu@meta.data, scds.meta[rownames(seu@meta.data),] )
seu@meta.data<-new.seu.meta

FeaturePlot(seu, features="cxds_score", order=T)&theme(legend.position=c(0.1,0.2))&ggtitle("SCDS_cxds_score")
FeaturePlot(seu, features="bcds_score", order=T)&theme(legend.position=c(0.1,0.2))&ggtitle("SCDS_bcds_score")
FeaturePlot(seu, features="hybrid_score", order=T)&theme(legend.position=c(0.1,0.2))&ggtitle("SCDS_hybrid_score")
