AddAugmentedAssay<- function(seu, 
            meta.data.features=c("bcds_score", "hybrid_score",
                                 "cxds_score", "DF_pANN","percent.mt","percent.disso") ){
    dge<-as.matrix(GetAssayData(object = seu, slot = "data",assay = "RNA"))
    meta<-t(as.matrix(seu@meta.data[colnames(dge), meta.data.features]))
    seu[["AugRNA"]]<-CreateAssayObject(data=rbind(dge, meta))
    DefaultAssay(seu)<-"AugRNA"
    
    hvg<-c(VariableFeatures(seu,assay="RNA"),meta.data.features)
    VariableFeatures(seu,assay="AugRNA")<-hvg
    
    seu <- ScaleData(seu, features =VariableFeatures(seu, assay="AugRNA")  )
    return(seu)
}

CentrifugeData <- function(sub,meta.data.features){
    tic("generating new assay")
    sub <- NormalizeData(sub, normalization.method="LogNormalize", scale.factor=10000)
    sub <- FindVariableFeatures(sub,  selection.method = "vst", nfeatures = 2000)
    sub <- ScaleData(sub, features=VariableFeatures(sub))
    sub <- AddAugmentedAssay(sub,meta.data.features)

    DefaultAssay(sub)<-"AugRNA"
    sub <- RunPCA(sub,  npcs=30, features=VariableFeatures(sub,assay="AugRNA"),verbose=FALSE)
    sub <- RunUMAP(object=sub, dims=1:30, verbose=FALSE)
    sub <- FindNeighbors(sub, dims = 1:30,graph.name ="sub")
    sub <- FindClusters(sub, resolution = c(0.5,0.7,0.9,1.2), algorithm=4, group.singletons=T,graph.name="sub")
    toc()
    return (sub)
}
