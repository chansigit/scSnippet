FindSubClusters<- function(seu, idents, clusters.to.refine, new.cluster.name, resolution){
    RenormalizeData<-function(seu, resolution, npcs=20){
        DefaultAssay(seu) <- "RNA"
        seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 1000,verbose = F)
        seu <- ScaleData(seu, features=VariableFeatures(object=seu),verbose = F)
        seu <- RunPCA(seu, npcs=npcs, features=VariableFeatures(object=seu),verbose=F)
        seu <- FindNeighbors(seu, dims=1:npcs, graph.name ="tmp",verbose = F)
        seu <- FindClusters(seu, algorithm=4, resolution=resolution, group.singletons=T, verbose=F, graph.name = "tmp")
        return (seu)
    }
    Idents(seu)<- idents
    seu@meta.data[ ,new.cluster.name] <-as.character(seu@meta.data[ ,idents] )
    message(paste("working on clustername=",idents))
    for (cl in clusters.to.refine){
        message(paste("processing cluster",cl))
        sub <- subset(seu, idents = cl)
        sub <- RenormalizeData(sub, resolution=1.0, npcs=20)
        tmp_cluster <- paste("tmp_res.",resolution,sep="")
        sub@meta.data[,tmp_cluster] <- paste(cl, sub@meta.data[,tmp_cluster] ,sep=".")
        
        metachanged<- seu@meta.data
        cellchanged<- rownames(sub@meta.data)
        
        metachanged[cellchanged, new.cluster.name]<-sub@meta.data[cellchanged,tmp_cluster]
        
        metachanged-> seu@meta.data
       
    }
    return(seu)
}

FindSubClusters(T.cells, idents="sub_res.1.5", clusters.to.refine=c(9,11), 
                new.cluster.name="sub_res.1.5_ref1.0",resolution=1.0)[[]]
