library(purrr)
GatherData <- function(object, ...) {
  UseMethod("GatherData")
}


GatherData.Seurat <- function(object,assay,slot_use,...) {
    assay    <-assay %||% "RNA"
    slot_use <-slot_use %||% "data"
    obj_data <-GetAssayData(object=object, assay=assay, slot=slot_use) %>%
               as.matrix()
    return(obj_data)
}


PseudoCell <- function(object,assay_use = NULL,slot_use = NULL,
                        cluster_use =NULL, pseudocell.size  =NULL){
    message("tips: Cluster_use : one col in metadata\npseudocell.size : how many cell will be pseudo")
            
    Inter<- GatherData(object = object,
                       assay = assay_use,
                       slot_use = slot_use) 
    Inter[Inter<0]=0
    idd<-object@meta.data
    Inter.id<-cbind(rownames(idd), as.vector(idd[,cluster_use]))

    rownames(Inter.id)<-rownames(idd)
    colnames(Inter.id)<-c("CellID","Celltype")
    Inter.id<-as.data.frame(Inter.id)
    Inter.id$Celltype <-as.factor(Inter.id$Celltype)
    Inter1<-Inter[,Inter.id$CellID]
    Inter<-as.matrix(Inter1)
    pseudocell.size = pseudocell.size ## 10 test
    new_ids_list = list()
    for (i in 1:length(levels(Inter.id$Celltype))) {
        cluster_id = levels(Inter.id$Celltype)[i]
        cluster_cells <- rownames(Inter.id[Inter.id$Celltype == cluster_id,])
        cluster_size <- length(cluster_cells)       
        pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
        pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
        names(pseudo_ids) <- sample(cluster_cells)  
        new_ids_list[[i]] <- pseudo_ids     
    }

    new_ids <- unlist(new_ids_list)
    new_ids <- as.data.frame(new_ids)
    new_ids_length <- table(new_ids)

    new_colnames <- rownames(new_ids)  ###add
    all.data<-Inter[,as.character(new_colnames)] ###add
    all.data <- t(all.data)###add

    new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
                        list(name=new_ids[,1]),FUN=mean)
    rownames(new.data)<-new.data$name
    new.data<-new.data[,-1]

    new_ids_length<-as.matrix(new_ids_length)##
    short<-which(new_ids_length< pseudocell.size -1 )##
    new_good_ids<-as.matrix(new_ids_length[-short,])##
    result<-t(new.data)[,rownames(new_good_ids)]
    rownames(result)<-rownames(Inter)

    Tool(object) <- list(pseuoData = result,metaData = new_ids)
    #object <- LogSeuratCommand(object, return.command = TRUE)
    return(object)
}

CreatePseudoCellSeuratObject <- function(object, cluster_use, 
                                         pseudocell.size = 10, 
                                         assay_use="RNA", slot_use="counts"){
    seu.pc<-PseudoCell(object=object, 
                       assay_use = assay_use,slot_use = slot_use,
                       cluster_use =cluster_use, pseudocell.size=20 )
    
    avg.counts<-seu.pc@tools$PseudoCell$pseuoData
    
    meta    <-seu.pc@tools$PseudoCell$metaData
    meta$bc <-rownames(meta)
    meta$ann<-object@meta.data[meta$bc, cluster_use]
    meta<-meta[match( unique(meta$new_ids),   meta$new_ids  ), ]
    rownames(meta)<-meta$new_ids
    meta<-meta[colnames(avg.counts), ]  # make pseudocells' order consistent with the matrix    
    
    seu <- CreateSeuratObject(counts = avg.counts,
                    project = "PseudoCells",
                    assay = "RNA",
                    meta.data = meta,
                    )
    return (seu)
}
