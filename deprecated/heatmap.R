my_heatmap<-function(object, genes, slot="scale.data", annotation_height=1,
                     group.by="RNA_snn_res.0.5", given.identity.order=NULL){
  require(ComplexHeatmap)  
  # filter undetected genes out
  gene.all  <- rownames(GetAssayData(object = object, slot = slot))
  gene.use  <- intersect(genes, gene.all)
  gene.miss <- setdiff  (genes, gene.all)
  warning(paste("Genes not found: ", paste(gene.miss,collapse="  "),"\n"))
    
  # transform sparse matrix into matrix
  dgeSparse <- GetAssayData(object = object, slot = slot)
  dge       <- as.matrix(dgeSparse[gene.use,])

  # reorder cells by (given) group
  meta <- object[[]] 
  if (is.null(given.identity.order)){
    given.identity.order<-levels(meta[,group.by])
  }
  cell.barcode.order    <- c()
  cell.identity.ordered <- c()
  for (identity in given.identity.order){ 
    
    subgroup <- meta[meta[,group.by]==identity,  ] 
    bc       <- rownames(subgroup)
    cell.identity.ordered <- c(cell.identity.ordered, as.character(subgroup[,group.by])) # convert to character, or it wrongly convert 0 to 1 due to thr factor datatype
    cell.barcode.order    <- c(cell.barcode.order, bc)
    
  }
  dge<-dge[, cell.barcode.order]
  
  
  # build heatmap annotations
  identities<-unique(given.identity.order)
  colors.used<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                 '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                 '#cab2d6','#6a3d9a','#ffff99','#b15928',
                 '#980043','#02818a','#00441b','#984ea3')[1:length(identities)]
  names(colors.used)<-identities
  ha <- HeatmapAnnotation(identity= cell.identity.ordered, 
                          col     =list(identity=colors.used),
                          annotation_height  = unit(annotation_height, "cm"))
    
  # draw heatmap
  require(circlize)
  hm <- Heatmap(dge, name=slot, col=colorRamp2(c(-3, 0, 3), c("#FF08BD", "#000000", "#FFF42F")),
              cluster_rows=FALSE, cluster_columns = FALSE, 
              show_column_names=FALSE,
              top_annotation = ha)
  draw(hm, heatmap_legend_side = "right")
}
