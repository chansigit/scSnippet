my_heatmap<-function(object, genes, slot="scale.data", 
                     group.by="RNA_snn_res.0.5", given.identity.order=NULL){
  # filter undetected genes out
  gene.all  <- rownames(GetAssayData(object = object, slot = slot))
  gene.use  <- intersect(markers, gene.all)
  gene.miss <- setdiff  (markers, gene.all)
  warning(paste("Genes not found: ", paste(gene.miss,collapse="  ")))
    
  # transform sparse matrix into matrix
  dgeSparse <- GetAssayData(object = object, slot = slot)
  dge       <- as.matrix(dgeSparse[gene.use,])

  # reorder cells by (given) group
  meta <- object[[]] 
  if (is.null(given.identity.order)){
    given.identity.order<-levels(meta[,group.by])
  }
  cell.barcode.order <- c()

  for (identity in given.identity.order){ 
    bc   <- rownames( meta[meta[,group.by]==identity,  ] )
    cell.barcode.order<-c(cell.barcode.order,bc)
  }

  dge<-dge[, cell.barcode.order]

  # draw heatmap
  hm<-Heatmap(dge,
              cluster_rows=FALSE, cluster_columns = FALSE, 
              show_column_names=FALSE)
  draw(hm, heatmap_legend_side = "right")
}

options(repr.plot.width=21,repr.plot.height=17,repr.plot.resolution=500)
my_heatmap(ASI, genes=markers, given.identity.order=0:11)

# see here for more info
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#clustering