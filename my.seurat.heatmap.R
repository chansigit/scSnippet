###############################################################################
#
# Author:  Sijie Chen    
# Version: 1.0.1
# Time:    2019/07/26
#
###############################################################################
my.seurat.heatmap<-function(object, genes, slot="scale.data", 
                     group.by="RNA_snn_res.0.5", 
                     given.identity.order=NULL,
                     annotation_height=1, 
                     identity.palette=c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                                        '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                                        '#cab2d6','#6a3d9a','#ffff99','#b15928',
                                        '#980043','#02818a','#00441b','#984ea3'),
                     hm.color.lower="#FF08BD", hm.color.mid="#000000", hm.color.upper="#FFF42F",
                     hm.color.lb=-1.5, hm.color.ub=1.5, ...){
    library(ComplexHeatmap)
    # ==========================================================================
    # Filter undetected genes out
    ## find existed given genes in the supplied slot
    gene.all  <- rownames(GetAssayData(object = object, slot = slot))
    gene.use  <- intersect(genes, gene.all)
    gene.miss <- setdiff  (genes, gene.all)
    if (length(gene.miss)>0){ # if any missing genes
        warning(paste("Genes not found: ", paste(gene.miss,collapse="  "),"\n"))	
    }


    # ==========================================================================
    # Transform sparse matrix into matrix
    dgeSparse <- GetAssayData(object = object, slot = slot)
    dge       <- as.matrix(dgeSparse[gene.use,])


    # ==========================================================================
    # Reorder cells by (given) group
    meta <- object[[]] 
    if (is.null(given.identity.order)){
        given.identity.order<-levels(as.factor(meta[,group.by]))
    }
    cell.barcode.order    <- c()
    cell.identity.ordered <- c()
    for (identity in given.identity.order){ 
        subgroup <- meta[meta[,group.by]==identity,  ] 
        bc       <- rownames(subgroup)
        cell.identity.ordered <- c(cell.identity.ordered, 
                                   as.character(subgroup[,group.by])) # convert to character, or it wrongly convert 0 to 1 due to thr factor datatype
        cell.barcode.order    <- c(cell.barcode.order, bc)
    }
    dge<-dge[, cell.barcode.order]


    # ==========================================================================
    # build heatmap annotations
    identities<-unique(given.identity.order)
    colors.used<-palette[1:length(identities)]
    names(colors.used)<-identities
    ha <- HeatmapAnnotation(identity= cell.identity.ordered, 
                            col     =list(identity=colors.used),
                            annotation_height  = unit(annotation_height, "cm"))


    # ========================================================================== 
    # draw heatmap
    require(circlize)
    hm <- Heatmap(dge, name=slot, 
    	        col=colorRamp2(c(hm.color.lb, 0, hm.color.ub),
    	                       c(hm.color.lower, hm.color.mid, hm.color.upper)),
                cluster_rows=FALSE, cluster_columns = FALSE, 
                show_column_names=FALSE,
                top_annotation = ha, ...)
    draw(hm, heatmap_legend_side = "right")
}
