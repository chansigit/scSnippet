DefaultAssay(epi.combined)<-"RNA"
# feature by sample matrix, with row/column orders not set
# should be a matrix, with rownames and colnames
mat.orig<- as.matrix( GetAssayData(epi.combined, assay = "RNA", slot = "data") )

# retrieve the original metadata 
meta.orig <-epi.combined@meta.data[, c("ann210320","ann210320_treatment","treatment")]






hm.disp.order<-c(
'AT1_HDM',
'AT2_PBS',
'AT2_HDM',
'Ciliated_PBS',
'Ciliated_HDM',
'Club_PBS',
'Club_HDM',
'Goblet_HDM',
'Basal_HDM',
'Neuroendocrine_HDM'
)



# choose one metadata column as the sorting column, set the sorting criterion as the given order
meta <- (meta.orig%>% arrange(factor(ann210320_treatment, levels = hm.disp.order)) )

# generate cell orders
ordered.cells<-rownames(meta)

# a vector of features to show, whose order will be obeyed in the final figure
features.sel <-  epi.genes
# rearrange original matrix
mat <- mat.orig[features.sel, ordered.cells]





palette=c(
'AT2'="#ffd319",
'Ciliated'="#5AAE61",
'Club'="#95F5F5",
'AT1'="#A50026",
'Goblet'="#B2ABD2",
'Basal'="#05709D",
'Neuroendocrine'="#D6604D")



suppressPackageStartupMessages({
library(ComplexHeatmap)
library(circlize)
})

# annotations for samples (columns) in the rearrange matrix (mat)
# must be a character/numerical vector with the same length as ncol(mat)
colann <- HeatmapAnnotation(
    cluster   = meta$ann210320,
    treatment = meta$treatment,
    col = list(cluster   = palette, 
               treatment = c("PBS"="lightgray", "HDM" = "darkblue") ),
    annotation_legend_param=list(
        cluster = list(nrow=1),
        treatment = list(nrow=1)
    )
)



hm<-Heatmap(mat, name = "Normalized expression", 
        cluster_rows = F,
        cluster_columns = F, show_column_names=FALSE,
        cluster_row_slices=F,
        column_split=meta$ann210320, cluster_column_slices=F,
        col= colorRamp2(c(0, 1.5, 3), c("#486E9E", "white", "#D84B59")),
        column_title_rot=90, column_gap=unit(2, "mm"),
        top_annotation=colann, heatmap_legend_param = list(direction = "horizontal"),
        use_raster = FALSE)

options(repr.plot.width=15, repr.plot.height=14)
draw(hm,
     padding = unit(c(2, 2, 30, 10), "mm"),
     merge_legend = TRUE,
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
