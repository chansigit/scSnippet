my.seurat.scatterplot<- function(object, group.by, reduction="umap",contour=FALSE, pt.size=0.1, label.pt.size=5,
                           label.size=10, label.face="plain",label.rectangle=FALSE, palette=c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                '#cab2d6','#6a3d9a','#00441b','#b15928',
                '#980043','#02818a','#984ea3','#ffff99'),...){   
  # 1. bind dimensionality reduction info with metainfo
  meta_dr <- cbind( object[[]], Embeddings(object = object, reduction = reduction))
  
  # 2. choose dimensionality reduction name
  xcoord.name<-""
  ycoord.name<-""
  if (reduction=="umap"){
    xcoord.name<-"UMAP_1";ycoord.name<-"UMAP_2"
  }
  if (reduction=="tsne"){
    xcoord.name<-"tSNE_1";ycoord.name<-"tSNE_2"
  }
  if (reduction=="pca"){
    xcoord.name<-"PC_1";  ycoord.name<-"PC_2"
  }
    
  # 3. bind barcode info stored in rownames in case of getting lost when rownames were discarded
  meta_dr <- cbind("barcode"=rownames(meta_dr), meta_dr) 
  meta_dr <- meta_dr[order(meta_dr[,group.by], meta_dr[,xcoord.name], meta_dr[,ycoord.name]) , ]
  
  # 4. add text labels for each group
  for (identity in levels(as.factor(meta_dr[,group.by]))){
    # select ordered dataframe by groups
    subgroup<-meta_dr[ meta_dr[,group.by]==identity, ]
    
    # obtain the middle data item's cell barcode
    as.character(subgroup[nrow(subgroup)/2,]$barcode) -> cell_barcode
  
    # access the middle point with obtained cell barcode
    meta_dr[cell_barcode,"text_annotation"]<-identity
  }
    
  # 5. draw
  plot <- ggscatter(meta_dr, 
    x = xcoord.name, y = ycoord.name,    # x and y coordinates of points
    size = pt.size,                      # scatter plot dot size
    ellipse = (contour==TRUE) , 
    label = "text_annotation",repel=TRUE,font.label = c(label.size, label.face, "black"),label.rectangle=label.rectangle,
    color = group.by, 
    palette = palette )  + 
  guides(colour = guide_legend(override.aes = list(size=label.pt.size)))+ # legend dot size
  theme(legend.position = "bottom", # legend themes
        legend.title    = element_text(color = "blue", size = 10),
        legend.text     = element_text(color = "blue", size = 6))
  
  plot
}
