Idents(seu)<-"RNA_snn_res.0.5"
tic("DE")
markers <- FindAllMarkers(seu, only.pos = FALSE, min.pct = 0.25, 
                          logfc.threshold = 0.25, return.thresh = 0.01)
toc()



# Filter
markers <- markers[ markers$p_val_adj<0.01, ]
# Sort
markers <- markers[ order(markers$cluster, -markers$avg_logFC), ]

# Annotate
## load gene list from github (THIS IS FOR MOUSE)
load(url("https://github.com/chansigit/SSAT/raw/master/mm.cellsurfacemarker.rda"))
load(url("https://github.com/chansigit/SSAT/raw/master/mm.secretory.rda"))
load(url("https://github.com/chansigit/SSAT/raw/master/mm.tf.rda"))
## load annotation functions
source("https://raw.github.com/chansigit/SSAT/master/annotate.genelist.R")
markers<-annotate.genelist(markers, tf=mm.tf, surface=mm.cellsurfacemarker, secretory=mm.secretory)

# Visualize
library(devtools)
markers.viz  <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
markers.viz2 <- markers[markers$is.tf=="TF"| markers$is.secretory=="Secretory"|markers$is.surface=="Surface", ]%>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)  
options(repr.plot.width=15, repr.plot.height=16,repr.plot.resolution=300)
DoHeatmap(seu, features = markers.viz2$gene) + NoLegend()
