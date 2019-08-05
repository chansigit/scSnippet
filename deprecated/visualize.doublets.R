DF.param     <-
seu.original <-
seu          <-
Idents(seu)<-"RNA_snn_res.0.5"


p1<-FeaturePlot(seu,reduction = "umap", features = "nCount_RNA")+
    theme(legend.position = c(0.1,0.75))+ggtitle("nCount_RNA")
p2<-DimPlot(seu.original, reduction ="umap")+
    theme(legend.position = c(0.1,0.55))+ggtitle("clusters")
p3<-FeaturePlot(seu.original,reduction ="umap", features = "doublet")+
    theme(legend.position = c(0.1,0.75))+ggtitle("doublet")
p4<-bcmvnPlot(bcmvn=DF.param$bcmvn, maxpos =DF.param$pK )+ggtitle("bcmvn")

w<-12; h<-4;
options(repr.plot.width = w, repr.plot.height = h, repr.plot.res=300)
ggarrange(p1,p2,p3,p4,nrow=1,ncol=4)
