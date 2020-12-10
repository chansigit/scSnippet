TryUMAPs<- function(seu, group.by, n.neighbors=c(5,15,30,50), min.dist=c(0.001,0.01,0.05,0.1,0.5) ){
    palette.ann0507=c(
    "Chondrocytes.1"="#cc3300",#0 red
    "Myeloid cells"="#ff9900",#1 blue
    "Chondrocytes.2"="#ffcc00",#2 yellow
    "Endothelial cells"="#0066ff",#3 blue    
    "Lymphoid cells"="#494ca2",#5 purple blue
    "ACTA2+ cells"="#009933",#7 green
    "CDK1+ cells"="#162447"#8 pink   
    )
    
    pics <- list()
    
    for (nn in n.neighbors){
        for (md in min.dist){
            param<-paste("nn=",nn,",mindist=",md,sep="")
            print(param)
            flush.console()
            tic(param)
            seu <- RunUMAP(seu, n.neighbors=nn, uwot.sgd = TRUE,
                         min.dist=md,dims = 1:20,
                         verbose=FALSE)
            p <- DimPlot(seu, group.by=group.by, reduction="umap",label=T,
                         cols = palette.ann0507)+
                       NoLegend()+xlab(param)+ylab(seu[[]]$orig.ident[1])
            pics[[param]]<-p
            toc()
            flush.console()
        }
    }
    return (pics)
}
