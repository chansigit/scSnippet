Create10xObject <- function(path, proj.name, min.cells=3, min.features=200, mito.regex="^MT-" ){
	raw <- CreateSeuratObject(counts= Read10X(data.dir= path), 
		                          project = proj.name, 
		                          min.cells = min.cells, min.features = min.features) 
	raw <- PercentageFeatureSet(object = raw, pattern = mito.regex, col.name = "percent.mt") 
	raw <- RenameCells(object= raw, add.cell.id=proj.name)
	return (raw)
}


MergeMetaData <- function(... ){
	metadata <- data.frame()
	objects  <- list(...)
	for (obj in objects){
		metadata <-rbind(metadata, obj[[]])
	}
	metadata 
}

QCPlot <- function(..., feature){
	metadata <- MergeMetaData(...)
	plot <-ggviolin(metadata,
         x="orig.ident",  xlab="",
         y=feature, ylab=feature,
         color="orig.ident",palette="jco",  width=1.2,
         add=c("jitter"), add.params=list(size=0.001, jitter=0.25))

	if (feature!="percent.mt"){
		plot<-plot+ rremove("legend") + 
		      scale_y_continuous(trans="log1p", limits=c(200, 7000),
                                 breaks=c(200,300,400,500,750,1000,1500,2000,3000,4000,5000,6000,7000)) + 
              stat_summary(fun.y=median, geom="point", shape=23, size=2)
	}else{
		plot<-plot+ rremove("legend") + 
		      scale_y_continuous(breaks=c(5,10,15,20,25,30), limits=c(0.0, 40))+
              stat_summary(fun.y=median, geom="point", shape=23, size=2)
	}
	return(plot)
}
