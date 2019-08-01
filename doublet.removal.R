DoubletRemovalParameters <- function(seu, PCs=1:30, use.SCT=FALSE, num.cores=1, quietly=TRUE){
    tic("sweep parameters")
    # calculate parameters
    if (quietly==TRUE){
        invisible(capture.output(sweep.res.list <- paramSweep_v3(seu, PCs = PCs, sct=use.SCT, num.cores=num.cores)))
        invisible(capture.output(sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)))
        ff <- tempfile()
        png(filename=ff)
        invisible(capture.output(bcmvn <- find.pK(sweep.stats)))
        dev.off()
        unlink(ff)
    }else{
        sweep.res.list <- paramSweep_v3(seu, PCs = PCs, sct=use.SCT, num.cores=num.cores)
        sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)
        ff <- tempfile()
        png(filename=ff)
        bcmvn <- find.pK(sweep.stats)
        dev.off()
        unlink(ff)
    }
    toc()
    # choose parameters
    maxBCmetric    <- max(bcmvn$BCmetric, na.rm = TRUE)
    pK.maxBCmetric <- as.numeric(as.character(bcmvn[bcmvn$BCmetric==maxBCmetric, ]$pK))
    return(list(pK=pK.maxBCmetric, bcmvn=bcmvn))
}

bcmvnPlot<-function(bcmvn, maxpos=0.08){
    data<-data.frame(pK=as.numeric(as.character(bcmvn$pK)),
              BCmvn=as.numeric(as.character(bcmvn$BCmetric)))
    ggplot(data=data, aes(x=pK, y=BCmvn)) +
    geom_point(color = "#20639B")+
    geom_line(color = "#3CAEA3")+
    theme(axis.text.x = element_text(angle = 45,size=6, hjust = 1))+ 
    scale_x_continuous(name="pK", limits = c(0.001,0.3),
        breaks=c(0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
                 0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,
                 0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3))+
    geom_vline(xintercept = maxpos, linetype="dashed", color = "red", size=0.5)->plot
    plot
}

DoubletRemoval <-function(seu, identity, pK, doublet.rate, pN=0.25, PCs=1:30, use.SCT=FALSE){
    # compute doublet scores
    tic("Removing doublets")
    annotations    <- seu@meta.data$identity  ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    homotypic.prop <- modelHomotypic(annotations) 
    nExp_poi       <- round(doublet.rate*length(colnames(x = seu)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj   <- round(nExp_poi*(1-homotypic.prop))
    seu.scored     <- doubletFinder_v3(seu, PCs =PCs, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = use.SCT)
    toc()
    
    # pick out doublets
    cname <-colnames(seu.scored[[]])
    DF<-cname[grep('^DF',cname)] 
    seu.scored[["doublet"]] <- as.numeric(seu.scored[[DF]]=="Doublet")

    # remove doublets
    seu.removed <- subset(seu.scored, subset = doublet != 1)
    return(list(removed=seu.removed, original=seu.scored))
}
