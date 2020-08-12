LRDB<-read.csv("LR_DASkelly.csv")


seu<-wti

#仅仅保留数据中有的LR pair
ligand.indata   <- unique(LRDB$ligand_symbol[LRDB$ligand_symbol %in%rownames(seu)])
receptor.indata <- unique(LRDB$receptor_symbol[LRDB$receptor_symbol %in%rownames(seu)])
lrpair.indata   <-(LRDB$ligand_symbol %in% ligand.indata)&(LRDB$receptor_symbol %in% receptor.indata)
LRDB            <- LRDB[lrpair.indata, ]



# 随机Permute细胞类型变量，生成cell type的空分布，存在Backgrounds变量里，供后面计算pvalue用
tic()
NPerm<-10000
Backgrounds<-list()
library(progress)
pb <- progress_bar$new(total = NPerm,format = "  generating [:bar] :percent in :elapsed eta: :eta",clear = T, width= 60,show_after=50)

for (i in 1:NPerm){
    pb$tick()
    #progress(i)
    # permute the cell identities
    identity<-seu$celltype0627
    names(identity)<-NULL
    set.seed(i)
    seu$random<-sample(identity)

    # compute the average expression after permutation
    Idents(seu)<-"random"
    expr <- AverageExpression(object = seu, assays = "RNA" ,slot = "data",verbose = F)$RNA
    expr[expr<expr.cutoff] <- 0
    expr <- as(as.matrix(expr), "dgCMatrix") 
    Backgrounds[[i]]<-list("expr.ligand"  =expr[LRDB$ligand_symbol,  ], 
                           "expr.receptor"=expr[LRDB$receptor_symbol,]) 
    if (i==NPerm) cat("Done!\n")
}
toc()


# 枚举细胞类型,计算互作分数，并使用背景分布计算P-value
expr.cutoff  <- 1.0
Idents(seu)<-"celltype0627"
expr <- AverageExpression(object = seu, assays = "RNA" ,slot = "data",verbose = F)$RNA
expr[expr<expr.cutoff] <- 0
expr <- as(as.matrix(expr), "dgCMatrix") 
expr.ligand  <-expr[LRDB$ligand_symbol,  ]
expr.receptor<-expr[LRDB$receptor_symbol,]


crosstalk<-list()
for (ct1 in levels(seu)){
    for (ct2 in levels(seu)){
        cat(ct1,"-> ⊃-",ct2,"\n")
        flush.console()
        # crosstalk between two cell types
        inner.prod <- expr.ligand[,ct1,drop=F] * expr.receptor[,ct2,drop=F]
        row.names(inner.prod) <- LRDB$pair
        LR.scores  <- sort(inner.prod[inner.prod[,1]>0,],decreasing = T)
        
        # compute p values
        LR.scores.pval<-rep(0,length(LR.scores))
        for(i in 1:NPerm){
            expr.ligand  <-Backgrounds[[i]]$expr.ligand
            expr.receptor<-Backgrounds[[i]]$expr.receptor
            # compute the background ligand-receptor interaction scores
            inner.prod <- expr.ligand[,ct1,drop=F] * expr.receptor[,ct2,drop=F]
            row.names(inner.prod) <- LRDB$pair
            LR.scores.null<-inner.prod[names(LR.scores),]
            LR.scores.pval<-LR.scores.pval+(as.numeric(LR.scores.null>LR.scores))
        }
        LR.scores.pval<-LR.scores.pval/NPerm
        
        
        crosstalk[[ct1]][[ct2]] <- data.frame(score=LR.scores, pval=LR.scores.pval)%>%arrange(pval)
    }
    
}
