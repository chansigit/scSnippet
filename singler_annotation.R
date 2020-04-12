tic()
load("/data/hca/SingleRReference/HumanPrimaryCellAtlasData.rda")
load("/data/hca/SingleRReference/BlueprintEncodeData.rda")
load("/data/hca/SingleRReference/DatabaseImmuneCellExpressionData.rda")
load("/data/hca/SingleRReference/NovershternHematopoieticData.rda")
load("/data/hca/SingleRReference/MonacoImmuneData.rda")
toc()
SingleR_Annotation <- function(seu, reference="HPCA", use_local=T){
    if (use_local==F){
        if       (reference=="HPCA"){
            ref <- HumanPrimaryCellAtlasData()   
        }else if (reference=="BED"){
            ref <- BlueprintEncodeData()
        }else if (reference=="DbImmExp"){
            ref <- DatabaseImmuneCellExpressionData()
        }else if (reference=="Hemato"){
            ref <- NovershternHematopoieticData()
        }else if (reference=="Monaco"){
            ref <- MonacoImmuneData()
        }else if (reference=="ImmGen"){
            ref <- ImmGenData()
        }else if (reference=="MouseRNA"){
            ref <- MouseRNAseqData()
        }else{
            ref <- HumanPrimaryCellAtlasData()   
        }
    }else{
        if       (reference=="HPCA"){
            ref <- hpca
        }else if (reference=="BED"){
            ref <- blueprint
        }else if (reference=="DbImmExp"){
            ref <- dbimmexp
        }else if (reference=="Hemato"){
            ref <- hemato
        }else if (reference=="Monaco"){
            ref <- monaco
        }else if (reference=="ImmGen"){
            ref <- immgen
        }else if (reference=="MouseRNA"){
            ref <- mmrna
        }else{
            ref <- hpca
        }
    }
    
        
    mat       <- GetAssayData(seu, slot="data")
    pred.fine <- SingleR(test = mat, ref = ref, labels = ref$label.fine)
    pred.main <- SingleR(test = mat, ref = ref, labels = ref$label.main)
       
    seu[[paste0(reference,"_main")]]<-pred.main$pruned.labels
    seu[[paste0(reference,"_fine")]]<-pred.fine$pruned.labels
    return(seu) 
}
