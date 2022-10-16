DEG.analysis <- function(mat, ident.1=NULL, ident.2=NULL){
    col_names = colnames(mat)
    gene_names = rownames(mat)
    sel = startsWith(col_names, ident.1)
    s1 = col_names[sel]
    
    if (is.null(ident.2)){
        sel = !sel
    }else{
        sel = startsWith(col_names, ident.2)
    }
    s2 = col_names[sel]
    
    
    pvals = c()
    diffs = c()
    for(i in 1:nrow(mat)){
        p.value = wilcox.test(mat[i, s1], mat[i, s2], exact=F)$p.value
        pvals= c(pvals, p.value)
        diff = mean(mat[i, s1]) - mean(mat[i, s2])
        diffs=c(diffs, diff)
    }
    
    padj = p.adjust(pvals, method ='BH')
    
    results=data.frame(diffs,pvals, padj)
    rownames(results)=gene_names
    
    
    return (results)
    
}
