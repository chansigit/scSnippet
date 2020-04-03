tic()
dge <- t(as.matrix(GetAssayData(seu, assay = "RNA", slot = "counts")))
var <- colnames(dge)
obs <- seu@meta.data
write.table(dge, sep = "\t", row.names =T, col.names = T, file="./*_dge.tsv")
write.table(var, sep = "\t" ,row.names =F, col.names = F, file="./*_var.tsv")
write.table(obs, sep = "\t", row.names =T, col.names = T, file="./*_obs.tsv")
toc()
