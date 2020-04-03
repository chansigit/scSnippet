tic()
dge <- t(as.matrix(GetAssayData(seu, assay = "RNA", slot = "counts")))
var <- colnames(dge)
obs <- seu@meta.data
write.table(dge, sep = "\t", row.names =T, col.names = T, file="./*_dge.tsv")
write.table(var, sep = "\t" ,row.names =F, col.names = F, file="./*_var.tsv")
write.table(obs, sep = "\t", row.names =T, col.names = T, file="./*_obs.tsv")
toc()


X   = pd.read_csv("./*_dge.tsv",delimiter='\t',encoding='utf-8').to_numpy()
var = pd.DataFrame(index=pd.read_csv("./*_var.tsv",delimiter='\t',encoding='utf-8',header=None, names=["gene"]).to_numpy().flatten() )
obs = pd.read_csv("./*_obs.tsv",delimiter='\t',encoding='utf-8')
adata = ad.AnnData(X, obs=obs, var=var, dtype='int32')
adata.var_names_make_unique() 
adata
