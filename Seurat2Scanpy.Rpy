#----------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------          simplified          ------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------

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

adata.obsm["X_umap"]=np.array(adata.obs[["UMAP_1","UMAP_2"]])
adata.obsm["X_pca"]=np.array(adata.obs[["PC_1","PC_2"]])
adata

#----------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------            updated           ------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------
sample<-"A2inner"
seu     <- meni.list[[sample]]

tic("embeddings")
  dr.umap <- Embeddings(seu,reduction ="umap")[rownames(seu@meta.data), 1:2]
  dr.pca  <- Embeddings(seu,reduction ="pca") [rownames(seu@meta.data) ,1:2]
  seu@meta.data <- cbind(seu@meta.data,dr.umap, dr.pca)
toc()

tic("extracting dge,var,obs")
  dge <- t(as.matrix(GetAssayData(seu, assay = "RNA", slot = "counts")))
  var <- colnames(dge)
  obs <- seu@meta.data
toc()


tic("writing dge,var,obs")
dge <- t(as.matrix(GetAssayData(seu, assay = "RNA", slot = "counts")))
var <- colnames(dge)
obs <- seu@meta.data
write.table(dge, sep = "\t", row.names =T, col.names = T, file=paste(sample, "_dge.tsv" ,sep=""))
write.table(var, sep = "\t" ,row.names =F, col.names = F, file=paste(sample, "_var.tsv" ,sep=""))
write.table(obs, sep = "\t", row.names =T, col.names = T, file=paste(sample, "_obs.tsv" ,sep=""))
toc()



X   = pd.read_csv("./*_dge.tsv",delimiter='\t',encoding='utf-8').to_numpy()
var = pd.DataFrame(index=pd.read_csv("./*_var.tsv",delimiter='\t',encoding='utf-8',header=None, names=["gene"]).to_numpy().flatten() )
obs = pd.read_csv("./*_obs.tsv",delimiter='\t',encoding='utf-8')
adata = ad.AnnData(X, obs=obs, var=var, dtype='int32')
adata.var_names_make_unique()
adata.obsm["X_umap"]=np.array(adata.obs[["UMAP_1","UMAP_2"]])
adata.obsm["X_pca"]=np.array(adata.obs[["PC_1","PC_2"]])
adata
