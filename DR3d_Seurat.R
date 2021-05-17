# Old version
selected.genes<-c("Flt3","Id2","Zbtb16","Irf8")
expression<-as.data.frame(t(seu.imputed@assays$alra@data[selected.genes,]))

options(repr.plot.width=8, repr.plot.height=7)
# for vigenettes of gg3D, see http://htmlpreview.github.io/?https://github.com/AckerDWM/gg3D/blob/master/gg3D-vignette.html
# To install the R package gg3D, run devtools::install_github("AckerDWM/gg3D")
library("gg3D")
theta=-70
phi=30


ggplot(cbind(seu[[]],
             seu@reductions$pca@cell.embeddings[,1:3],
             expression
            ),
       aes(x=PC_1, y=PC_2, z=PC_3, color=Zbtb16)
      ) + 
axes_3D(theta=theta, phi=phi) +
stat_3D(theta=theta, phi=phi) +
axis_labs_3D(theta=theta, phi=phi, size=3, 
             hjust=c(1,1,1.2,1.2,1.2,1.2), 
             vjust=c(-.5,-.5,-.2,-.2,1.2,1.2)) +
labs_3D(theta=theta, phi=phi, 
        hjust=c(1,0,0), vjust=c(1.5,1,-.2),
        labs=c("PC1", "PC2", "PC3")) +
theme_void()

# New version
seu = RunUMAP(seu, dims = 1:20, n.components = 3L, reduction.name = "umap3d", verbose=FALSE)

# feature plot
o(20,20)
library(plotly)
dim1 <-(seu@reductions$umap3d@cell.embeddings)[,"umap3d_1"] #c(0.04499301	,2,3)
dim2 <-(seu@reductions$umap3d@cell.embeddings)[,"umap3d_2"]#c(1.93572454,  1.3,4)
dim3 <-(seu@reductions$umap3d@cell.embeddings)[,"umap3d_3"] #c(-2.129847,   4,5.5)
plot_ly(x=as.vector(dim1), y=as.vector(dim2), z=as.vector(dim3), 
        color=seu@assays$RNA@data["CD8A",],
        type="scatter3d",
        mode="markers",size=0.1)

# dimplot
library(plotly)
dim1 <-(seu@reductions$umap3d@cell.embeddings)[,"umap3d_1"] #c(0.04499301	,2,3)
dim2 <-(seu@reductions$umap3d@cell.embeddings)[,"umap3d_2"]#c(1.93572454,  1.3,4)
dim3 <-(seu@reductions$umap3d@cell.embeddings)[,"umap3d_3"] #c(-2.129847,   4,5.5)
plot_ly(x=as.vector(dim1), y=as.vector(dim2), z=as.vector(dim3), 
        color=seu@meta.data$RNA_snn_res.0.7,
        type="scatter3d",
        mode="markers",size=0.1)
