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
