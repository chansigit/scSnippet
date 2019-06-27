w<-6
h<-5
options(repr.plot.width=w,repr.plot.height=h,repr.plot.resolution=400)
ggbarplot(
    as.data.frame(table(LTi[["LTi_identity"]])),
    x = "Var1", y = "Freq", label=TRUE,
    xlab="identity", ylab="#Cell",
    fill="Var1",
    palette = c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                '#cab2d6','#6a3d9a','#ffff99','#b15928',
                '#980043','#02818a','#00441b','#984ea3'))+
guides(fill=guide_legend(title="Identities:"))