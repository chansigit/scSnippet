# obtain DEGs
plan("sequential")
plan("multiprocess", workers = 12)
options(future.globals.maxSize = 3*1000 * 1024^2)

Idents(LTi)<-"LTi_identity"
tic("Two-group DE analysis")
deg <- FindMarkers(LTi, ident.1 = "FLv1_5", ident.2 = "FSI1_9", min.pct = 0.25)
toc()

plan("sequential")



# filtering DEGs
adjpval.thres<-0.005
topNGenes<-100
markers <- deg
markers["gene"]<-rownames(deg)
markers.sig                  <- markers[      markers$p_val_adj<adjpval.thres , ]
markers.sig["abs_avg_logFC"] <-     abs(      markers.sig["avg_logFC"]          )
markers.sigpos               <- markers.sig[  markers.sig$"avg_logFC">0,        ]
markers.signeg <- markers.sig[ markers.sig$"avg_logFC"<0, ]


markers.sig.top    <- markers.sig[ order(-markers.sig$abs_avg_logFC),   ] %>%
                         top_n(n = topNGenes, wt = abs_avg_logFC)
markers.sig.top    <- markers.sig.top[order(-markers.sig.top$"avg_logFC"),]

markers.sigpos.top <- markers.sigpos[ order(-markers.sigpos$avg_logFC), ] %>%
                         top_n(n = topNGenes, wt = avg_logFC)

markers.signeg.top <- markers.signeg[ order(-markers.signeg$avg_logFC), ] %>%
                           top_n(n = topNGenes, wt = avg_logFC)
