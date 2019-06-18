# replace __object__ with other names

__object__.markers.sig   <- __object__.markers[__object__.markers$p_val_adj<0.005,]
__object__.markers.sig["abs_avg_logFC"] <- abs(__object__.markers.sig["avg_logFC"])
__object__.markers.sigpos <- __object__.markers.sig[ __object__.markers.sig$"avg_logFC">0, ]
__object__.markers.signeg <- __object__.markers.sig[ __object__.markers.sig$"avg_logFC"<0, ]

__object__.markers.sig.top   <- __object__.markers.sig[ order( __object__.markers.sig$cluster,
                                                -__object__.markers.sig$abs_avg_logFC), ] %>%
                         group_by(cluster) %>% top_n(n = 40, wt = abs_avg_logFC)

__object__.markers.sigpos.top<- __object__.markers.sigpos[ order( __object__.markers.sigpos$cluster,
                                                   -__object__.markers.sigpos$avg_logFC), ] %>%
                         group_by(cluster) %>% top_n(n = 40, wt = avg_logFC)

__object__.markers.signeg.top<- __object__.markers.signeg[ order( __object__.markers.signeg$cluster,
                                                   -__object__.markers.signeg$avg_logFC), ] %>%
                         group_by(cluster) %>% top_n(n = 40, wt = avg_logFC)