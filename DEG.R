my_DE <-function(object, identity="RNA_snn_res.0.5", name="sample",
                 outpath.prefix, 
                 min.pct=0.25, 
                 logfc.threshold=0.25, 
                 pval.thres=0.01, 
                 adjpval.thres=0.005,...
                  ){
  # 0. setting clustering information and do DEG analysis
  Idents(object)<-identity
  tic("Perform differential expression analysis")
  markers <- FindAllMarkers(object, only.pos = FALSE, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, return.thresh = pval.thres)
  toc()

  # 1. filtering for significant marker genes (all, positive only, and negative only)
  markers.sig                  <- markers[      markers$p_val_adj<adjpval.thres , ]
  markers.sig["abs_avg_logFC"] <-     abs(      markers.sig["avg_logFC"]          )
  markers.sigpos               <- markers.sig[  markers.sig$"avg_logFC">0,        ]
  markers.signeg               <- markers.sig[  markers.sig$"avg_logFC"<0,        ]

  # 2. store all marker genes to file 
  fp_allmarkerlist <- paste0(outpath.prefix, 
                             "/MarkerList.", name, ".", identity, "_", 
                             "all.csv",      sep = "")
  write.csv(markers.sig,    file = fp_allmarkerlist)

  fp_posmarkerlist <- paste0(outpath.prefix, 
                             "/MarkerList.", name, ".", identity, "_", 
                             "pos.csv",      sep = "")
  write.csv(markers.sigpos, file = fp_posmarkerlist)

  fp_negmarkerlist <- paste0(outpath.prefix, 
                             "/MarkerList.", name, ".", identity, "_", 
                             "neg.csv",      sep = "")
  write.csv(markers.signeg, file = fp_negmarkerlist)

  # 3. prepare a list object to return
  DE.result                       <- list()  
  DE.result[["all_markers"]]      <- markers.sig
  DE.result[["positive_markers"]] <- markers.sigpos
  DE.result[["negative_markers"]] <- markers.signeg

  return(DE.result)
}

my_DEplot <- function(object, DE.result, w, h, outpath.prefix,
                      identity="RNA_snn_res.0.5", name="sample",
                      given.identity.order=NULL, topNGenes=40L, 
                      min.pct=0.25, logfc.threshold=0.25, pval.thres=0.01, adjpval.thres=0.005,...){
  DE.result[["all_markers"]]      -> markers.sig
  DE.result[["positive_markers"]] -> markers.sigpos
  DE.result[["negative_markers"]] -> markers.signeg

  markers.sig.top    <- markers.sig[ order( markers.sig$cluster       ,
                                           -markers.sig$abs_avg_logFC),   ] %>%
                         group_by(cluster) %>% top_n(n = topNGenes, wt = abs_avg_logFC)

  markers.sigpos.top <- markers.sigpos[ order( markers.sigpos$cluster   ,
                                              -markers.sigpos$avg_logFC), ] %>%
                         group_by(cluster) %>% top_n(n = topNGenes, wt = avg_logFC)

  markers.signeg.top <- markers.signeg[ order( markers.signeg$cluster   ,
                                              -markers.signeg$avg_logFC), ] %>%
                         group_by(cluster) %>% top_n(n = topNGenes, wt = avg_logFC)
  

  # ---------- Draw all markers on heatmap --------------
  fp <- paste0(outpath.prefix,
               "/Heatmap.", name, ".", identity, "_", 
               "allMarkers.pdf",  sep = "")
  pdf(fp, width=w, height=h)
  my_heatmap(object, genes=markers.sig.top, given.identity.order=given.identity.order, 
             annotation_height=2, group.by = identity)
  dev.off()

  # ---------- Draw positive markers on heatmap ----------
  fp <- paste0(outpath.prefix,
               "/Heatmap.", name, ".", identity, "_", 
               "posMarkers.pdf",  sep = "")
  pdf(fp, width=w, height=h)
  my_heatmap(object, genes=markers.sigpos.top, given.identity.order=given.identity.order, 
             annotation_height=2, group.by = identity)
  dev.off()

  # ---------- Draw negative markers on heatmap ----------
  fp <- paste0(outpath.prefix,
               "/Heatmap.", name, ".", identity, "_", 
               "negMarkers.pdf",  sep = "")
  pdf(fp, width=w, height=h)
  my_heatmap(object, genes=markers.signeg.top, given.identity.order=given.identity.order, 
             annotation_height=2, group.by = identity)
  dev.off()

}