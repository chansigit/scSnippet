dissociation.genes.hs<-c("ACTG1","ANKRD1","ARID5A","ATF3","ATF4","BAG3","BHLHE40",
"CCNL1","CCRN4L","CEBPB","CEBPD","CEBPG","CSRNP1","CXCL1","CYR61","DCN","DDX3XX","DDX5",
"DES","DNAJA1","DNAJB1","DNAJB4","DUSP1","DUSP8","EGR1","EGR2","EIF1","EIF5",
"ERF","ERRFI1","FAM132B","FOS","FOSB","FOSL2","GADD45A","GADD45G","BRD2","BTG1","BTG2",
"GCC1","GEM","H3F3B","HIPK3","HSP90AA1","HSP90AB1","HSPA1A","HSPA1B","HSPA5",
"HSPA8","HSPB1","HSPE1","HSPH1","ID3","IDI1","IER2","IER3","IER5","IFRD1","IL6",
"IRF1","IRF8","ITPKC","JUN","JUNB","JUND","KCNE4","KLF2","KLF4","KLF6","KLF9",
"LITAF","LMNA","MAFF","MAFK","MCL1","MIDN","MIR22HG","MT1",
"MT2","MYADM","MYC","MYD88","NCKAP5L","NCOA7","NFKBIA","NFKBIZ","NOP58","NPPC",
"NR4A1","ODC1","OSGIN1","OXNAD1","PCF11","PDE4B","PER1","PHLDA1","PNP","PNRC1",
"PPP1CC","PPP1R15A","PXDC1","RAP1B","RASSF1","RHOB","RHOH","RIPK1","SAT1X","SBNO2",
"SDC4","SERPINE1","SKIL","SLC10A6","SLC38A2","SLC41A1","SOCS3","SQSTM1","SRF",
"SRSF5","SRSF7","STAT3","TAGLN2","TIPARP","TNFAIP3","TNFAIP6","TPM3","TPPP3",
"TRA2A","TRA2B","TRIB1","TUBB4B","TUBB6","UBC","USP2","WAC","ZC3H12A",
"ZFAND5","ZFP36","ZFP36L1","ZFP36L2","ZYX")

dissociation.genes.mm<-c("Actg1", "Ankrd1", "Arid5a", "Atf3", "Atf4", "Bag3", "Bhlhe40", 
	"Brd2", "Btg1", "Btg2", "Ccnl1", "Ccrn4l", "Cebpb", "Cebpd", "Cebpg", 
	"Csrnp1", "Cxcl1", "Cyr61", "Dcn", "Ddx3xX", "Ddx5", "Des", "Dnaja1", 
	"Dnajb1", "Dnajb4", "Dusp1", "Dusp8", "Egr1", "Egr2", "Eif1", "Eif5", 
	"Erf", "Errfi1", "Fam132b", "Fos", "Fosb", "Fosl2", "Gadd45a", "Gcc1", 
	"Gem", "H3f3b", "Hipk3", "Hsp90aa1", "Hsp90ab1", "Hspa1a", "Hspa1b", 
	"Hspa5", "Hspa8", "Hspb1", "Hsph1", "Id3", "Idi1", "Ier2", "Ier3", 
	"Ifrd1", "Il6", "Irf1", "Irf8", "Itpkc", "Jun", "Junb", "Jund", 
	"Klf2", "Klf4", "Klf6", "Klf9", "Litaf", "Lmna", "Maff", "Mafk", 
	"Mcl1", "Midn", "Mir22hg", "Mt1", "Mt2", "Myadm", "Myc", "Myd88", 
	"Nckap5l", "Ncoa7", "Nfkbia", "Nfkbiz", "Nop58", "Nppc", "Nr4a1", 
	"Odc1", "Osgin1", "Oxnad1", "Pcf11", "Pde4b", "Per1", "Phlda1", 
	"Pnp", "Pnrc1", "Ppp1cc", "Ppp1r15a", "Pxdc1", "Rap1b", "Rassf1", 
	"Rhob", "Rhoh", "Ripk1", "Sat1X", "Sbno2", "Sdc4", "Serpine1", 
	"Skil", "Slc10a6", "Slc38a2", "Slc41a1", "Socs3", "Sqstm1", "Srf", 
	"Srsf5", "Srsf7", "Stat3", "Tagln2", "Tiparp", "Tnfaip3", "Tnfaip6", 
	"Tpm3", "Tppp3", "Tra2a", "Tra2b", "Trib1", "Tubb4b", "Tubb6", 
	"Ubc", "Usp2", "Wac", "Zc3h12a", "Zfand5", "Zfp36", "Zfp36l1", 
	"Zfp36l2", "Zyx", "Gadd45g", "Hspe1", "Ier5", "Kcne4")

DefaultAssay(meni.integrated) <- "RNA"
gset <- dissociation.genes.hs
gset <- gset[gset %in% rownames(meni.integrated)]
meni.integrated[["percent.disso"]]<-PercentageFeatureSet(meni.integrated,features = gset)
