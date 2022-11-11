from collections import OrderedDict
cell_phase_genes = OrderedDict()

cell_phase_genes["G1-S"] = pd.Series(
    ["Arglu1","Brd7","Cdc6","Clspn","Esd","Gins2","Gmnn",
    "Luc7l3","Mcm5","Mcm6","Nasp","Pcna","Pnn",
    "Slbp","Srsf7","Ssr3","Zranb2"])

cell_phase_genes["S"] = pd.Series(
    [
"Asf1b","Calm2","Calm3","Calm1","Cdc45","Cdca5",
"Cenpm","Dhfr","Ezh2","Fen1","Pkmyt1","Prim1",
"Rfc2","Rpa2","Rrm2","Rsrc2","Srsf5","Svip","Top2a",
"Tyms","Ube2t","Zwint"])

cell_phase_genes["G2-M"] = pd.Series(
    ["Aurkb","Bub3","Ccna2","Ccnf","Cdca2","Cdca3",
    "Cdca8","Cdk1","Ckap2","Dcaf7","Hmgb2","Kif5b",
    "Kif20b","Kif22","Kif23","Kifc5b","Kifc1","Kpna2",
    "Gm10184","Lbr","Mad2l1","Mnd1","Ndc80","Nucks1",
    "Nusap1","Pif1","Psmd11","Psrc1","Smc4","Timp1",
    "Top2a","Tubb5","Tubb4b","Vps25"])

cell_phase_genes["M"] = pd.Series(
    ["Anp32b","Anp32e","Arl6ip1","Aurka","Birc5","Bub1",
    "Ccna2","Ccnb2","Cdc20","Cdc27","Cdc42ep1","Cdca3",
    "Cenpa","Cenpe","Cenpf","Ckap2","Ckap5","Cks1b",
    "Cks2","Depdc1a","Dlgap5","Dnaja1","Dnajb1","Grk6",
    "Gtse1","Hmg20b","Hmgb3","Hmmr","Hspa8","Kif2c",
    "Kif5b","Kif20b","Lbr","Mki67","Mzt1","Nuf2",
    "Nusap1","Pbk","Plk1","Prr11","Psmg3","Pwp1",
    "Rad51c","Rbm8a","Rbm8a2","Rnf126","Rnps1","Rrp1","Sfpq",
    "Smarcb1","Srsf3","Tacc3","Thrap3","Tpx2","Tubb4b","Usp16","Ywhah","Zfp207"])

cell_phase_genes["M-G1"] = pd.Series(
    ["Amd2","Amd1","Anp32e","Cbx3","Cdc42","Cnih4",
    "Cwc15","Dkc1","Dnajb6","Eif4e","Fxr1","Grpel1",
    "Gspt1","Hmg20b","Hspa8","Ilf2","Kif5b","Kpnb1",
    "Larp1","Lyar","Morf4l2","Mrpl19","Mrps2","Mrps18b",
    "Nucks1","Prc1","Ptms","Pttg1","Ran","Rheb",
    "Srsf3","Syncrip","Taf9","Ak6","Tmem138","Top1","Troap","Zfp593"])


dyn.pp.cell_cycle_scores(alp1,gene_list=cell_phase_genes)
