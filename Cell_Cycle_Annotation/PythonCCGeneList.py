import urllib
response = urllib.request.urlopen("https://raw.githubusercontent.com/chansigit/scSnippet/master/Cell_Cycle_Annotation/mm.g2m.tsv")
mm_g2m_genes = [eval(str.strip(str(line, encoding = "utf-8"))) for line in response]

response = urllib.request.urlopen("https://raw.githubusercontent.com/chansigit/scSnippet/master/Cell_Cycle_Annotation/mm.s.tsv")
mm_s_genes = [eval(str.strip(str(line, encoding = "utf-8"))) for line in response]



cell_cycle_genes = [x.strip() for x in open('./regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
#cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]