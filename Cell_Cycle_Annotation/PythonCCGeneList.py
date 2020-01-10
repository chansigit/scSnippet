import urllib
response = urllib.request.urlopen("https://raw.githubusercontent.com/chansigit/scSnippet/master/Cell_Cycle_Annotation/mm.g2m.tsv")
mm_g2m_genes = [eval(str.strip(str(line, encoding = "utf-8"))) for line in response]

response = urllib.request.urlopen("https://raw.githubusercontent.com/chansigit/scSnippet/master/Cell_Cycle_Annotation/mm.s.tsv")
mm_s_genes = [eval(str.strip(str(line, encoding = "utf-8"))) for line in response]
