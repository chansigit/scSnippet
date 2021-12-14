library(biomaRt)
x = as.character(ref.ImmGen@NAMES)

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
rat   = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")

genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                 values = x , mart = mouse, 
                 attributesL = c("rgd_symbol"), martL = rat, uniqueRows=T)


# query what databases are available 
mart = useMart('ensembl')
listDatasets(mart) 

# query what columns are available
listAttributes(mart = rat)
