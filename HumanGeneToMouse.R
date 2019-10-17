# Basic function to convert human to mouse gene names
# https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
ConvertHumanGeneListToMM <- function(x){
  
  # Load human ensembl attributes
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # Load mouse ensembl attributes
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  # Link both datasets and retrieve mouse genes from the human genes
  genes.list = biomaRt::getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)
  
  # Get unique names of genes (in case gene names are duplicated)
  mouse.gene.list <- unique(genes.list[, 2])
  
  # # Print the first 6 genes found to the screen
  # print(head(mouse.gene.list))
  return(mouse.gene.list)
}
