library(RColorBrewer)
n <- 35
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',] # one of 'div', 'qual', 'seq'
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector=sample(col_vector)
pie(rep(1,n), col=sample(col_vector, n))
