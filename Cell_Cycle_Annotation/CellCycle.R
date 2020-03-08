# mouse  mm.s.genes  mm.g2m.genes
load(url("https://raw.githubusercontent.com/chansigit/scSnippet/master/Cell_Cycle_Annotation/mm.g2m.rda"))
load(url("https://raw.githubusercontent.com/chansigit/scSnippet/master/Cell_Cycle_Annotation/mm.s.rda"))
CellCycleScoring(seu, s.features=mm.s, g2m.features=mm.g2m)

# human
CellCycleScoring(seu, s.features=Seurat::cc.genes$s.genes, g2m.features=Seurat::cc.genes$g2m.genes)
