# phylogeny - genus level

library(V.PhyloMaker)
sp.list <- read.csv("gbif_png_all_spp_names_corrected.csv", header = T, fileEncoding="UTF-8-BOM")
sp.list <- read.csv("PNG_GBIF_6924_phy.csv", header = T)
sptree <- phylo.maker(sp.list, tree = GBOTB.extended, nodes = nodes.info.1, output.sp.list = TRUE, output.tree = FALSE, scenarios = "S3", r = 1)

# outputs appear in R terminal - including unmatched names
# these must be corrected and the process repeated

phy <- sptree$scenario.3
is.ultrametric.phylo(phy)

library(phytools)
# get list of genera
tips <- phy$tip.label
genera <- unique(sapply(strsplit(tips,"_"),function(x) x[1]))
# drop all but one of each
ii <- sapply(genera, function(x,y) grep(x,y)[1], y=tips)
tree <- drop.tip(phy, setdiff(phy$tip.label, tips[ii]))
# check tree works
plotTree(tree, ftype="i")
# rename tips with the genus names
tree$tip.label<-sapply(strsplit(tree$tip.label,"_"),function(x) x[1])
length(unique(tree$tip.label))
writeNexus(tree, file="genus_tree.nex")