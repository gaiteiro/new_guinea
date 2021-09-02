# prepare data for BIODIVERSE

library(tidyr)
library(dplyr)
coords <- read.csv('50Kgridcoords511.csv', header=T)
pred_mean_nocost <- readRDS("pred2costzeromean.rds") # matrix of mean predicted probabilities
preds_nc <- as.data.frame(pred_mean_nocost)
pred_sd_nocost <- readRDS("pred2costzerosd.rds") # matrix of predictions SDs
pred_CIs <- pred_mean_nocost + (pred_sd_nocost * 2)
preds_ci <- as.data.frame(pred_CIs)
preds_ci$ID <- 1:nrow(preds_ci) # keep track of rows
# get rid of records whose probability +2 standard deviations does not exceed 0.5
ci2 <- tidyr::pivot_longer(data = preds_ci, cols = !ID, values_to = 'count')
ci3 <- ci2[ci2$count > 0.5,]
ci4 <- ci3 %>% unite(ID, name, col=uniqueID, sep="_") # merge ID and genus name into 1 column
nocost <- cbind(coords, preds_nc)
# get a long (many rows) list of each genus / each site
nocost2 <- tidyr::pivot_longer(data = nocost, cols = !(xcoords|ycoords|ID), values_to = 'count')
nocost3 <- nocost2 %>% unite(ID, name, col=uniqueID, sep="_")
nocost4 <- cbind(nocost2, nocost3$uniqueID)
# not every genus is present in every grid square
# remove all items whose grid ID + genus name is not found in ci4:
nocost5 <- inner_join(nocost4, ci4, by = c("nocost3$uniqueID" = "uniqueID"))
# select columns needed for BIODIVERSE
nocost5 <- nocost5[, c("ID", "xcoords", "ycoords", "name", "count.x")]
length(unique(nocost5$ID))
# add coordinates of all 549 grid squares just in case there are blank spaces in the BIODIVERSE map
extras <- read.csv('50Kgridcoords.csv', header=T)
names(extras) <- c("xcoords", "ycoords")
extras$name <- "x"
extras$count.x <- 0
extras$ID <- 1:nrow(extras)
nocost5 <- rbind(nocost5, extras)
write.table(nocost5, 
            file = file.path(getwd(), "biodiverse_input.csv"), 
            sep = ",", row.names = FALSE, col.names = TRUE)

# also for BIODIVERSE prepare actual values in 293 grids
x <- read.csv('50Kgrid_4var_ecoreg_cop99_PCA_nozeros_293.csv', header=T) 
y <- read.csv('Ygrid50_293_1155.csv', header=T)
str(x)
coords <- x[,c("xcoord", "ycoord")]
act1 <- cbind(coords, y)
act2 <- tidyr::pivot_longer(data = act1, cols = !("xcoord"|"ycoord"), values_to = 'count')
act2a <- act2[act2$count >0,]
# add coordinates of all 549 grids so there are no spaces in the Biodiverse map
extras <- read.csv('C:/Users/smith/Documents/nnMaxEnt3/50Kgridcoords.csv', header=T)
extras$name <- "x"
extras$count <- 0
act3 <- rbind(act2a, extras)
write.table(act3, 
            file = file.path(getwd(), "actuals_for_biodiverse.csv"), 
            sep = ",", row.names = FALSE, col.names = TRUE)

# load this into BIODIVERSE along with .nex phylogeny; build spatial index
# run spatial analysis: subroutine calc_pe_single (checkbox 'Phylogenetic Endemism single'):
# and export as .csv
# see https://github.com/shawnlaffan/biodiverse/wiki/Indices#phylogenetic-endemism-indices