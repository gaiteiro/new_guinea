# setwd(".")
require(INLA)
require(rgdal)
require(rgeos)
library(ggplot2)
library(sp)
library(raster)
library(spdep)
library(spatialEco)
library(broom)
library(lme4)
library(caret)
library(spmoran)
library(tidyr)
library(dplyr)
library(ncf)

# Moran's I
grid3 <- readOGR("grid_10km_with_vars.shp")
ww <-  nb2listw(lattice_temp, style='B', zero.policy=TRUE)
moran(grid3@data$def_count, ww, n=length(ww$neighbours), S0=Szero(ww), zero.policy=TRUE)
moran.test(grid3@data$def_count, ww, randomisation=FALSE, zero.policy=TRUE)
moran.mc(grid3@data$def_count, ww, nsim=999, zero.policy=TRUE)
wm <- nb2mat(lattice_temp, style='B', zero.policy=TRUE) # spatial weights matrix
rwm <- mat2listw(wm, style='W')
y <- grid3@data$def_count
ybar <- mean(y)
n <- length(grid3@data$def_count)
ms <- cbind(id=rep(1:n, each=n), y=rep(y, each=n), value=as.vector(wm * y))
ms <- ms[ms[,3] > 0, ]
ams <- aggregate(ms[,2:3], list(ms[,1]), FUN=mean)
ams <- ams[,-1]
colnames(ams) <- c('y', 'spatially lagged y')
plot(ams)
# Moran's scatterplot
reg <- lm(ams[,2] ~ ams[,1])
abline(reg, lwd=2)
abline(h=mean(ams[,2]), lt=2)
abline(v=ybar, lt=2)
local_moran <- localmoran(grid3$def_count, ww)
grid3@data$lmr <- local_moran[, 'Ii'] # local observed Moran's I
summary(grid3@data$lmr)
plot(grid3@data$lmr, cex=0.6, pch=20)
histogram(grid3@data$lmr)
spplot(grid3, "lmr")
# mgimond.github.io/Spatial
coo <- coordinates (grid3)
sdist <- dnearneigh(coo, 0, 50000) # all polygons within 50km are neighbours
lw <- nb2listw(sdist, style="W", zero.policy=T)
MI <- moran.mc(grid3@data$def_count, lw, nsim=999, zero.policy=TRUE)
# 0.18743
plot(MI, main="", las=1)

# investigate spatial autocorrelation
nb.gab <- spdep::graph2nb(spdep::gabrielneigh(coords, nnmult=10), sym=TRUE)
par(mar=c(0,0,0,0))
plot(nb.gab, coords)
listw.gab <- spdep::nb2listw(nb.gab)
spdep::moran.test(grid10$c20, listw.gab)
spdep::moran.test(grid10$cost_mean, listw.gab)
spdep::moran.test(grid10$elev_mean, listw.gab)
spdep::moran.test(grid10$PAs_mean, listw.gab)

# correlogram
m1 <- meigen(cmat = nb2mat(nb.gab, style='W'))
data_sf <- st_as_sf(grid2)
queen <- dnearneigh(coords, d1=0, d2=sqrt(2)*10100)
card_queen <- card(queen)
grid2@data$card_queen <- card_queen
queenk <- knearneigh(coords, k=8)
queenw <- nb2listw(knn2nb(queenk),style='W')
grid2@data$x <- coords[,1]
grid2@data$y <- coords[,2]
cor1 <- with(grid2@data, correlog(x, y, c20, increment = 10000, resamp = 0))
saveRDS(cor1, file = "10km_correlogram.rds")
plot(cor1)

# 10km grid with deforestation counts added

# load shapefile with covariates for regression
grid2 <- readOGR("grid_10km_with_vars.shp")
# files exported from QGIS - calculations in Excel, 10066 rows = grid squares
# origin is Global Forest Watch / GLAD dataset
events <- read.csv('10Kgrid_def_calcs.csv', header=T)
grid2@data$pop00 <- events$pop00
grid2@data$pop01 <- events$pop01
grid2@data$pop02 <- events$pop02
grid2@data$pop03 <- events$pop03
grid2@data$pop04 <- events$pop04    
grid2@data$pop05 <- events$pop05    
grid2@data$pop06 <- events$pop06    
grid2@data$pop07 <- events$pop07    
grid2@data$pop08 <- events$pop08    
grid2@data$pop09 <- events$pop09    
grid2@data$pop10 <- events$pop10    
grid2@data$pop11 <- events$pop11    
grid2@data$pop12 <- events$pop12   
grid2@data$pop13 <- events$pop13   
grid2@data$pop14 <- events$pop14   
grid2@data$pop15 <- events$pop15   
grid2@data$pop16 <- events$pop16   
grid2@data$pop17 <- events$pop17   
grid2@data$pop18 <- events$pop18
grid2@data$pop19 <- events$pop19
grid2@data$pop20 <- events$pop20 
# cumulative deforestation
# deforestation 2001-20: start with column 02 = deforestation to 2002 (2001 is zero)
grid2@data$c02 <- events$c02
grid2@data$c03 <- events$c03 
grid2@data$c04 <- events$c04 
grid2@data$c05 <- events$c05 
grid2@data$c06 <- events$c06 
grid2@data$c07 <- events$c07 
grid2@data$c08 <- events$c08 
grid2@data$c09 <- events$c09 
grid2@data$c10 <- events$c10 
grid2@data$c11 <- events$c11 
grid2@data$c12 <- events$c12 
grid2@data$c13 <- events$c13 
grid2@data$c14 <- events$c14 
grid2@data$c15 <- events$c15 
grid2@data$c16 <- events$c16 
grid2@data$c17 <- events$c17 
grid2@data$c18 <- events$c18 
grid2@data$c19 <- events$c19 
grid2@data$c20 <- events$c20 
# check for NAs
sum(is.na(grid2@data$elev_mean))
sum(is.na(grid2@data$cost_mean))
sum(is.na(grid2@data$PAs_mean))
sum(is.na(grid2@data$cop00_majo))
sum(is.na(grid2@data$pop00_mean))
# merge some vegetation categories with few examples
grid2@data$cop00 <- ifelse(grid2@data$cop00_majo == 0, NA, grid2@data$cop00_majo)
grid2@data$cop00 <- ifelse(grid2@data$cop00_majo == 11, 10, grid2@data$cop00)
grid2@data$cop00 <- ifelse(grid2@data$cop00_majo == 12, 10, grid2@data$cop00)
grid2@data$cop00 <- ifelse(grid2@data$cop00_majo == 20, 10, grid2@data$cop00)
grid2@data$cop00 <- ifelse(grid2@data$cop00_majo == 80, 50, grid2@data$cop00)
grid2@data$cop00 <- ifelse(grid2@data$cop00_majo == 100, 110, grid2@data$cop00)
grid2@data$cop00 <- ifelse(grid2@data$cop00_majo == 130, 110, grid2@data$cop00)
grid2@data$cop00 <- ifelse(grid2@data$cop00_majo == 150, 110, grid2@data$cop00)
# remove NAs
grid2 <- sp.na.omit(grid2, col.name = "cop00_majo") # there must be at least 1 NA or this will return grid2=null!
grid2 <- sp.na.omit(grid2, col.name = "cost_mean") # ditto - so check NAs first
# check again
sum(is.na(grid2@data$elev_mean))
sum(is.na(grid2@data$PAs_mean))
sum(is.na(grid2@data$pop00_mean))
# scale elevation and cost
grid2@data[c("elev_mean", "cost_mean")] <- scale(grid2@data[c("elev_mean", "cost_mean")])
grid2 <- grid2[grid2@data$cop00 != "190",] # not enough cop = 190 data (= urban)
grid2@data$cop00 <- as.factor(grid2@data$cop00)
grid2@data$id2 <- 1:nrow(grid2@data)
saveRDS(grid10, file = "inla_shp_grid10km.rds")
grid10 <- grid2@data
saveRDS(grid10, file = "inla_cml_grid10km.rds")

# tibble for cumulative counts for use in regression - using c02-c20, ignoring pop00-20
grid10b <- grid10[,44:85]
grid10b <- cbind(grid10b, grid10[,(c(5, 6, 7, 27))])
colnames(grid10b)[43] <- ("meancost") # so column title doesn't start with 'c'
colnames(grid10b)[41] <- ("veg")
grid10b <- grid10b %>%
  tidyr::pivot_longer(cols = starts_with('c'), values_to = 'def') %>% # should just apply to c02-c20 
  group_by(id2) %>%
  mutate(tstart = row_number(),
         tstop = tstart+1) %>%
  select(-name)
# result is tibble with 185079 rows 
# equals 9741 grid squares times 1 row for each of 19 years
# year 1 = change from 2001 to 2002 i.e. starting value
# year 19 = change from 2019 to 2020
# tibble is sorted by year, then by grid square ID

# same but with population, if this variable is wanted for regression
grid10b <- grid10[,64:87]
grid10b <- cbind(grid10b, grid10[,(c(5, 6, 27))])
grid10c <- pivot_longer(grid10b, cols = -c('id2', 'cop00', 'elev_mean', 'PAs_mean', 'cost_mean'),
                        names_to = c(".value", "year"),
                        names_pattern = "(X|pop)(.*)",
                        names_transform = list(year = as.integer))

# regression with zero-inflated binomial model

# use grid10b saved earlier
grid10b$id3 <- grid10b$id2
grid10f <- grid10b[order(grid10b$tstart),] # make sure it's ordered by year
grid10f <- grid10f[,22:29]
grid2 <- readOGR(paste(getwd(), "/grid_10km_with_vars.shp", sep = ""))
 # prepare spatial lattice for INLA
lattice_temp <- poly2nb(grid2, row.names = grid2@data$id2)
nb2INLA(paste(getwd(), "/lattice.graph", sep = ""), lattice_temp)
lattice.adj <- paste(getwd(), "/lattice.graph", sep = "")

# family = zeroinflatednbinomial1 - needs lots of RAM

formula <- def2 ~ meancost + elev_mean + veg + PAs_mean +
  f(id2, model = "bym", graph = lattice.adj, scale.model=TRUE) + # spatial term
  f(id3, model = "iid", tstart) + tstart # year iid + year

inlanb2 <- inla(formula, data = grid10f, family = "zeroinflatednbinomial1", control.family=list(link='log'),
                control.predictor=list(link=1, compute=TRUE), 
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE), verbose=FALSE)
saveRDS(inlanb2, file = paste(getwd(), "/inla_zinbfull.rds", sep = ""))

# hold out 5 years' data
grid10f$def2 <- ifelse(grid10f$tstart > 14, NA, grid10f$def) # set last 5 years to NA
# in INLa, instead of fitting a model and then predicting with new data input,
# we fit the model with 'new data' response set to NA - 
# so fit and prediction are simultaneous
inlanb2 <- inla(formula, data = grid10f, family = "zeroinflatednbinomial1", control.family=list(link='log'),
                control.predictor=list(link=1, compute=TRUE), 
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE), verbose=FALSE)
saveRDS(inlanb2, file = paste(getwd(), "/inla_zinb5yr.rds", sep = ""))

# hold out 10 years' data
grid10f$def2 <- ifelse(grid10f$tstart > 9, NA, grid10f$def) # set last 10 years to NA
inlanb2 <- inla(formula, data = grid10f, family = "zeroinflatednbinomial1", control.family=list(link='log'),
                control.predictor=list(link=1, compute=TRUE), 
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE), verbose=FALSE)
saveRDS(inlanb2, file = paste(getwd(), "/inla_zinb10yr.rds", sep = ""))

# project for 2021-25

gridnew21 <- grid10c[grid10c$year == 20,] # same as previous years (except for outcome)
gridnew21$def2 <- NA # set outcome (deforestation) as NA - model will predict this
gridnew21$year <- 21 # change year variable
gridnew22 <- gridnew21
gridnew22$def2 <- NA
gridnew22$year <- 22
gridnew23 <- gridnew21
gridnew23$def2 <- NA
gridnew23$year <- 23
gridnew24 <- gridnew21
gridnew24$def2 <- NA
gridnew24$year <- 24
gridnew25 <- gridnew21
gridnew25$def2 <- NA
gridnew25$year <- 25
gridnew2125 <- rbind(gridnew21, gridnew22, gridnew23, gridnew24, gridnew25)
grid25test <- rbind(grid10c, gridnew2125)
inlanb2 <- inla(formula, data = grid10f, family = "zeroinflatednbinomial1", control.family=list(link='log'),
                control.predictor=list(link=1, compute=TRUE), 
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE), verbose=FALSE)
saveRDS(inlanb2, file = paste(getwd(), "/inla_zinbnew.rds", sep = ""))

# calculate increase in deforestation from 2015-20 (actual and predicted)
grid10f <- readRDS("inla_grid10f.rds")
inlanb2 <- readRDS( file = "inla_zinb.rds")
grid2 <- readOGR(paste(getwd(), "/grid_10km_with_vars.shp", sep = ""))
grid1620 <- cbind(grid10f, inlanb2$summary.fitted.values) # all years
grid1620b <- grid1620[, c("id2", "tstart", "def", "mean")]
grid19 <- grid1620b[grid1620b$tstart == 19,] # year 2020
grid14 <- grid1620b[grid1620b$tstart == 14,] # year 2015 (base is 0)
colnames(grid19) <- c("id2b", "tstartb", "def19", "mean19")
grid1419 <- cbind(grid14, grid19)
grid1419 <- mutate (grid1419, incact = def19-def, incpred = mean19-mean) # change in deforestation 2015-20
# add to SpatialPolygons
grid10 <- readRDS("inla_shp_grid10km.rds")
grid10$ID <- as.numeric(grid10$ID) # 9741 grid square IDs corresponding to model
grid2@data$ID <- as.numeric(grid2@data$ID) # 1066 grid square IDs
idlist <- as.data.frame(grid2@data$ID)
idlist2 <- as.data.frame(grid10$ID)
gridpred <- grid1419[, c('incact', 'incpred')]
gridpred <- cbind(gridpred, idlist2) # add back original IDs before grids with no predictor values were excluded
colnames(gridpred)[3] <- "ID"
colnames(idlist) <- "ID"
idjoin <- full_join(idlist, gridpred, by="ID") # 10066 items including NAs
grid2@data$actual <- as.numeric(idjoin$incact)
grid2@data$preds <- idjoin$incpred
writeOGR(obj=grid2, dsn="grid_10km_inla_1419.shp", layer="1", driver="ESRI Shapefile")
cor(gridpred$incact,gridpred$incpred, method="spearman") # 0.6555299
cor(gridpred$incact,gridpred$incpred, method="kendall") # 0.4868637

# predicted vs actual deforestation over 2010-2020
grid10f <- readRDS("grid10f_09.rds")
inlanb2 <- readRDS( file = "inla_zinb10yr.rds")
grid2 <- readOGR(paste(getwd(), "/grid_10km_inla_1419.shp", sep = ""))
grid1620 <- cbind(grid10f, inlanb2$summary.fitted.values) # all years
grid1620b <- grid1620[, c("id2", "tstart", "def", "mean")]
grid19 <- grid1620b[grid1620b$tstart == 19,] # year 2020
grid09 <- grid1620b[grid1620b$tstart == 9,] # year 2010 (base is 0)
colnames(grid19) <- c("id2b", "tstartb", "def19", "mean19")
grid0919 <- cbind(grid09, grid19)
grid0919 <- mutate (grid0919, incact = def19-def, incpred = mean19-mean) # change in deforestation 2010-20
gridpred <- grid0919[, c('incact', 'incpred')]
cor(gridpred$incact,gridpred$incpred, method="spearman") # 0.7070741
cor(gridpred$incact,gridpred$incpred, method="kendall") # 0.5315181

# predictions to 2025
grid10f <- readRDS("grid10f_new.rds")
inlanb2 <- readRDS( file = "inla_zinbnew.rds")
grid2 <- readOGR(paste(getwd(), "/grid_10km_inla_1419.shp", sep = ""))
grid1620 <- cbind(grid10f, inlanb2$summary.fitted.values) # all years
grid1620b <- grid1620[, c("id2", "tstart", "def", "mean")]
grid24 <- grid1620b[grid1620b$tstart == 24,] # year 2025
grid19 <- grid1620b[grid1620b$tstart == 19,] # year 2020
colnames(grid24) <- c("id2b", "tstartb", "def24", "mean24")
grid1924 <- cbind(grid19, grid24)
grid1924 <- mutate (grid1924, incact = def24-def, incpred = mean24-mean) # change in deforestation 2020-25
# add to SpatialPolygons shapefile so we can plot it in QGIS
grid10 <- readRDS("inla_shp_grid10km.rds")
grid10$ID <- as.numeric(grid10$ID) # 9741 grid square IDs corresponding to model
grid2@data$ID <- as.numeric(grid2@data$ID) # 1066 grid square IDs
idlist <- as.data.frame(grid2@data$ID)
idlist2 <- as.data.frame(grid10$ID)
gridpred <- grid1924[, c('incact', 'incpred')]
gridpred <- cbind(gridpred, idlist2) # add back original IDs before grids with no predictor values were excluded
colnames(gridpred)[3] <- "ID"
colnames(idlist) <- "ID"
idjoin <- full_join(idlist, gridpred, by="ID") # 10066 items including NAs
grid2@data$actual <- as.numeric(idjoin$incact)
grid2@data$preds <- idjoin$incpred
writeOGR(obj=grid2, dsn="grid_10km_inla_5yr.shp", layer="1", driver="ESRI Shapefile")

