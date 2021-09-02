rm(list=ls()) 

library(sf)
library(raster)
library(RQGIS)
library(mlr)
library(dplyr)
library(RStoolbox)
library(stringr)
library(rgeos)
library(INLA)
library(spdep)
library(sp)
library(rgdal)

# stack variables
REPO_HOME <- paste0(pth, '/NG_rasters/out/')
files <- list.files(REPO_HOME, pattern=".tif")
outfiles <- paste0(REPO_HOME, files)
stack <- stack(outfiles)
writeRaster(stack, 'env_variables_46_plus_DEM.tif', format='GTiff')
write.csv(files,"variable_list.csv",row.names = F) # keep track of which rasters in which order!

# convert lat/long in degrees/minutes/seconds  to decimal
pngcoords <- read.csv("latlong.csv", header = T)
pngcoords$Latitude = measurements::conv_unit(pngcoords$Latitude, from = 'deg_dec_min', to = 'dec_deg')
pngcoords$Longitude = measurements::conv_unit(pngcoords$Longitude, from = 'deg_dec_min', to = 'dec_deg')
pngcoords2 <- cbind(pngcoords$Longitude, pngcoords$Latitude)

# Moran's I calculations: environmental variables
data6a <- read.csv(paste0(getwd(), "/25kmgrids_1820.csv"), header = T)
data6 <- st_as_sf(data6a, coords=c("latitude", "longitude"))
cellsize<-25000
rook <- dnearneigh(data6, d1=0, d2=cellsize + 2000)
data6$card_rook <- card(rook)
data6 <- subset(data6, card_rook > 0) # remove cells with no neighbours
# now must recalculate
rook <- dnearneigh(data6, d1=0, d2=cellsize + 2000)
data6$card_rook <- card(rook)
# calculate weights
rookw <- nb2listw(rook, style='W')
moranlist <- list()
data7 <- st_drop_geometry(data6)
for (i in 2:48){
  x3 <- moran.test(data7[i], zero.policy=TRUE, na.action=na.exclude, rookw)
  moranlist[i-1] <- x[3]
}
saveRDS(moranlist, file = "1820_grids_env_vars_moranI.rds")

# INLA preliminary tests to investigate most influential environmental variables
# a spatial regression with response as simply the number of records (of any species) per grid square

# read in data as shapefile (i.e. having used QGIS to check boundaries, extents etc.)
lattice <- readOGR("grid_with_points.shp")
data <- lattice@data
data$ID <- as.numeric(data$ID)
Lattice_Temp <- poly2nb(lattice, row.names = data$ID) # construct the neighbour list
# create the adjacency matrix in INLA format
nb2INLA("Lattice.graph", Lattice_Temp) 
Lattice.adj <- "Lattice.graph" # name object and save

# spatial effect: ID is a numeric identifier for each area in the lattice
# bym not CAR or besag for area data
# start off with all variables
formula <- NUMPOINTS ~ NGcost_mea + elev_mean + elev_min + elev_max + 
  X1_mean + X2_mean + X3_mean + X4_mean + X5_mean + X6_mean + X7_mean + 
  X8_mean + X9_mean + X10_mean + X11_mean + X12_mean + X13_mean + X14_mean + 
  X15_mean + X16_mean + X17_mean + X18_mean + X19_mean + X21_mean + X22_mean + 
  X23_mean + X24_mean + X25_mean + X26_mean + X27_mean + X28_mean + X29_mean + 
  X30_mean + X31_mean + X32_mean + X33_mean + X34_mean + X35_mean + X36_mean + 
  X37_mean + X38_mean + X39_mean + X40_mean + X41_mean + X42_mean + X43_mean + 
  X44_mean + X45_mean + X46_mean + f(ID, model = "bym", graph = Lattice.adj)
Mod_Lattice <- inla(formula,     
                    family = "poisson", 
                    data = data,
                    control.compute = list(mlik = T, cpo = T, dic = T, waic = T)) 

# try zero-inflated: no improvement
Mod_Lattice2 <- inla(formula,     
                     family = "zeroinflatedpoisson1", 
                     data = data,
                     control.compute = list(mlik = T, cpo = T, dic = T, waic = T))         

# set priors on hyperparameters - log gamma since it's Poisson
formula_p <- NUMPOINTS ~ NGcost_mea + elev_mean + elev_min + elev_max + 
  X1_mean + X2_mean + X3_mean + X4_mean + X5_mean + X6_mean + X7_mean + 
  X8_mean + X9_mean + X10_mean + X11_mean + X12_mean + X13_mean + X14_mean + 
  X15_mean + X16_mean + X17_mean + X18_mean + X19_mean + X21_mean + X22_mean + 
  X23_mean + X24_mean + X25_mean + X26_mean + X27_mean + X28_mean + X29_mean + 
  X30_mean + X31_mean + X32_mean + X33_mean + X34_mean + X35_mean + X36_mean + 
  X37_mean + X38_mean + X39_mean + X40_mean + X41_mean + X42_mean + X43_mean + 
  X44_mean + X45_mean + X46_mean + f(ID, model = "bym", graph = Lattice.adj, scale.model = TRUE,
                                     hyper = list(
                                       prec.unstruct = list(prior = "loggamma", param = c(1,0.001)),   # precision for the unstructured effect (residual noise)
                                       prec.spatial =  list(prior = "loggamma", param = c(1,0.001))    # precision for the spatial structured effect
                                     ))
Mod_Lattice_p <- inla(formula_p,
                      family = "poisson",
                      data = data,
                      control.compute = list(cpo = T)
)

# I will skip the backwards selection steps
# the least significant variables (by WAIC) were dropped one by one
# final formula:
formula_p4 <- NUMPOINTS ~ NGcost_mea + elev_max + 
  X6_mean  + X11_mean + X12_mean + X15_mean + X17_mean + X19_mean + 
  X38_mean + X43_mean + X44_mean + X45_mean + X46_mean + 
  f(ID, model = "bym", graph = Lattice.adj, scale.model = TRUE,
    hyper = list(
      prec.unstruct = list(prior = "loggamma", param = c(1,0.001)),  
      prec.spatial =  list(prior = "loggamma", param = c(1,0.001))   
    ))
Mod_Lattice_p4 <- inla(formula_p4,
                       family = "poisson",
                       data = data,
                       control.compute = list(cpo = T)
)

# Hmsc models using cost + 11 environmental variables
# 716 25km grid sites, spp with records in >=30 grids - analysis in QGIS:
# QGIS/Vector analysis/Count Points in Polygon (use unique class field so only one species counted in each square)
# then export to csv: QGIS/MMQGIS plugin/Attributes Export to CSV file
# and use Excel to eliminate spp with records in <30 grids
# Same procedure was followed to get map of 1156 genera
y <- read.csv('/species_by_site_716_pa.csv', header=T)
XData <- read.csv('env_vars_by_site_716.csv', header=T)
Y=as.matrix(y)
library(Hmsc)
# set MCMC variables
nChains = 2
thin = 5
samples = 1000
transient = 500*thin
verbose = 50*thin
n <- nrow(XData)
studyDesign = data.frame(sample = as.factor(1:n))
xycoords <- XData[c("xcoord", "ycoord")]
rL.spatial = HmscRandomLevel(sData = xycoords)
rL.spatial = setPriors(rL.spatial,nfMin=1,nfMax=1)

m = Hmsc(Y=Y, XData=XData, XFormula=~NG_costmean + elev_range + cecmean + t10mean + 
           thrmmea + trimean + fparmea + laimean + soc1mea + sq3mean + sq5mean + bio7mea,
         studyDesign=studyDesign, ranLevels=list("sample"=rL.spatial), distr="probit")
# this model will be used for further backwards selection of variables 
# for final Hmsc model of individual environmental predictors

# extract PCA variables

pth <- getwd()
REPO_HOME <- file.path(pth, "rasters")
files <- list.files(REPO_HOME, pattern=".tif")
outfiles <- paste0(REPO_HOME, files)
stack <- stack(outfiles)
library(RStoolbox)
rpca <- rasterPCA(stack, nComp = 15, spca = TRUE, maskCheck = TRUE)
pcamap <- rpca$map
saveRDS(pcamap, file = "rpca_46_plus_DEM.rds")

# copernicus vegetation type data needs to be extracted from ncdf4 format

library(ncdf4)
pth <- getwd()
nc_data <- nc_open('ESACCI-LC-L4-LCCS-Map-300m-P1Y-2000-v2.0.7cds.nc')
r2 <- raster(pth, 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2000-v2.0.7cds.nc', 
             varname="lccs_class")
writeRaster(r2, 'copernicus15.tif', format="GTiff")
nc_close(nc_data)
cop <- raster(paste(getwd(), "/copernicus_UTM_1000m_res.tif", sep = ""))
m <- aggregate(cop, fact = 3, fun = modal, na.rm = TRUE)

# extract deforestation values 2000-20, by 10km grid square

library(exactextractr) # works fast with large files
# first compress huge .tif file
def <- raster(paste(getwd(), "/2000_2020_UTM_samepixels_cubic.tif", sep = ""))
tifoptions <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")
writeRaster(def, "2000_2020_UTM_samepixels_cubic_compressed.tif",
                          options = tifoptions, overwrite = FALSE)
def <- raster(paste(getwd(), "/2000_2020_UTM_samepixels_cubic_compressed.tif", sep = ""))
grid <- readOGR(paste(getwd(), "/grid_10km.shp", sep = ""))
results <- exact_extract(def, grid, function(value, count) {table(value)})
write.csv(results,"10Kgrid_def_calcs.csv")

# extract mean population values by 10km grid square

grid50 <- readOGR("grid_10km.shp")
lenval <- length(grid50@data$left)
lenval <- length(grid10@data$left)
popnresults <- data.frame(matrix(data=NA, nrow = length(seq(1:lenval)), 
                             ncol = length(seq(1:21))))
folder <- list.files(pth, "WorldPop2020/merged")
for (i in 1:(length(folder))){
  rp <- raster(paste(getwd(), "/", folder[i], sep = ""))
  message(i)
  z1 <- exact_extract(rp, grid2, fun="mean")
  popnresults[,i] <- z1
}
write.csv(popnresults,"10Kgrid_pop_calcs.csv") # this was added to 10Kgrid_def_calcs.csv

# test significance of 27K collection records against cost-distance raster

cost <- raster(paste(getwd(), "/NG_cost_distance_map.tif", sep = ""))
points <- readOGR("27K_envs_clipped_tobuffer_UTM_singlepoint.shp") 
# had to convert the above in QGIS from multipoint to singlepoint
# using: QGIS/SAGA/Vector tools/Convert multipoints to points
( rmean <- cellStats(cost, stat="mean", na.rm=TRUE, asSample=FALSE) )  
( rvar <- cellStats(cost, stat="sd", na.rm=TRUE, asSample=FALSE)^2 )     
( rquant <- quantile(cost, probs = c(0.25, 0.50, 0.75)) )
extracted <- raster::extract(cost, points, method='bilinear')
( smean <- mean(extracted, na.rm=TRUE) )
( svar <- var(extracted, na.rm=TRUE) )
( squant <- quantile(extracted, probs = c(0.25, 0.50, 0.75), na.rm=TRUE) )
e <- na.omit(extracted)
r <- as.vector(as.matrix(cost))
r <- na.omit(r)
t.test(e, y = r)
# independent 2-group Mann-Whitney U Test is better - non-parametric
wilcox.test(r,e)

