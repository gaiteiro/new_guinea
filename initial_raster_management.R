library(raster)
library(rgeos)
library(sp)
library(sf)
library(stringr)
library(vegan)

# setwd(".")
REPO_HOME <- paste(getwd(), "/NG_rasters/in/", sep = "") # set up source directory
REPO_DEST <- paste(getwd(), "/NG_rasters/out/", sep = "")
REPO_SHP <- paste(getwd(), "/NG_rasters/shpfiles/", sep = "")
BOUNDARY <- paste(getwd(), "/NG_rasters/NG_boundary_UTM_10km_buffer.shp", sep = "")
BOUNDARY2 <- paste(getwd(), "/NG_rasters/NG_boundary_UTM_10km_buffer_1000sq.tif", sep = "")
res <- 1000

# rasterize the UTM boundary shapefile
boundary <- st_read(BOUNDARY)
ext <- extent(boundary)
emptyr <- raster(ext, res=res)
boundaryr <- rasterize(boundary, emptyr)
crs(boundaryr) <- 32554
UTM_extent <- extent(boundaryr)
emptyrUTM <- raster(UTM_extent, res=res)
writeRaster(boundaryr, BOUNDARY2, format="GTiff")
# will use this lat/long extent
crop_extent <- extent(128.6, 157, -12.7, 2.1)

# list shapefiles to rasterize
filesshp <- list.files(REPO_SHP, pattern=".shp")
# set the output directory
outfilesshp <- paste0(REPO_HOME, filesshp)

for(i in 1:length(filesshp)) {
  mypolygon <- st_read(paste(REPO_SHP, filesshp[i], sep = ""))
  mypolygonUTM <- st_transform(mypolygon, 32554)
  message(extent(mypolygonUTM))
  mypolygonUTM <- gBuffer(mypolygonUTM, byid=TRUE, width=0)
  mypolygoncropped <- st_crop(mypolygonUTM, UTM_extent)
  message(extent(mypolygoncropped))
  rasterizedpolygon <- rasterize(mypolygoncropped, emptyrUTM)
  finalraster = resample(rasterizedpolygon, boundaryr, "bilinear")
  outname <- str_sub(outfilesshp[i],1, nchar(outfilesshp[i])-4)
  writeRaster(finalraster, outname, format="GTiff")
}

# find tif rasters
files <- list.files(REPO_HOME, pattern=".tif")
# set the output directory
outfiles <- paste0(REPO_DEST, files)

# NG boundary raster used as standard for resampling
for(i in 1:length(files)) {
  r <- raster(paste(REPO_HOME, files[i], sep = ""))
  message(files[i])
  cropped <- crop(r, crop_extent) # crop to 1 degree wider than boundary
  message(extent(cropped))
  rUTM <- projectRaster(cropped, crs= 32554) # convert raster into UTM 54S
  message(extent(rUTM))
  rUTM1000m = raster::resample(rUTM, boundaryr, "bilinear") # resample to same grid as in UTM boundary tif
  # note: if categorial should be 'nbg' - nearest neighbour - instead of 'bilinear'
  message(extent(rUTM1000m))
  clipped <- mask(rUTM1000m, boundaryr) # mask to remove very small islands
  message(extent(clipped))
  message(outfiles[i])
  writeRaster(clipped, outfiles[i], format="GTiff")
}

# stack and save results
raster_stack = stack(outfiles)
writeRaster(raster_stack,"raster_stack.tif", format="GTiff")
message(outfiles)rp = dplyr::select(rp, -id, -spri) %>%
  st_drop_geometry()

# alternatively
#files <- list.files(REPO_DEST, pattern=".tif")
#outfiles <- paste0(REPO_DEST, files)
#raster_stack = stack(outfiles)

# convert FAO soil files from asc
# NG boundary raster used as standard for resampling
boundaryr <- raster(paste(getwd(), "/NG_rasters/NG_boundary_UTM_10km_buffer_1000sq.tif", sep = ""))
FAO <- paste(getwd(), "/NG_rasters/shpfiles/FAO_soils/", sep = "") 
files <- list.files(FAO, pattern=".asc")
for(i in 1:length(files)) {
  r <- raster(paste(FAO, files[i], sep = ""))
  crs(r) <- 4326
  cropped <- crop(r, crop_extent) # crop to 1 degree wider than boundary
  rUTM <- projectRaster(cropped, crs= 32554) # convert raster into UTM 54S
  message(extent(rUTM))
  rUTM1000m = raster::resample(rUTM, boundaryr, "bilinear") # resample to same grid as in UTM boundary tif
  # note: if categorial should be 'nbg' - nearest neighbour - instead of 'bilinear'
  message(extent(rUTM1000m))
  clipped <- mask(rUTM1000m, boundaryr) # mask to remove very small islands
  message(extent(clipped))
  outname <- str_sub(outfiles[i],1, nchar(outfiles[i])-4)
  outname <- paste0(outname, ".tif")
  message(outfiles[i])
  writeRaster(clipped, outname, format="GTiff")
}


