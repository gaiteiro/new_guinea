# new_guinea
Biodiversity and Conservation Priority Setting: Vascular Plants of New Guinea

Workflow

PNG_scrape.py

- I had a list of the genera in the pngplants.org website (Lae herbarium website) but the records, which are searchable by genus, had to be called up and downloaded with this script

initial_raster_management.R

- Requires  NG boundary shapefile .This boundary excludes a myriad of tiny islands, for which environmental data is often lacking, as well as Bougainville and a number of north-westerly islands belonging to non-Papuan Indonesian provinces

- This and other shapefiles (e.g. land systems, ecoregions) are rasterized so as to check extent

- Many files already available as geotiff (but often in a WGS84 projection, or of global extent)

- A boundary raster of New Guinea is used as model for resampling to the same extent and to UTM, 1km resolution

- Rasters are stacked and saved

- FAO soil files must be converted from ASC format

- Copernicus vegetation type data needs to be extracted from ncdf4 format

- A 10km grid square map of New Guinea is used as a basis for extracting population and deforestation data


data_preliminaries_various.R

-Some lines of code are provided to convert lat/long coordinates from degrees/min/sec into decimal

-The GBIF download was converted into Excel for processing, and the Lae data added (both have unique ID per record). Use of Excel allowed rapid sorting, manipulation of pivot tables and visual inspection of results. The data were also processed through QGIS to allow grid square ID and coordinates to be added. 

-Moran’s I calculations: these were performed on 25km2 grid squares – the code for extraction of environmental variables by grid square is also given

-INLA preliminary tests to investigate most influential environmental variables – a spatial regression with response as simply the number of records (of any species) per grid square. The same model was run in HMSC. 

-QGIS was used to prepare the data; at this stage I was using species rather than genera

-Extraction of species and genera per grid square followed the same procedure; joining attributes of species records (as a set of points), with grid squares (as a set of polygons) – tool used was QGIS/Vector analysis/Count Points in Polygon (using a unique class field so only one count of each species/genus was made in each grid square). The results were exported to csv: QGIS/MMQGIS plugin/Attributes Export to CSV file – and Excel was used to eliminate e.g. species with records in <30 grids

-PCA extraction made from (continuous, scaled) environmental rasters

-The cost-distance values for 27K collection records are compared with the overall values of the cost-distance raster

norberg_SDM_bakeoff.R

-This is not intended to be run all at once. It has been altered so as to perform only the SDM tests which were successful when applied to the dataset of 1,156 genera. Other SDMs are commented out. The relevant lists ‘mod_names’ and ‘pred_names’ have also had the unused SDMs commented out, and the file feat.csv has had relevant rows removed

-Various other scripts in the package would have to be altered slightly to adapt to my procedure

-The two Matlab scripts will need to be run in Matlab at the same time as model fitting is performed. The Norberg paper predates the release of HMSC in R, so HMSC must be tested in Matlab

-As well as the dataset of 1,156 genera, I also tested a dataset of 180 more common species (recorded over 716 25km2 grid square sites)

biodiversity_Hmsc.R

-This is the HMSC biodiversity modelling stage

-Input data came from QGIS analysis – comments include details of the QGIS tools used

-Note that 299 grid squares had PCA variables available, but only 293 grid squares had environmental variables available – this is because the second model includes ecoregion (category variable not used for PCA) and a few grid squares have missing ecoregion data

-Various models were tested (e.g. with and without the spatial term) with 2-fold cross-validation, and where practicable with 5-fold cross-validation

-2-fold cross-validation can be performed here – but note that the results may vary according to data split. For 5-fold cross-validation, code is give here to set up the 5 folds and to analyse the results after Matlab processing

-For a ‘single-species stacked ensemble model’ (with genera not species) I used Matlab, and 2-fold cross-validation only – code is give here to set up the 2 folds and to analyse the results after Matlab processing

fit_predict_hmsc_matlab_5cv.m

-Matlab code for 5-fold cross-validation using HMSC

hmsc_single_species_ensemble_predictions.m

-Matlab code for 2-fold cross-validation using HMSC

phylogeny.R

-Builds phylogeny and reduces species tree to genus tree

-The phylogeny is ultrametric, so whichever species is retained (and relabelled with genus name) the total tree length, and length from genus node to tip, is the same

-biodiverse_data_preparation.R

-Input data is a version of that used for the HMSC modelling with unnecessary columns removed

-BIODIVERSE requires an input file with 1 row per record per location, which is the output here

deforestation.R

-Local Moran’s I calculations for deforestation clustering – with correlogram etc.

-Deforestation and population variables are added to other variables and a dataframe is prepared for use in INLA regression model

-INLA regression: full model, models with 5yr and 10yr data withheld, and model with dummy inputs for 2021–25

-Calculations of predicted vs actual values

-Prepare data for plotting in QGIS (export to shapefile)

-QGIS was used to combine biodiversity and deforestation figures (Field calculator was used to take logarithm, normalize, and add the 2 datasets)
