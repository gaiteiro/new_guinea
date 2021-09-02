rm(list = ls(all = TRUE)); gc()

# fit and predict models with 5 PCA axes plus cost variable

# 1156 genera - first import collection records into QGIS as points (Add Layer/Add Delimited Text Layer), 
# Build polygons (50km grid squares) set to extent of NG boundary file (MMGIS plugin/Create Grid Layer)
# Clip grid polygons by NG boundary file (Geoprocessing tools/Clip) so not all are square - sea excluded
# Use Vector geometry/Centroids to get grid centroid coordinates
# Use Vector analysis/Count Points in Polygon to count genera (use unique class field so only 1 genus counted in each square)
# then export to csv: MMQGIS plugin/Attributes Export to CSV file
# NB: need an ID field to keep track of grid squares and their coordinates
X <- read.csv('50Kgrid_pca_511_sites.csv', header=T) # covariates for HMSC

# Use QGIS to extract grid square site ID of collection records (QGIS/Join attributes by location - use one-to-one)
# Export to Excel and use pivot table  to arrange genera (columns) and grid square sites (rows)
# Need 1 column for each genus: presence/absence 1/0, so convert all figures >=1 to 1
y <- read.csv('Y_all.csv', header=T) # response data for HMSC
# if presence records are not 1/0 use this:
# y = as.matrix(da[,-1])>0
# y = apply(y,MARGIN = 2,FUN = as.numeric)

# exclude site IDs not in both files
X1 <- X[X$ID %in% y$ID, , drop = FALSE]
y <- y[y$ID %in% X$ID, , drop = FALSE]
# drop ID column
y1 <- y[, -1]
# Hmsc requires a matrix as Y input (response)
Y1 <- as.matrix(y1)
X2 <- X[!X$ID %in% y$ID, , drop = FALSE] # keep excluded site IDs to later make predictions across ALL sites

library(Hmsc)
# set MCMC variables
nChains = 2
thin = 5
samples = 1000
transient = 500*thin
verbose = 50*thin

studyDesign <- data.frame(sample = as.factor(rownames(X1))) # keeps track of row IDs
xycoords <- X[c("xcoords", "ycoords")] # get coordinates of sites
rL = HmscRandomLevel(sData = xycoords) # spatial term
rL = setPriors(rL,nfMin=1,nfMax=10) # set min and max no. of latent variables
# set formula
fit1 = Hmsc(Y=Y1, XData=X1, XFormula= ~ env + env.1 + env.2 + env.3 + env.4 + cost,
            studyDesign=studyDesign, ranLevels=list("sample"=rL), distr="probit")
# run MCMC chains
fit1 = sampleMcmc(fit1, thin = thin, samples = samples, transient = transient,
                  nChains = nChains, verbose = verbose, updater=list(GammaEta=FALSE))
saveRDS(fit1, file = "fit1cost.rds")

# variance partitioning, convergence checks, evaluate model fit for internal consistency
VP = computeVariancePartitioning(fit1)
plotVariancePartitioning(fit1, VP=VP, las=2, horiz=F)
mpost <- convertToCodaObject(fit1)
effectiveSize(mpost$Beta) # effective size
gelman.diag(mpost$Beta)$psrf # PSRF
plot(mpost$Alpha[[1]])
summary(mpost$Alpha[[1]])$quantiles
preds = computePredictedValues(fit1)
MF = evaluateModelFit(hM=fit1, predY=preds3)
hist(MF$AUC, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$AUC),2))) # average AUC per site

# rebuild study design from grids across all New Guinea
X3 <- rbind(X1, X2)
studyDesign <- data.frame(sample = as.factor(rownames(X3)))
xycoords <- X3[c("xcoords", "ycoords")]
rL = HmscRandomLevel(sData = xycoords)
rL = setPriors(rL,nfMin=1,nfMax=10)
# make predictions from fitted model
pred1 = Hmsc:::predict.Hmsc(fit1, XData = data.frame(X3), type = "response", studyDesign = studyDesign, 
                            ranLevels=list("sample"=rL), expected = TRUE)
saveRDS(pred1, file = "pred1cost.rds")

pred_mean <- apply(simplify2array(pred1), 1:2, mean)
pred_sd <- apply(simplify2array(pred1), 1:2, sd)
saveRDS(pred_mean, file = "pred1costmean.rds")
saveRDS(pred_sd , file = "pred1costsd.rds")
mean_summed <- rowSums(pred_mean) # total number of genera per site
write.csv(mean_summed , "pcacost.csv", row.names=FALSE)

#### bias correction ####
# set cost variable to zero and re=predict from fitted model

X3$cost <- 0

studyDesign <- data.frame(sample = as.factor(rownames(X3)))
xycoords <- X3[c("xcoords", "ycoords")]
rL = HmscRandomLevel(sData = xycoords)
rL = setPriors(rL,nfMin=1,nfMax=10)

fit1 <- readRDS("fit2cost.rds")

pred1 = Hmsc:::predict.Hmsc(fit1, XData = data.frame(X3), type = "response", studyDesign = studyDesign, 
                            ranLevels=list("sample"=rL), expected = TRUE)
saveRDS(pred1, file = "pred2costzero.rds")

pred_mean <- apply(simplify2array(pred1), 1:2, mean)
pred_sd <- apply(simplify2array(pred1), 1:2, sd)
saveRDS(pred_mean, file = "pred2costzeromean.rds")
saveRDS(pred_sd , file = "pred2costzerosd.rds")
# export predictions as .csv
mean_summed <- rowSums(pred_mean) # total number of genera per site
write.csv(mean_summed , "pcacostnozeros.csv", row.names=FALSE) # graph this output in QGIS

# model with cost plus 4 environmental variables inc. ecoregion 
# I started off with 11 environmental variables - backwards selection by % of variance explained

X <- read.csv('50Kgrid_env_vars_293.csv', header=T) 
y <- read.csv('Ygrid50_293.csv', header=T)
Y <- as.matrix(y)
X$ecoreg <- as.factor(X$ecoreg)
library(Hmsc)
nChains = 2
thin = 5
samples = 1000
transient = 500*thin
verbose = 50*thin
# set up spatial variable
n <- nrow(X)
studyDesign = data.frame(sample = as.factor(1:n))
xycoords <- X[c("xcoord", "ycoord")]
rL.spatial = HmscRandomLevel(sData = xycoords)
rL.spatial = setPriors(rL.spatial,nfMin=1,nfMax=10)
fit3 = Hmsc(Y=Y, XData=X, XFormula= ~ cost + fpar + lai + therm + ecoreg,
            studyDesign=studyDesign, ranLevels=list("sample"=rL.spatial), distr="probit")
fit3 = sampleMcmc(fit3, thin = thin, samples = samples, transient = transient,
                  nChains = nChains, verbose = verbose, updater=list(GammaEta=FALSE))
saveRDS(fit3, file = "hmsc_50Kgrid_4x1env_ecoreg_293_spatial_fit.rds")
mpost <- convertToCodaObject(fit3)
effectiveSize(mpost$Beta)
gelman.diag(mpost$Beta)$psrf
plot(mpost$Alpha[[1]])
summary(mpost$Alpha[[1]])$quantiles
VP = computeVariancePartitioning(fit3)
plotVariancePartitioning(fit3, VP=VP, las=2, horiz=F)

# 2-fold cross-validation of PCA model, the hard way
# Hmsc can do this more simply but takes too much time and RAM; process would be simply:
# partition = createPartition(fit1, nfolds = 2)
# preds2 = computePredictedValues(fit1, partition=partition)

X <- read.csv('50Kgrid_pca_511_sites.csv', header=T)
y <- read.csv('Y_all.csv', header=T)
X <- X[X$ID %in% y$ID, , drop = FALSE]
y <- y[y$ID %in% X$ID, , drop = FALSE] # leaves 299 sites present in both X and Y 
dummy <- sample(rep(0:1,each=150)) # set up a random partition
dummy <- dummy[-1] # 149 + 150 = 299
X1 <- X[dummy == 0, ] # split into 2 folds
X2 <- X[dummy == 1, ]
y1 <- y[dummy == 0, ]
y2 <- y[dummy == 1, ]
rownames(y1) <- X1$ID
rownames(y2) <- X2$ID
y1a <- y1[, -1] # remove ID column
y2a <- y2[, -1]
Y1 <- as.matrix(y1a) # convert to matrix Y
Y2 <- as.matrix(y2a)
library(Hmsc)
nChains = 2
thin = 5
samples = 1000
transient = 500*thin
verbose = 50*thin
# set up spatial variable
studyDesign <- data.frame(sample = as.factor(rownames(X1))) # keeps track of row IDs
xycoords <- X1[c("xcoords", "ycoords")]
rL.spatial = HmscRandomLevel(sData = xycoords) # prepare spatial term
rL.spatial = setPriors(rL.spatial,nfMin=1,nfMax=10) # set min and max no. of latent variables
# set formula and fit model
fit1 = Hmsc(Y=Y1, XData=X1, XFormula= ~ env + env.1 + env.2 + env.3 + env.4 + cost,
            studyDesign=studyDesign, ranLevels=list("sample"=rL.spatial), distr="probit")
fit1 = sampleMcmc(fit1, thin = thin, samples = samples, transient = transient,
                  nChains = nChains, verbose = verbose, updater=list(GammaEta=FALSE))
saveRDS(fit1, file = "2cvpcaC299_fit1.rds")
# must re-build spatial element of study design
X3 <- rbind(X1, X2)
studyDesign = data.frame(sample = as.factor(rownames(X3)))
xycoords <- X3[c("xcoords", "ycoords")]
rL.spatial = HmscRandomLevel(sData = xycoords)
rL.spatial = setPriors(rL.spatial,nfMin=1,nfMax=10)
# now predict for all sites, using a model fitted with half the sites
pred1 = Hmsc:::predict.Hmsc(fit1, XData = data.frame(X3), type = "response", studyDesign = studyDesign, 
                            ranLevels=list("sample"=rL.spatial), expected = TRUE)
saveRDS(pred1, file = "2cvpcaC299_pred1.rds")
rm(fit1, pred1)
# set up spatial variable for 2nd fold (swap training and validation sets)
studyDesign <- data.frame(sample = as.factor(rownames(X2)))
xycoords <- X2[c("xcoords", "ycoords")]
rL.spatial = HmscRandomLevel(sData = xycoords)
rL.spatial = setPriors(rL.spatial,nfMin=1,nfMax=10)
# set formula and fit model
fit1 = Hmsc(Y=Y2, XData=X2, XFormula= ~ env + env.1 + env.2 + env.3 + env.4 + cost,
            studyDesign=studyDesign, ranLevels=list("sample"=rL.spatial), distr="probit")
fit1 = sampleMcmc(fit1, thin = thin, samples = samples, transient = transient,
                  nChains = nChains, verbose = verbose, updater=list(GammaEta=FALSE))
saveRDS(fit1, file = "2cvpcaC299_fit2.rds")
# re-build spatial element
studyDesign = data.frame(sample = as.factor(rownames(X3)))
xycoords <- X3[c("xcoords", "ycoords")]
rL.spatial = HmscRandomLevel(sData = xycoords)
rL.spatial = setPriors(rL.spatial,nfMin=1,nfMax=10)
pred1 = Hmsc:::predict.Hmsc(fit1, XData = data.frame(X3), type = "response", studyDesign = studyDesign, 
                            ranLevels=list("sample"=rL.spatial), expected = TRUE)
saveRDS(pred1, file = "2cvpcaC299_pred2.rds")

# compare actual and predicted
r <- Reduce("+", pred1) / length(pred1)
# check reduce works
y <- 0
for (i in 1:2000){
  p <- pred1[[i]][1,1]
  y <- y + p
}
y/2000
# this works on a system with more RAM:
# r <- apply(simplify2array(pred1), 1:2, mean)
saveRDS(r, file = "2cvpcaC299_pred2_matrix.rds")
Y3 <- rbind(Y1, Y2)
actsum <- rowSums(Y3)
pred1 <- readRDS("2cvpcaC299_pred1.rds")
r2 <- Reduce("+", pred1) / length(pred1)
saveRDS(r2, file = "2cvpcaC299_pred1_matrix.rds")
# predicted values are in first half of r and first half of r2
r1a <- r2[150:299, ]
r2a <- r[1:149, ]
r3 <- rbind(r1a, r2a)
predsum <- rowSums(r3)
cor(actsum, predsum , method="spearman")

# prepare folds for 5-fold cross-validation

X <- read.csv('50Kgrid_pca_511_sites.csv', header=T)
y <- read.csv('Y_all.csv', header=T)
X <- X[X$ID %in% y$ID, , drop = FALSE] # remove sites in X not in Y and vice versa
y <- y[y$ID %in% X$ID, , drop = FALSE]
y <- y[, -1] # remove ID column
s <- X[, c("xcoords", "ycoords")]
x <- X[, 10:14, 6] # remove superfluous columns
# x <- cbind(intercept = 1, x) # only if using a categorical variable - this avoids a Matlab problem
x$ecoreg <- as.factor(x$ecoreg)
dummy <- sample(rep(c(1,2,3,4,5),each=60)) # set up random partition (300)
dummy <- dummy[-300] # 299 sites - CHECK that this removes an item from dummy=1, if not, repeat until it does
x1 <- x[dummy == 1, ] # split into 5 folds
x2 <- x[dummy == 2, ]
x3 <- x[dummy == 3, ]
x4 <- x[dummy == 4, ]
x5 <- x[dummy == 5, ]
y1 <- y[dummy == 1, ]
y2 <- y[dummy == 2, ]
y3 <- y[dummy == 3, ]
y4 <- y[dummy == 4, ]
y5 <- y[dummy == 5, ]
s1 <- s[dummy == 1, ]
s2 <- s[dummy == 2, ]
s3 <- s[dummy == 3, ]
s4 <- s[dummy == 4, ]
s5 <- s[dummy == 5, ]
# save folds
write.table(y1,"Y_1_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(y2,"Y_2_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(y3,"Y_3_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(y4,"Y_4_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(y5,"Y_5_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(s1,"S_1_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(s2,"S_2_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(s3,"S_3_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(s4,"S_4_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(s5,"S_5_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(x1,"X_1_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(x2,"X_2_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(x3,"X_3_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(x4,"X_4_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")
write.table(x5,"X_5_plant.csv", row.names = FALSE, col.names = FALSE, sep=",")

############################################
# now do 5-fold cross-validation in Matlab #
############################################

# 5-fold cross-validation, Matlab outputs - check Matlab results vs actual

pth <- "C:/Users/etc" # alter this as needed
WD <- file.path(pth, "bakeoff", "pipeline") # alter this as needed
setwd(WD)	
PD3 <- file.path(WD,"HMSC") # where the Matlab inputs and outputs were stored
Y1 <- read.csv(paste0(PD3, "/", "Y_1_plant", ".csv"), header=FALSE)
Y2 <- read.csv(paste0(PD3, "/", "Y_2_plant", ".csv"), header=FALSE)
Y3 <- read.csv(paste0(PD3, "/", "Y_3_plant", ".csv"), header=FALSE)
Y4 <- read.csv(paste0(PD3, "/", "Y_4_plant", ".csv"), header=FALSE)
Y5 <- read.csv(paste0(PD3, "/", "Y_5_plant", ".csv"), header=FALSE)
y_all <- rbind(Y5, Y1, Y2, Y3, Y4)
yspp <- rowSums(y_all) # actual values: total no. of genera per site
# I renamed matlab output files first eg preds_plant_hmsc3_d5_60.csv becomes d_5.csv
d1 <- read.csv(paste0(PD3, "/", "d1", ".csv"), header=FALSE)
d2 <- read.csv(paste0(PD3, "/", "d2", ".csv"), header=FALSE)
d3 <- read.csv(paste0(PD3, "/", "d3", ".csv"), header=FALSE)
d4 <- read.csv(paste0(PD3, "/", "d4", ".csv"), header=FALSE)
d5 <- read.csv(paste0(PD3, "/", "d5", ".csv"), header=FALSE)
d1a <- d1[(239+1):299,] 
d2a <- d2[(240+1):299,] # only 59 rows in X1, Y1, S1 which are validation sets for Matlab loop dTyp=2 (2 of 1-5)
d3a <- d3[(239+1):299,]
d4a <- d4[(239+1):299,]
d5a <- d5[(239+1):299,]
d_all <- rbind(d1a, d2a, d3a, d4a, d5a)
dmatrix <- array(as.matrix(d_all), dim=list(299, 1157, 100)) # sites x genera x predictions
# get rid of intercept column, leaving us with 1156 columns = genera
dmatrix2 <- dmatrix[,2:1157,]
corvec <- vector()
for (i in 1:100) {
  t1 <- dmatrix2[,,i] # for each prediction matrix
  t1a <- rowSums(t1)
  corvec[i] <- cor(yspp, t1a, method="spearman")
}
folds_5_mean <- apply(dmatrix2, c(1,2), mean)
folds_5_totals <- rowSums(folds_5_mean)
write.csv(folds_5_totals, "pcanewpluscost5cv.csv", row.names=FALSE)
mean(corvec)
sd(corvec)

