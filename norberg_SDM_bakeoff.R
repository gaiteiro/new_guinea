### code adapted from: ###
# Norberg, A., Abrego, N., Blanchet, F. G., Adler, F. R., Anderson, B. J., Anttila, J., Araújo, M. B.,
# Dallas, T., Dunson, D., Elith, J., Foster, S. D., Fox, R., Franklin, J., Godsoe, W., Guisan, A., 
# O'Hara, B., Hill, N. A., Holt, R. D., Hui, F. K. C., . Ovaskainen, O. (2019). 
# A comprehensive evaluation of predictive performance of 33 species distribution models at species and 
# community levels. Ecological Monographs, 89, e01370. https://doi.org/10.1002/ecm.1370

# we need a directory: /bakeoff/pipeline with subdirectories:
# /DATA
# FITS/plant
# /MODELS
# /PREDICT
# /PREDICTIONS/plant
# /RESULTS_final
# /SCRIPTS

# we need the following files in folder 'DATA':
# Yt_1_plant.csv, Yv_1_plant.csv: training and validation sets - genus presence/absence (1/0) by site
# Xt_1_plant.csv, Xv_1_plant.csv: training and validation sets - 5 independent variables for model input
# St_1_plant.csv, Sv_1_plant.csv: training and validation sets - 2 columns, X and Y coordinate respectively
# Don't put row or column names in any of these - no. of rows should equal no. of sites - must be in correct order
# There must be equal number of sites/rows in training and validation sets

# Preliminaries
#-----------------------------------------------------------------------------------------
rm(list = ls(all = TRUE)); gc()	# clear workspace

pth <- "C:/Users/thisUser/myPath" 	# set root directory

WD <- file.path(pth, "bakeoff", "pipeline") # define working directory to pth/bakeoff/pipeline

setwd(WD)	# set working directory to the pipeline folder
SD <- file.path(WD, "SCRIPTS")  # define path to pipeline scripts

OS <- "win" # set your OS

# all these directories must be created

# data directory
DD <- file.path(WD,"DATA")

# model fitting scripts directory
MD <- file.path(WD,"MODELS")

# model fits directory
FD <- file.path(WD,"FITS")
FD2 <- FD

# prediction scripts directory
PD <- file.path(WD,"PREDICT")

# predictions output directory
PD2 <- file.path(WD,"PREDICTIONS")

# final results directory
RDfinal <- file.path(WD,"RESULTS_final")

# models: names list
#-----------------------------------------------------------------------------------------

# all except 10 are commented out as they couldn't be used with my dataset of genera
mod_names <- list("GAM1","GAM2",
                  #"GLM1","GLM8","GLM12",
                  "GLM10",#"GLM11",
                  #"GLM2","GLM3",
                  #"GLM6",
                  #"GLM9",
                  "MRTS1",
                  "GNN1",
                  #"RF1",
                  "BRT1",
                  "XGB1",
                  "SVM1",
                  #"MARS1","MARS2",
                  #"GJAM1",
                  #"SAM1",
                  #"MISTN1", "MISTN2",
                  #"BORAL1", "BORAL2",#"GLM7","BORAL1",
                  #"BC1","BC2",
                  #"GLM4","GLM5","GLM13",
                  #"SSHMSC1",
                  "SSHMSC2",
                  #"HMSC1",
                  "HMSC2","HMSC3","HMSC4"
                  )

# predictions: names list
#-----------------------------------------------------------------------------------------
# this list will be used in naming output files
pred_names	<-	list("gam1_PAs_","gam_spat1_PAs_",
                   #"glm1_PAs_","glm1b_PAs_","glm1_intXs_PAs_",
                   "glmnet1_PAs_",#"glmnet1b_PAs_",
                   #"glmmPQL1_PAs_","glmmPQLspat1_PAs_",
                   #"manyglm1_PAs_",
                   #"traitglm1_PAs_",
                   "mrt1_PAs_",
                   "gnn1_PAs_",
                   #"rf1_PAs_",
                   "brt1_PAs_",
                   "xgb1_PAs_",
                   "svm1_PAs_",
                   #"mars1_PAs_","mars2_PAs_",
                   #"gjam1_PAs_",
                   #"sam1_PAs_",
                   #"mstnt1_PAs_", "mstnt2_PAs_",
                   #"boral1_PAs_","boral2_PAs_",
                   #"bc1_PAs_","bc2_PAs_",
                   #"ss_hmsc1_PAs_",
                   "ss_hmsc2_PAs_",#"ss_hmsc1_intXs_PAs_",
                   "hmsc1_PAs_",
                   "hmsc2_PAs_",
                   "hmsc3_PAs_"#,
                   #"hmsc1_intXs_PAs_"
                   )

names(pred_names)<-mod_names

if (length(mod_names)!=length(pred_names)) {
  stop("Prediction objects and their names are of different size")
}

#-----------------------------------------------------------------------------------------

require(doParallel)
no_cores <- detectCores() - 1 
registerDoParallel(cl = makeCluster(no_cores))
sz <-1 # was list of 3 sizes - for dataset subsetting - we will use the full dataset only
set_no <- "plant" # label for dataset used - was originally list of 6
Sets<- c("plant") # for use in labelling output files
setwd(DD)
# the variables below were used to switch between SDM alternatives
# we will ignore these options but need to set variables for some subroutines
MCMC2 <- FALSE
commSP <- FALSE
intXs <- FALSE


# select random samples from the whole training set = half of sites
#-----------------------------------------------------------------------------------------

# subsetting species present at least once in the training set

absentSp<-list()
for (i in sz) {
  y_tmp <- read.csv(paste("Yt_",i,"_",set_no,".csv",sep=""),header=FALSE)
  absentSp[[i]] <- which(colSums(y_tmp)==0) # identifies empty columns in training set - i.e. absent species
}
absentSp<-as.numeric(unlist(absentSp))
spSel<-1:ncol(y_tmp)
if (sum(absentSp)!=0) { 
  spSel<-spSel[-absentSp]
}
# spSel is a list of column numbers, for each remaining species. absentSp refers to the other columns

write.table(spSel, file=paste("spSel_",Sets[1],".csv",sep=""),sep=",",row.names=F,col.names=F)
save(spSel, file=paste("spSel_",Sets[1],".RData",sep=""))

samprow <- nrow(y_tmp) # another marker for a subsetting exercise which we will not use
samp = 1:samprow
setwd(DD)

# training sets
#-----------------------------------------------------------------------------------------
y_train <- list() 
x_train <- list()
s_train <- list()

for (i in sz) {
  y_train[[i]] <- read.csv(paste("Yt_",i,"_",set_no,".csv",sep=""),header=FALSE)
  y_train[[i]] <- apply(y_train[[i]], 2, as.numeric)[samp,spSel] # all rows, spSel cols only

  x_train[[i]] <- as.matrix(read.csv(paste("Xt_",i,"_",set_no,".csv",sep=""),header=FALSE))
  
  s_train[[i]] <- read.csv(paste("St_",i,"_",set_no,".csv",sep=""),header=FALSE)
  
  colnames(s_train[[i]])<-paste('Rand',1:ncol(s_train[[i]]),sep='')
  
  ncovar<-ncol(x_train[[i]]) # number of predictors, which for this exercise should be 5
  
  for (k in 1:ncovar) {
    x_train[[i]]<-cbind(x_train[[i]],x_train[[i]][,k]^2) # copy each column and square it
  }
  
  x_train[[i]]<-apply(x_train[[i]],2,scale)
  x_train[[i]]<-cbind(1,x_train[[i]]) # column of 1s, for intercept, call it IC
  colnames(x_train[[i]])<-c('IC',paste('V',1:ncovar,sep=''),paste('V',1:ncovar,'_2',sep=''))
  
}


# validation sets: same procedure as for training
#-----------------------------------------------------------------------------------------
y_valid<-list()
x_valid<-list()
s_valid<-list()

for (i in sz) {
  y_valid[[i]] <- read.csv(paste("Yv","_",i,"_",set_no,".csv",sep=""),header=FALSE)
  y_valid[[i]] <- apply(y_valid[[i]], 2, as.numeric)[,spSel]
  
  
  x_valid[[i]] <- as.matrix(read.csv(paste("Xv","_",i,"_",set_no,".csv",sep=""),header=FALSE))

  s_valid[[i]] <- read.csv(paste("Sv","_",i,"_",set_no,".csv",sep=""),header=FALSE)

  colnames(s_valid[[i]])<-paste('Rand',1:ncol(s_valid[[i]]),sep='')
  
  ncovar<-ncol(x_valid[[i]])
  
  for (k in 1:ncovar) {
    x_valid[[i]]<-cbind(x_valid[[i]],x_valid[[i]][,k]^2)
  }
  
  x_valid[[i]]<-apply(x_valid[[i]],2,scale)
  x_valid[[i]]<-cbind(1,x_valid[[i]])
  colnames(x_valid[[i]])<-c('IC',paste('V',1:ncovar,sep=''),paste('V',1:ncovar,'_2',sep=''))
  
}

# lists: for some SDMs we will need a list of lists with 1 sublist per species
#-----------------------------------------------------------------------------------------
DD_t <- list()
DD_v <- list()

for (i in sz) {
  nsp <- ncol(y_train[[i]]) # no. of remaining species
  
  dd_t <- list()
  for (j in 1:nsp) {
    dd_t[[j]] <- data.frame(cbind(y_train[[i]][,j],x_train[[i]],s_train[[i]]))
    colnames(dd_t[[j]]) <- c('sp',colnames(x_train[[i]]),colnames(s_train[[i]]))
  }	
  
  dd_v <- list()
  for (j in 1:nsp) {
    dd_v[[j]] <- data.frame(cbind(y_valid[[i]][,j],x_valid[[i]],s_valid[[i]]))
    colnames(dd_v[[j]]) <- c('sp',colnames(x_valid[[i]]),colnames(s_valid[[i]]))
  }	
  
  DD_t[[i]]<-dd_t # = produces a list each with n=nsp sublists; one for each species
  # columns as follows: species, intercept, 5 variables, 5 squared variables, X & Y coordinates
  # "sp" "IC" "V1" "V2" "V3" "V4" "V5" "V1_2" "V2_2" "V3_2" "V4_2" "V5_2" "Rand1" "Rand2"
  # each sublist has n=no. of sites rows
  # 1st column (sp) = is the species present at this site (1=yes)
  DD_v[[i]]<-dd_v
  
}

#-----------------------------------------------------------------------------------------

setwd(WD)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIT MODELS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# I have commented out most of these - leaving the ones that worked on my dataset
# some never converged, others required too much RAM e.g. BayesComm

# model fitting will take seconds to days, depending on SDM type
# some of these packages may need installing of course

MCMC2 <- FALSE
commSP <- FALSE
intXs <- FALSE
dataN <- samprow # no. of sites
REPs <- 100		# no. of prediction replicates
dszs <- c(1) # index not used so set to 1
Sets <- "NGplants"# label
betaInds<-c("sim", "nest", "sor") # predictive measures

MD <- file.path(WD,"MODELS") # which is where all the model-fitting R files should be

#Sys.time()
#print('mistnet')
#source(file.path(MD,"fit.mistnet.r"))  
#### mistnet will fail but 'source' loads environment which we can correct....
#### ...the following trace call opens the mistnet_fit_optimx subroutine 
#trace(mistnet_fit_optimx, edit=TRUE)
#### we need to change change one word in the mistnet_fit_optimx file: xtimes to xtime
#### then try again
#source(file.path(MD,"fit.mistnet.r"))  	# mistnet

#Sys.time()
#print('boral')
#source(file.path(MD,"fit.boral.r"))  	# BORAL

Sys.time()
print('xgb')
source(file.path(MD,"fit.xgb.r"))		# XGB

#Sys.time()
#print('mars')
#source(file.path(MD,"fit.mars.r"))		# MARS

# Sys.time()
# print('traitglm')
# source(file.path(MD,"fit.traitglm.r"))	# MVABUND
# not much point using this without trait data

# Sys.time()
#print('brt')
#source(file.path(MD,"fit.brt.r"))		# BRT

#Sys.time()
#print('rf')
#source(file.path(MD,"fit.rf.r"))		# RFs

#Sys.time()
#intXs <- TRUE
#source(file.path(MD,"fit.glm.r"))  # basic glm

#Sys.time()
#intXs <- FALSE
#print('gjam')
#source(file.path(MD,"fit.gjam.r"))  	# GJAMS

#Sys.time()
#print('BayesComm')
#source(file.path(MD,"fit.bc.r"))		# BayesComm
# very slow and huge memory requirements

#Sys.time()
#print('glm')
#source(file.path(MD,"fit.glm.r"))		# single-species GLM 1

#Sys.time()
#print('glmmPQL')
#source(file.path(MD,"fit.glmmPQL.r"))	# ssGLM 2

Sys.time()
print('mrt')
source(file.path(MD,"fit.mrt.r"))		# MRTs

Sys.time()
print('svm')
source(file.path(MD,"fit.svm.r"))		# SVM

Sys.time()
print('gnn')
source(file.path(MD,"fit.gnn.r"))		# GNNs

#Sys.time()
#print('manyglm')
#source(file.path(MD,"fit.manyglm.r"))	# MVABUND

Sys.time()
print('glmnet')
source(file.path(MD,"fit.glmnet.r"))		# GLMNET
 
#Sys.time()
#print('sam')
#source(file.path(MD,"fit.sam.r"))		# SAMs

Sys.time()
print('gam')
source(file.path(MD,"fit.gam.r"))		# GAM
Sys.time()

# must run the 2 Matlab files separately - in Matlab, these are in MODELs folder:
# fit_predict_hmsc.m - this will output 3 models; basic, grid ID as iid, spatial term as iid
# need to check variable for no. or rows in data e.g. dSizes=[146,292];
# single-species ensemble - this will output 2 models; basic, spatial term as iid
# fit_predict_hmsc_singleSDM.m


stopCluster(cl) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GET PREDICTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FD2 <- FD
PD <- file.path(WD,"PREDICT") # which is where all these R files should be

# ss GLM 1
#source(file.path(PD, "predict.glm.r"))

# ss GLM 2
#source(file.path(PD, "predict.glmmPQL.r"))

# BRT
#source(file.path(PD, "predict.brt.r"))

# GNNs
source(file.path(PD, "predict.gnn.r"))

# MARS
#source(file.path(PD, "predict.mars.r"))

# MVABUND
#source(file.path(PD, "predict.manyglm.r", sep = ""))

# traitGLM
#source(file.path(PD, "predict.traitglm.r"))

# GAM
source(file.path(PD, "predict.gam.r"))

# RFs
#source(file.path(PD, "predict.rf.r"))

# MRTs
source(file.path(PD, "predict.mrt.r"))

# XGB
source(file.path(PD, "predict.xgb.r"))

# GJAMS
#source(file.path(PD, "predict.gjam.r"))

# SAMs
#source(file.path(PD, "predict.sam.r"))

# mistnet
#source(file.path(PD, "predict.mistnet.r"))

# BayesComm
#source(file.path(PD, "predict.bc.r"))

# BORAL
#source(file.path(PD, "predict.boral.r"))

# SVM
source(file.path(PD, "predict.svm.r"))

# GLMnet
source(file.path(PD, "predict.glmnet.r"))

# ssHMSC
source(paste(PD, "predict.hmsc.all.r", sep = ""))

# HMSC

source(paste(PD, "predict.hmsc.ss.r", sep = ""))

#-----------------------------------------------------------------------------------------

# calculate performance measures
#-----------------------------------------------------------------------------------------

setwd(WD)

MCMC2 <- FALSE

source(file.path(SD, "list_preds.r"))

ens <- list( 'ALL' = 'all') # we only use one ensemble measure to pool results of tested SDMs
predtypes <- c("interpol", "extrapol", "extrapol2") # we only use the first of these - simple interpolation
NsitePairs <- dataN[1] # rows in training set = rows in validation set
Nsamples <- 100

pms_res <- NA

pred_list <- vector("list", length(mod_names))
names(pred_list) <- mod_names # from list of models used
dsz <- dataN[1] # this also refers to rows in training set = rows in validation set
Nsamples <- 100

j <- 1
for (i in 1:length(pred_list)) {
  pred_list[[i]] <- vector("list", 1)
  tmp <- paste0(pred_names[[i]], 
                paste(j, dsz, sep = "_"))
  names(pred_list[[i]]) <- tmp
}
preds <- list()
directory = file.path

# these are various subroutines which should be in pipeline/SCRIPTS

source(file.path(SD, "list_preds.r"))
preds$predictions <- list_preds(directory = file.path(PD2, set_no), 
                                pred_list = pred_list)

source(file.path(SD, "modify_predictions.r"))
mod_preds <- modify_predictions(predictions = preds, 
                                yvalid = y_valid[[j]],
                                n_of_site_pairs = NsitePairs)

source(file.path(SD, "create_ensembles.r"))
ens_preds <- create_ensembles(ensembles = ens,
                              predictions = preds,
                              nsamples = Nsamples)

source(file.path(SD, "modify_predictions.r"))
mod_ens_preds <- modify_predictions(predictions = ens_preds, 
                                    yvalid = y_valid[[j]],
                                    n_of_site_pairs = NsitePairs)

source(file.path(SD, "calculate_pms.r"))
tmp <- cbind(predtypes[j],
             dsz,
             set_no,
             calculate_pms(modified_predictions = mod_preds))  

source(file.path(SD, "calculate_pms_ens.r"))
tmp_ens <- cbind(predtypes[j],
                 dsz,
                 set_no,
                 calculate_pms(modified_predictions = mod_ens_preds)) 

pms_res <- rbind(pms_res,
                 tmp,
                 tmp_ens)

pms_res_o <- pms_res[-1,]
pms_res <- pms_res_o
head(pms_res)
dim(pms_res)
pms_res[, c(1, 3:5)] <- apply(pms_res[, c(1, 3:5)], 2, as.character)
pms_res$modelling.framework[which(pms_res_o$modelling.framework=="ENS" & pms_res_o$model.variant=="ALL")] <- "ENS_all"

filebody <- "pms_res"
save(pms_res, file = file.path(RDfinal, paste0(filebody,".RData")))


# combine performance measures result tables and model features
#-----------------------------------------------------------------------------------------

setwd(WD)
commSP <- FALSE

filebody <- "pms_res"
RDfinal <- file.path(WD,"RESULTS_final")
load(file = file.path(RDfinal, "pms_res.RData"))
pms_res_allsp <- pms_res
rm(pms_res)

pms_res_allsp <- cbind(pms_res_allsp, 0)
colnames(pms_res_allsp)[ncol(pms_res_allsp)] <- "prev_threshold"
pms_res_comb <- pms_res_allsp 

write.table(pms_res_comb, 
            file = file.path(RDfinal, "pms_res_tbl.csv"), 
            sep = ",", 
            row.names = FALSE,
            col.names = TRUE)

feats <- read.csv(file.path(DD, "feats.csv"), header = TRUE)
# this file should come with the Norberg data. It categorizes all the SDM types' model features
# NOTE you will have to delete all rows in 'feats' not referring to an SDM type you have tested successfully
# i.e. PAs .RData files must be available with predictions 

head(pms_res_comb)
sum(is.na(pms_res_comb))
sum(is.nan(unlist(pms_res_comb)))

pms_res_tbl <- NA
for (m in 1:nrow(feats)) {
  #message(feats[m,])
  tmp <- cbind(feats[m,], 
               pms_res_comb[which(pms_res_comb$modelling.framework==feats$Abbreviation[m]),])
  pms_res_tbl <- rbind(pms_res_tbl, tmp)
}
pms_res_tbl <- pms_res_tbl[-1,]
head(pms_res_tbl)

pms_res_tbl <- pms_res_tbl[,-c(14:15)]
colnames(pms_res_tbl)[which(colnames(pms_res_tbl)=="predtypes[j]")] <- "predType"
colnames(pms_res_tbl)[which(colnames(pms_res_tbl)=="dsz")] <- "dataSize"
colnames(pms_res_tbl)[which(colnames(pms_res_tbl)=="set_no")] <- "dataSet"

pms_res_tbl$predType[which(pms_res_tbl$predType=="interpol")] <- "i"
pms_res_tbl$predType[which(pms_res_tbl$predType=="extrapol")] <- "e"
pms_res_tbl$predType[which(pms_res_tbl$predType=="extrapol2")] <- "e2"
pms_res_tbl <- pms_res_tbl[,c(1:10, 12, 11, 13, 14:ncol(pms_res_tbl))]

pms_res_tbl_names <- colnames(pms_res_tbl)
pms_res_tbl_names[which(pms_res_tbl_names=="accuracy_expected_me")] <- "accuracy1"
pms_res_tbl_names[which(pms_res_tbl_names=="discrimination_auc")] <- "discrimination1"
pms_res_tbl_names[which(pms_res_tbl_names=="calibration_prob_bin_mse")] <- "calibration1"
pms_res_tbl_names[which(pms_res_tbl_names=="precision_sqrt_probs")] <- "precision1"
pms_res_tbl_names[which(pms_res_tbl_names=="sp_rich_rmse")] <- "accuracy2site"
pms_res_tbl_names[which(pms_res_tbl_names=="beta_sim_rmse")] <- "accuracy3beta1"
pms_res_tbl_names[which(pms_res_tbl_names=="beta_sne_rmse")] <- "accuracy3beta2"
pms_res_tbl_names[which(pms_res_tbl_names=="beta_sor_rmse")] <- "accuracy3beta3"
pms_res_tbl_names[which(pms_res_tbl_names=="sp_rich_spear")] <- "discrimination2site"
pms_res_tbl_names[which(pms_res_tbl_names=="beta_sim_spear")] <- "discrimination3beta1"
pms_res_tbl_names[which(pms_res_tbl_names=="beta_sne_spear")] <- "discrimination3beta2"
pms_res_tbl_names[which(pms_res_tbl_names=="beta_sor_spear")] <- "discrimination3beta3"
pms_res_tbl_names[which(pms_res_tbl_names=="sp_rich_predint")] <- "calibration2site"
pms_res_tbl_names[which(pms_res_tbl_names=="beta_sim_predint")] <- "calibration3beta1"
pms_res_tbl_names[which(pms_res_tbl_names=="beta_sne_predint")] <- "calibration3beta2"
pms_res_tbl_names[which(pms_res_tbl_names=="beta_sor_predint")] <- "calibration3beta3"
pms_res_tbl_names[which(pms_res_tbl_names=="sp_rich_sd")] <- "precision2site"
pms_res_tbl_names[which(pms_res_tbl_names=="beta_sim_sd")] <- "precision3beta1"
pms_res_tbl_names[which(pms_res_tbl_names=="beta_sne_sd")] <- "precision3beta2"
pms_res_tbl_names[which(pms_res_tbl_names=="beta_sor_sd")] <- "precision3beta3"
cbind(pms_res_tbl_names, colnames(pms_res_tbl))
colnames(pms_res_tbl) <- pms_res_tbl_names

save(pms_res_tbl, file = file.path(RDfinal, "pms_res_tbl.RData"))
write.table(pms_res_tbl, 
            file = file.path(RDfinal, "pms_res_tbl.csv"), 
            sep = ",", 
            row.names = FALSE,
            col.names = TRUE)


#-----------------------------------------------------------------------------------------

# plot performance measures
#-----------------------------------------------------------------------------------------

setwd(WD)

commSP <- FALSE
mod_names2 <- mod_names

filebody <- "pms_res"
if (commSP) {
  filebody <- paste0(filebody, "_commSP")
}
load(file = file.path(RDfinal, paste0(filebody, ".RData")))
resTBL <- data.frame(pms_res)
head(resTBL)

pms <- colnames(resTBL[,which(colnames(resTBL)=="accuracy_expected_me"):(ncol(resTBL))])
pms_ls <- as.list(pms)
names(pms_ls)[which(pms_ls=="accuracy_expected_me")] <- "PM 1 a"
names(pms_ls)[which(pms_ls=="discrimination_auc")] <- "PM 1 b"
names(pms_ls)[which(pms_ls=="calibration_prob_bin_mse")] <- "PM 1 c"
names(pms_ls)[which(pms_ls=="precision_sqrt_probs")] <- "PM 1 d"
names(pms_ls)[which(pms_ls=="sp_rich_rmse")] <- "PM 2 a"
names(pms_ls)[which(pms_ls=="beta_sim_rmse")] <- "PM 3sim a"
names(pms_ls)[which(pms_ls=="beta_sne_rmse")] <- "PM 3nest a"
names(pms_ls)[which(pms_ls=="beta_sor_rmse")] <- "PM 3sor a"
names(pms_ls)[which(pms_ls=="sp_rich_spear")] <- "PM 2 b"
names(pms_ls)[which(pms_ls=="beta_sim_spear")] <- "PM 3sim b"
names(pms_ls)[which(pms_ls=="beta_sne_spear")] <- "PM 3nest b"
names(pms_ls)[which(pms_ls=="beta_sor_spear")] <- "PM 3sor b"
names(pms_ls)[which(pms_ls=="sp_rich_predint")] <- "PM 2 c"
names(pms_ls)[which(pms_ls=="beta_sim_predint")] <- "PM 3sim c"
names(pms_ls)[which(pms_ls=="beta_sne_predint")] <- "PM 3nest c"
names(pms_ls)[which(pms_ls=="beta_sor_predint")] <- "PM 3sor c"
names(pms_ls)[which(pms_ls=="sp_rich_sd")] <- "PM 2 d"
names(pms_ls)[which(pms_ls=="beta_sim_sd")] <- "PM 3sim d"
names(pms_ls)[which(pms_ls=="beta_sne_sd")] <- "PM 3nest d"
names(pms_ls)[which(pms_ls=="beta_sor_sd")] <- "PM 3sor d"

filebody <- paste0(RDfinal, "/raw_res_fig") 
if (commSP) {
  filebody <- paste0(filebody, '_commSP')
}

### plot raw results
pdf(file = paste0(filebody,".pdf"), 
    bg = "transparent", 
    width = 15, 
    height = 15)
par(family = "serif", 
    mfrow = c(5,4), 
    mar = c(7,3,2,1))
for (p in 1:length(pms)) {
  plot(0,0,
       xlim = c(0,length(mod_names2)),
       ylim = c(min(resTBL[,pms[p]], na.rm=T), max(resTBL[,pms[p]], na.rm = TRUE)),
       type = "n",
       xaxt = "n",
       xlab = "",
       ylab = "",
       main = pms[p])
  for (m in 1:length(mod_names)) {
    lines(resTBL[which(resTBL$modelling.framework==mod_names[[m]]), pms[p]],
          x = rep(m, 
                  times = length(resTBL[which(resTBL$modelling.framework==mod_names[[m]]),pms[p]])),
          lwd = 2)
    points(mean(resTBL[which(resTBL$modelling.framework==mod_names[[m]]),pms[p]]),
           x = m, 
           pch = 21, 
           col = "black", 
           bg = "red3",
           cex = 2)
  }
  axis(1, 1:length(mod_names), unlist(mod_names), las=2)
}
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# End of pipeline
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



