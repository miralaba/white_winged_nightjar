
library(maptools)
library(raster)
library(rgeos)
library(rgdal)
library(biomod2)
library(dismo)
library(virtualspecies)
#library(MuMIn)



#### importing occurrence data ####
occur <- read.csv("data/occurrence.txt", sep="\t", header = T)
occur <- occur[complete.cases(occur),]
# subsampling to minimize spatial autocorrelation

points<-SpatialPoints(occur[,c("decimalLongitude", "decimalLatitude")])
proj4string(points) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

points_matrix <- gWithinDistance(points, dist = 0.1, byid = TRUE)
points_matrix[lower.tri(points_matrix, diag=TRUE)] <- NA

occur$stay <- colSums(points_matrix, na.rm=TRUE) == 0

occur2 <- occur[occur$stay==T, ]

#write.csv(occur2, "occurrence_filter.csv", row.names = F)

rm(list=ls()[!(ls() %in% c("occur2"))])

#### creating background by convex hull + buffer ####
pts<-SpatialPoints(occur2[,c("decimalLongitude", "decimalLatitude")])
pres.bkg <- gBuffer(pts, width = 1)
bkg <- gBuffer(pres.bkg, width = 1, byid = T)


#### importing climatic variables ####
bio_curr_list <- list.files("../GIS/wc2.1_30s_bio", pattern = ".tif", full.names = T, recursive = T)
bio_curr <- stack(bio_curr_list)
rm("bio_curr_list")
names(bio_curr)<-gsub(pattern = "wc2.1_30s_", replacement = "", names(bio_curr))


#### Crop the climate layers and clip to background boundary
bio_curr_crop <- crop(bio_curr, extent(bkg))
#plot(bio_curr_crop[[1]])
#bio_curr_bkg <- mask(bio_curr_crop, bkg)
#bio_curr_bkg <- stack(bio_curr_bkg)
#plot(bio_curr_bkg[[1]])
#points(pts)

intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}

bio_curr_bkg <- bio_curr_crop
bio_curr_bkg <- stack(mask(bio_curr_bkg, intersect_mask(bio_curr_bkg))) 


#### removing collinearity ####
removeCollinearity(bio_curr_bkg, multicollinearity.cutoff = 0.7, select.variables = T, plot = T, method = "spearman")
sel.var <- c("bio_1", "elev", "bio_12", "bio_15", "bio_2", "bio_4")
#write.csv(sel.var, "data/selvar.csv")

bio_curr_bkg <- bio_curr_bkg[[grep("bio|elev", paste(sel.var), value = T)]]
#writeRaster(bio_curr_bkg, "raster/selvar_model.grd", format="raster", overwrite=TRUE)



rm(list= ls()[!(ls() %in% c("sel.var", "bkg", "intersect_mask", "bio_curr_bkg", "occur2", "pts"))])
#### Building data matrix ####
myBiomodData <- BIOMOD_FormatingData(resp.var = rep.int(1, times = nrow(occur2)),
                                     expl.var = bio_curr_bkg,
                                     resp.xy = occur2[,c("decimalLongitude", "decimalLatitude")],
                                     resp.name = 'Hcandicans',
                                     PA.nb.rep = 3,
                                     PA.nb.absences = as.numeric(nrow(occur2)*10),
                                     PA.strategy = "disk",
                                     PA.dist.min = 100000)



myBiomodData
nrow(myBiomodData@coord)
head(myBiomodData@coord, 15)
plot(bio_curr_bkg[[1]])
points(myBiomodData@coord[c(1:13),], col="red")
points(myBiomodData@coord[c(14:403),])

#default
myBiomodOption <- BIOMOD_ModelingOptions(GLM = list(control = glm.control(maxit = 15000)),
                                         MAXENT = list(path_to_maxent.jar = '/script', memory_allocated = 2048))
#edited
#myBiomodOption <- BIOMOD_ModelingOptions(GLM = list(type = 'simple', 
#                                                    interaction.level = 3, 
#                                                    control = glm.control(maxit = 1000)),
#                                         
#                                         GAM = list(algo = 'GAM_mgcv', 
#                                                    method = 'REML', 
#                                                    control = list(maxit = 1000)),
#                                         
#                                         RF = list(ntree = 10000),
#                                         
#                                         MAXENT.Phillips = list(path_to_maxent.jar = '/home/leonardo/Documents/app',
#                                                                memory_allocated = 2048))
#


############ modeling #############
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    models = c('GLM', 'MAXENT.Phillips'),
                                    models.options = myBiomodOption,
                                    NbRunEval = 10,
                                    DataSplit = 80,
                                    Yweights = NULL,
                                    VarImport = 5,
                                    models.eval.meth = c('TSS', 'ROC'),
                                    SaveObj = T,
                                    rescal.all.models = T,
                                    do.full.models = F,
                                    modeling.id = 'Hcandicans')





#saving
capture.output(myBiomodData, file = "Hcandicans/Hcandicans_BiomodData.txt")
capture.output(myBiomodModelOut, file = "Hcandicans/Hcandicans_BiomodModelOut.txt")
capture.output(get_evaluations(myBiomodModelOut), file = "Hcandicans/Hcandicans_eval_BiomodModelOut.txt")
capture.output(get_variables_importance(myBiomodModelOut), file = "Hcandicans/Hcandicans_var_importance__BiomodModelOut.txt")

#ev<-get_evaluations(myBiomodModelOut, as.data.frame=T)




#myBiomodModelOut <- load("Cmuelleri/Cmuelleri.Cmuelleri.models.out")
#myBiomodModelOut <- get(myBiomodModelOut)
#rm(Cmuelleri.Cmuelleri.models.out)

#ensemble
myBiomodEM_all <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                          chosen.models =  'all',
                                          em.by = 'all',
                                          eval.metric = c('TSS'),
                                          eval.metric.quality.threshold = 0.7,
                                          prob.mean = F,
                                          prob.cv = F,
                                          prob.ci = F,
                                          prob.ci.alpha = 0.05,
                                          prob.median = F,
                                          committee.averaging = F,
                                          prob.mean.weight = T,
                                          prob.mean.weight.decay = 'proportional',
                                          VarImport = 0)




# saving
capture.output(myBiomodEM_all, file = "Hcandicans/Hcandicans_EM_all.txt")
capture.output(get_evaluations(myBiomodEM_all), file = "Hcandicans/Hcandicans_eval_EM_all.txt")


#### Current ####
bio_curr_list <- list.files("../../../GIS/clima/current_2-5min", pattern = ".tif", full.names = T, recursive = T)
bio_curr <- stack(bio_curr_list)
rm("bio_curr_list")
names(bio_curr)<-gsub(pattern = "wc2.1_2.5m_", replacement = "", names(bio_curr))


#### Crop the climate layers and clip to background boundary
bio_curr_proj <- crop(bio_curr, extent(bkg))
#plot(bio_curr_proj[[1]])
#points(pts)

bio_curr_proj <- stack(mask(bio_curr_proj, intersect_mask(bio_curr_proj))) 

bio_curr_proj <- bio_curr_proj[[grep("bio|elev", paste(sel.var), value = T)]]

#projecting current
myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                  new.env = bio_curr_proj,
                                  proj.name = 'current',
                                  xy.new.env = NULL,
                                  selected.models = 'all',
                                  binary.meth = 'TSS',
                                  compress = F,
                                  build.clamping.mask = F,
                                  do.stack = F,
                                  output.format = '.img')



#projecting current ensemble
myBiomodProj_EM_all <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all,
                                                   projection.output = myBiomodProj,
                                                   new.env = NULL,
                                                   xy.new.env = NULL,
                                                   selected.models = 'all',
                                                   proj.name = 'current_ensemble',
                                                   binary.meth = 'TSS',
                                                   filtered.meth = NULL,
                                                   compress = NULL,
                                                   output.format = '.img',
                                                   total.consensus = T)





rm(list= ls()[(ls() %in% c("bio_curr", "bio_curr_proj", "bio_curr_bkg", "myBiomodProj", "myBiomodProj_EM_all"))])
#plot(raster("Hcandicans/proj_current_ensemble/individual_projections/Hcandicans_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img"))
#

#### uncertainty ####
models.list <- list.files(path = "Hcandicans/proj_current/individual_projections", pattern = ".img$", full.names = T, recursive = T)
models.list <- grep(pattern = "_TSSbin", models.list, value = T, invert = T)
individual_models<-stack(models.list)
#names(individual_models)
#plot(individual_models[[1:4]])
rm(models.list)

Mtotal <- mean(individual_models)
Mglm <- mean(individual_models[[grep("GLM", names(individual_models), value = T)]])
Mmaxent <- mean(individual_models[[grep("MAXENT.Phillips", names(individual_models), value = T)]])


SST <- sum((individual_models-Mtotal)^2)

SSMET <- sum((individual_models[[grep("GLM", names(individual_models), value = T)]]-Mglm)^2)+
         sum((individual_models[[grep("MAXENT.Phillips", names(individual_models), value = T)]]-Mmaxent)^2)
between.met <- SSMET/SST


SSglm <-sum((individual_models[[grep("GLM", names(individual_models), value = T)]]-Mglm)^2)
within.glm <- SSglm/SST
SSmaxent <- sum((individual_models[[grep("MAXENT.Phillips", names(individual_models), value = T)]]-Mmaxent)^2)
within.maxent <- SSmaxent/SST

uncert.part<-stack(between.met, within.glm, within.maxent)
names(uncert.part)<-c("between_algo", "within_glm", "within_maxent")
plot(uncert.part)

writeRaster(uncert.part, "Hcandicans/current_uncertainty.grd", format="raster")



rm(list= ls()[(ls() %in% c("individual_models", "Mtotal", "Mglm", "Mgam", "Mrf", "Mmaxent", "SST", "SSMET", "SSglm", "SSmaxent", "between.met", "within.glm", "within.maxent", "uncert.part"))])


#### PAST:: Holocene ~6kbp  ccsm ####
bio_holcc_list <- list.files("../../../GIS/clima/hol_ccsm/ccmidbi_2-5m", full.names = T, recursive = T)
bio_holcc <- stack(bio_holcc_list)
alt <- raster("../../../GIS/clima/current_2-5min/wc2.1_2.5m_elev.tif")
names(alt)<-"elev"

#### Crop the climate layers and clip to background boundary
bio_holcc_proj <- crop(bio_holcc, extent(bkg))
alt_crop <- crop(alt, extent(bkg))

bio_holcc_proj <- stack(bio_holcc_proj, alt_crop)
bio_holcc_proj <- stack(mask(bio_holcc_proj, intersect_mask(bio_holcc_proj))) 

bio_holcc_proj <- bio_holcc_proj[[grep("bio|elev", paste(sel.var), value = T)]]

#### transform temp variables == °C/10 
bio_holcc_proj[["bio_1"]]<-bio_holcc_proj[["bio_1"]]/10
bio_holcc_proj[["bio_2"]]<-bio_holcc_proj[["bio_2"]]/10
bio_holcc_proj[["bio_4"]]<-bio_holcc_proj[["bio_4"]]/10

rm(list= ls()[(ls() %in% c("bio_holcc_list", "bio_holcc", "bio_holcc_crop", "alt", "alt_crop", "alt_crop_res"))])

#projecting PAST                                                                                                                               
myBiomodProj_holcc <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = bio_holcc_proj,
                                        proj.name = 'holocene_ccsm',
                                        xy.new.env = NULL,
                                        selected.models = "all",
                                        binary.meth = 'TSS',
                                        compress = F,
                                        build.clamping.mask = F,
                                        do.stack=F,
                                        output.format = '.img')

#ensemble projecting PAST
myBiomodProj_EM_holcc_all <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all,
                                                        projection.output = myBiomodProj_holcc,
                                                        new.env = NULL,
                                                        xy.new.env = NULL,
                                                        selected.models = 'all',
                                                        proj.name = 'holocene_ccsm_ensemble',
                                                        binary.meth = 'TSS',
                                                        filtered.meth = NULL,
                                                        compress = NULL,
                                                        output.format = '.img',
                                                        total.consensus = TRUE)


rm(list= ls()[(ls() %in% c("bio_holcc_proj", "myBiomodProj_holcc", "myBiomodProj_EM_holcc_all"))])
#plot(raster("Hcandicans/proj_holocene_ccsm_ensemble/individual_projections/Hcandicans_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img"))
#



#### PAST:: Holocene ~6kbp  miroc ####
bio_holmr_list <- list.files("../../../GIS/clima/hol_miroc/mrmidbi_2-5m", full.names = T, recursive = T)
bio_holmr <- stack(bio_holmr_list)
names(bio_holmr)<-gsub(pattern = "mrmidbi", replacement = "bio_", names(bio_holmr))
alt <- raster("../../../GIS/clima/current_2-5min/wc2.1_2.5m_elev.tif")
names(alt)<-"elev"

#### Crop the climate layers and clip to background boundary
bio_holmr_proj <- crop(bio_holmr, extent(bkg))
alt_crop <- crop(alt, extent(bkg))

bio_holmr_proj <- stack(bio_holmr_proj, alt_crop)
bio_holmr_proj <- stack(mask(bio_holmr_proj, intersect_mask(bio_holmr_proj))) 

bio_holmr_proj <- bio_holmr_proj[[grep("bio|elev", paste(sel.var), value = T)]]

#### transform temp variables == °C/10 
bio_holmr_proj[["bio_1"]]<-bio_holmr_proj[["bio_1"]]/10
bio_holmr_proj[["bio_2"]]<-bio_holmr_proj[["bio_2"]]/10
bio_holmr_proj[["bio_4"]]<-bio_holmr_proj[["bio_4"]]/10

rm(list= ls()[(ls() %in% c("bio_holmr_list", "bio_holmr", "bio_holcc_crop", "alt", "alt_crop", "alt_crop_res"))])

#projecting PAST                                                                                                                               
myBiomodProj_holmr <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = bio_holmr_proj,
                                        proj.name = 'holocene_miroc',
                                        xy.new.env = NULL,
                                        selected.models = "all",
                                        binary.meth = 'TSS',
                                        compress = F,
                                        build.clamping.mask = F,
                                        do.stack=F,
                                        output.format = '.img')

#ensemble projecting PAST
myBiomodProj_EM_holmr_all <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all,
                                                        projection.output = myBiomodProj_holmr,
                                                        new.env = NULL,
                                                        xy.new.env = NULL,
                                                        selected.models = 'all',
                                                        proj.name = 'holocene_miroc_ensemble',
                                                        binary.meth = 'TSS',
                                                        filtered.meth = NULL,
                                                        compress = NULL,
                                                        output.format = '.img',
                                                        total.consensus = TRUE)


rm(list= ls()[(ls() %in% c("bio_holmr_proj", "myBiomodProj_holmr", "myBiomodProj_EM_holmr_all"))])
#plot(raster("Hcandicans/proj_holocene_miroc_ensemble/individual_projections/Hcandicans_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img"))
#


#### uncertainty ####
holcc.list <- list.files(path = "Hcandicans/proj_holocene_ccsm/individual_projections", pattern = ".img$", full.names = T, recursive = T)
holcc.list <- grep(pattern = "_TSSbin", holcc.list, value = T, invert = T)
holmr.list <- list.files(path = "Hcandicans/proj_holocene_miroc/individual_projections", pattern = ".img$", full.names = T, recursive = T)
holmr.list <- grep(pattern = "_TSSbin", holmr.list, value = T, invert = T)
models.list <- c(holcc.list, holmr.list)

individual_models<-stack(models.list)
#names(individual_models)
#plot(individual_models[[1:4]])
rm(list= ls()[(ls() %in% c("models.list", "holcc.list", "holmr.list"))])

Mtotal <- mean(individual_models)
Mcc <- mean(individual_models[[grep("ccsm", names(individual_models), value = T)]])
Mmr <- mean(individual_models[[grep("miroc", names(individual_models), value = T)]])
Mccglm <- mean(individual_models[[grep("ccsm.*GLM", names(individual_models), value = T)]])
Mccmaxent <- mean(individual_models[[grep("ccsm.*MAXENT.Phillips", names(individual_models), value = T)]])
Mmrglm <- mean(individual_models[[grep("miroc.*GLM", names(individual_models), value = T)]])
Mmrmaxent <- mean(individual_models[[grep("miroc.*MAXENT.Phillips", names(individual_models), value = T)]])


SST <- sum((individual_models-Mtotal)^2)

SSGCM <- sum((individual_models[[grep("ccsm", names(individual_models), value = T)]]-Mcc)^2)+
         sum((individual_models[[grep("miroc", names(individual_models), value = T)]]-Mmr)^2)
between.gcm <- SSGCM/SST


SSMET <- sum((individual_models[[grep("ccsm.*GLM", names(individual_models), value = T)]]-Mccglm)^2)+
         sum((individual_models[[grep("ccsm.*MAXENT.Phillips", names(individual_models), value = T)]]-Mccmaxent)^2)+
         sum((individual_models[[grep("miroc.*GLM", names(individual_models), value = T)]]-Mmrglm)^2)+
         sum((individual_models[[grep("miroc.*MAXENT.Phillips", names(individual_models), value = T)]]-Mmrmaxent)^2)
between.algo.within.gcm <- SSMET/SST


SSccglm <-sum((individual_models[[grep("ccsm.*GLM", names(individual_models), value = T)]]-Mccglm)^2)
within.ccglm <- SSccglm/SST
SSccmaxent <- sum((individual_models[[grep("ccsm.*MAXENT.Phillips", names(individual_models), value = T)]]-Mccmaxent)^2)
within.ccmaxent <- SSccmaxent/SST
SSmrglm <-sum((individual_models[[grep("miroc.*GLM", names(individual_models), value = T)]]-Mmrglm)^2)
within.mrglm <- SSmrglm/SST
SSmrmaxent <- sum((individual_models[[grep("miroc.*MAXENT.Phillips", names(individual_models), value = T)]]-Mmrmaxent)^2)
within.mrmaxent <- SSmrmaxent/SST

uncert.part<-stack(between.gcm, between.algo.within.gcm, within.ccglm, within.ccmaxent, within.mrglm, within.mrmaxent)
names(uncert.part)<-c("between_gcm", "between_algo_within_gcm", "within_ccglm", "within_ccmaxent", "within_mrglm", "within_mrmaxent")
plot(uncert.part)

writeRaster(uncert.part, "Hcandicans/holocene_uncertainty.grd", format="raster")



rm(list= ls()[!(ls() %in% c("bkg", "myBiomodData", "myBiomodEM_all", "myBiomodModelOut", "myBiomodOption", "occur2", "pres.bkg", "sel.var", "intersect_mask"))])


#### PAST:: LGM ~21kbp ccsm ####
bio_lgmcc_list <- list.files("../../../GIS/clima/lgm_ccsm/cclgmbi_2-5m", full.names = T, recursive = T)
bio_lgmcc <- stack(bio_lgmcc_list)
names(bio_lgmcc)<-gsub(pattern = "cclgmbi", replacement = "bio_", names(bio_lgmcc))
alt <- raster("../../../GIS/clima/current_2-5min/wc2.1_2.5m_elev.tif")
names(alt)<-"elev"

#### Crop the climate layers and clip to background boundary
bio_lgmcc_proj <- crop(bio_lgmcc, extent(bkg))
alt_crop <- crop(alt, extent(bkg))
bio_lgmcc_proj <- stack(bio_lgmcc_proj, alt_crop)

bio_lgmcc_proj <- stack(mask(bio_lgmcc_proj, intersect_mask(bio_lgmcc_proj))) 

bio_lgmcc_proj <- bio_lgmcc_proj[[grep("bio|elev", paste(sel.var), value = T)]]

#### transform temp variables == °C/10 
bio_lgmcc_proj[["bio_1"]]<-bio_lgmcc_proj[["bio_1"]]/10
bio_lgmcc_proj[["bio_2"]]<-bio_lgmcc_proj[["bio_2"]]/10
bio_lgmcc_proj[["bio_4"]]<-bio_lgmcc_proj[["bio_4"]]/10

rm(list= ls()[(ls() %in% c("bio_lgmcc_list", "bio_lgmcc", "alt", "alt_crop"))])

#projecting PAST                                                                                                                               
myBiomodProj_lgmcc <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = bio_lgmcc_proj,
                                        proj.name = 'lgm_ccsm',
                                        xy.new.env = NULL,
                                        selected.models = "all",
                                        binary.meth = 'TSS',
                                        compress = F,
                                        build.clamping.mask = F,
                                        do.stack=F,
                                        output.format = '.img')

#ensemble projecting PAST
myBiomodProj_EM_lgmcc_all <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all,
                                                        projection.output = myBiomodProj_lgmcc,
                                                        new.env = NULL,
                                                        xy.new.env = NULL,
                                                        selected.models = 'all',
                                                        proj.name = 'lgm_ccsm_ensemble',
                                                        binary.meth = 'TSS',
                                                        filtered.meth = NULL,
                                                        compress = NULL,
                                                        output.format = '.img',
                                                        total.consensus = TRUE)


rm(list= ls()[(ls() %in% c("bio_lgmcc_proj", "myBiomodProj_lgmcc", "myBiomodProj_EM_lgmcc_all"))])
#plot(raster("Hcandicans/proj_lgm_ccsm_ensemble/individual_projections/Hcandicans_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img"))
#


#### PAST:: LGM ~21kbp miroc ####
bio_lgmmr_list <- list.files("../../../GIS/clima/lgm_miroc/mrlgmbi_2-5m", full.names = T, recursive = T)
bio_lgmmr <- stack(bio_lgmmr_list)
names(bio_lgmmr)<-gsub(pattern = "mrlgmbi", replacement = "bio_", names(bio_lgmmr))
alt <- raster("../../../GIS/clima/current_2-5min/wc2.1_2.5m_elev.tif")
names(alt)<-"elev"

#### Crop the climate layers and clip to background boundary
bio_lgmmr_proj <- crop(bio_lgmmr, extent(bkg))
alt_crop <- crop(alt, extent(bkg))
bio_lgmmr_proj <- stack(bio_lgmmr_proj, alt_crop)

bio_lgmmr_proj <- stack(mask(bio_lgmmr_proj, intersect_mask(bio_lgmmr_proj))) 

bio_lgmmr_proj <- bio_lgmmr_proj[[grep("bio|elev", paste(sel.var), value = T)]]

#### transform temp variables == °C/10 
bio_lgmmr_proj[["bio_1"]]<-bio_lgmmr_proj[["bio_1"]]/10
bio_lgmmr_proj[["bio_2"]]<-bio_lgmmr_proj[["bio_2"]]/10
bio_lgmmr_proj[["bio_4"]]<-bio_lgmmr_proj[["bio_4"]]/10

rm(list= ls()[(ls() %in% c("bio_lgmmr_list", "bio_lgmmr", "alt", "alt_crop"))])

#projecting PAST                                                                                                                               
myBiomodProj_lgmmr <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = bio_lgmmr_proj,
                                        proj.name = 'lgm_miroc',
                                        xy.new.env = NULL,
                                        selected.models = "all",
                                        binary.meth = 'TSS',
                                        compress = F,
                                        build.clamping.mask = F,
                                        do.stack=F,
                                        output.format = '.img')

#ensemble projecting PAST
myBiomodProj_EM_lgmmr_all <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all,
                                                        projection.output = myBiomodProj_lgmmr,
                                                        new.env = NULL,
                                                        xy.new.env = NULL,
                                                        selected.models = 'all',
                                                        proj.name = 'lgm_miroc_ensemble',
                                                        binary.meth = 'TSS',
                                                        filtered.meth = NULL,
                                                        compress = NULL,
                                                        output.format = '.img',
                                                        total.consensus = TRUE)


rm(list= ls()[(ls() %in% c("bio_lgmmr_proj", "myBiomodProj_lgmmr", "myBiomodProj_EM_lgmmr_all"))])
#plot(raster("Hcandicans/proj_lgm_miroc_ensemble/individual_projections/Hcandicans_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img"))
#


#### uncertainty ####
lgmcc.list <- list.files(path = "Hcandicans/proj_lgm_ccsm/individual_projections", pattern = ".img$", full.names = T, recursive = T)
lgmcc.list <- grep(pattern = "_TSSbin", lgmcc.list, value = T, invert = T)
lgmmr.list <- list.files(path = "Hcandicans/proj_lgm_miroc/individual_projections", pattern = ".img$", full.names = T, recursive = T)
lgmmr.list <- grep(pattern = "_TSSbin", lgmmr.list, value = T, invert = T)
models.list <- c(lgmcc.list, lgmmr.list)

individual_models<-stack(models.list)
#names(individual_models)
#plot(individual_models[[1:4]])
rm(list= ls()[(ls() %in% c("models.list", "lgmcc.list", "lgmmr.list"))])

Mtotal <- mean(individual_models)
Mcc <- mean(individual_models[[grep("ccsm", names(individual_models), value = T)]])
Mmr <- mean(individual_models[[grep("miroc", names(individual_models), value = T)]])
Mccglm <- mean(individual_models[[grep("ccsm.*GLM", names(individual_models), value = T)]])
Mccmaxent <- mean(individual_models[[grep("ccsm.*MAXENT.Phillips", names(individual_models), value = T)]])
Mmrglm <- mean(individual_models[[grep("miroc.*GLM", names(individual_models), value = T)]])
Mmrmaxent <- mean(individual_models[[grep("miroc.*MAXENT.Phillips", names(individual_models), value = T)]])


SST <- sum((individual_models-Mtotal)^2)

SSGCM <- sum((individual_models[[grep("ccsm", names(individual_models), value = T)]]-Mcc)^2)+
         sum((individual_models[[grep("miroc", names(individual_models), value = T)]]-Mmr)^2)
between.gcm <- SSGCM/SST


SSMET <- sum((individual_models[[grep("ccsm.*GLM", names(individual_models), value = T)]]-Mccglm)^2)+
         sum((individual_models[[grep("ccsm.*MAXENT.Phillips", names(individual_models), value = T)]]-Mccmaxent)^2)+
         sum((individual_models[[grep("miroc.*GLM", names(individual_models), value = T)]]-Mmrglm)^2)+
         sum((individual_models[[grep("miroc.*MAXENT.Phillips", names(individual_models), value = T)]]-Mmrmaxent)^2)
between.algo.within.gcm <- SSMET/SST


SSccglm <-sum((individual_models[[grep("ccsm.*GLM", names(individual_models), value = T)]]-Mccglm)^2)
within.ccglm <- SSccglm/SST
SSccmaxent <- sum((individual_models[[grep("ccsm.*MAXENT.Phillips", names(individual_models), value = T)]]-Mccmaxent)^2)
within.ccmaxent <- SSccmaxent/SST
SSmrglm <-sum((individual_models[[grep("miroc.*GLM", names(individual_models), value = T)]]-Mmrglm)^2)
within.mrglm <- SSmrglm/SST
SSmrmaxent <- sum((individual_models[[grep("miroc.*MAXENT.Phillips", names(individual_models), value = T)]]-Mmrmaxent)^2)
within.mrmaxent <- SSmrmaxent/SST

uncert.part<-stack(between.gcm, between.algo.within.gcm, within.ccglm, within.ccmaxent, within.mrglm, within.mrmaxent)
names(uncert.part)<-c("between_gcm", "between_algo_within_gcm", "within_ccglm", "within_ccmaxent", "within_mrglm", "within_mrmaxent")
plot(uncert.part)

writeRaster(uncert.part, "Cmuelleri/lgm_uncertainty.grd", format="raster")



rm(list= ls()[!(ls() %in% c("bkg", "myBiomodData", "myBiomodEM_all", "myBiomodModelOut", "myBiomodOption", "occur2", "pres.bkg", "sel.var", "intersect_mask"))])


#### PAST:: LIG ~120kbp ####
bio_lig_list <- list.files("../../../GIS/clima/lig/lig_30s_bio", pattern = ".bil$", full.names = T, recursive = T)
bio_lig <- stack(bio_lig_list)
names(bio_lig)<-gsub(pattern = "lig_30s_", replacement = "", names(bio_lig))
alt <- raster("../../../GIS/clima/current_2-5min/wc2.1_2.5m_elev.tif")
names(alt)<-"elev"
#### Crop the climate layers and clip to background boundary
bio_lig_proj <- crop(bio_lig, extent(bkg))
alt_crop <- crop(alt, extent(bkg))

bio_lig_proj <- resample(bio_lig_proj, alt_crop, method="bilinear")
bio_lig_proj <- stack(bio_lig_proj, alt_crop)
bio_lig_proj <- stack(mask(bio_lig_proj, intersect_mask(bio_lig_proj))) 

bio_lig_proj <- bio_lig_proj[[grep("bio|elev", paste(sel.var), value = T)]]

#### transform temp variables == °C/10 
bio_lig_proj[["bio_1"]]<-bio_lig_proj[["bio_1"]]/10
bio_lig_proj[["bio_2"]]<-bio_lig_proj[["bio_2"]]/10
bio_lig_proj[["bio_4"]]<-bio_lig_proj[["bio_4"]]/10

rm(list= ls()[(ls() %in% c("bio_lig_list", "bio_lig", "alt", "alt_crop"))])

#projecting PAST                                                                                                                               
myBiomodProj_lig <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                      new.env = bio_lig_proj,
                                      proj.name = 'lig',
                                      xy.new.env = NULL,
                                      selected.models = "all",
                                      binary.meth = 'TSS',
                                      compress = F,
                                      build.clamping.mask = F,
                                      do.stack=F,
                                      output.format = '.img')

#ensemble projecting PAST
myBiomodProj_EM_lig_all <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all,
                                                      projection.output = myBiomodProj_lig,
                                                      new.env = NULL,
                                                      xy.new.env = NULL,
                                                      selected.models = 'all',
                                                      proj.name = 'lig_ensemble',
                                                      binary.meth = 'TSS',
                                                      filtered.meth = NULL,
                                                      compress = NULL,
                                                      output.format = '.img',
                                                      total.consensus = TRUE)


rm(list= ls()[(ls() %in% c("bio_lig_proj", "myBiomodProj_lig", "myBiomodProj_EM_lig_all"))])
#plot(raster("Hcandicans/proj_lig_ensemble/individual_projections/Hcandicans_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img"))
#


#### uncertainty ####
models.list <- list.files(path = "Hcandicans/proj_lig/individual_projections", pattern = ".img$", full.names = T, recursive = T)
models.list <- grep(pattern = "_TSSbin", models.list, value = T, invert = T)
individual_models<-stack(models.list)
#names(individual_models)
#plot(individual_models[[1:4]])
rm(models.list)

Mtotal <- mean(individual_models)
Mglm <- mean(individual_models[[grep("GLM", names(individual_models), value = T)]])
Mmaxent <- mean(individual_models[[grep("MAXENT.Phillips", names(individual_models), value = T)]])


SST <- sum((individual_models-Mtotal)^2)

SSMET <- sum((individual_models[[grep("GLM", names(individual_models), value = T)]]-Mglm)^2)+
         sum((individual_models[[grep("MAXENT.Phillips", names(individual_models), value = T)]]-Mmaxent)^2)
between.met <- SSMET/SST


SSglm <-sum((individual_models[[grep("GLM", names(individual_models), value = T)]]-Mglm)^2)
within.glm <- SSglm/SST
SSmaxent <- sum((individual_models[[grep("MAXENT.Phillips", names(individual_models), value = T)]]-Mmaxent)^2)
within.maxent <- SSmaxent/SST

uncert.part<-stack(between.met, within.glm, within.maxent)
names(uncert.part)<-c("between_algo", "within_glm", "within_maxent")
plot(uncert.part)

writeRaster(uncert.part, "Cmuelleri/lig_uncertainty.grd", format="raster")



rm(list= ls()[(ls() %in% c("individual_models", "Mtotal", "Mglm", "Mgam", "Mrf", "Mmaxent", "SST", "SSMET", "SSglm", "SSmaxent", "between.met", "within.glm", "within.maxent", "uncert.part"))])


#### FUTURE:: 2050, SSP5-8.5, CNRM-ESM2-1  ####
bio_fut1 <- stack("../../../GIS/clima/2050_ssp585_cnrm-esm2/wc2.1_2.5m_bioc_CNRM-ESM2-1_ssp585_2041-2060.tif")
names(bio_fut1)<-gsub(pattern = "wc2.1_2.5m_bioc_CNRM.ESM2.1_ssp585_2041.2060.", replacement = "bio_", names(bio_fut1))
alt <- raster("../../../GIS/clima/current_2-5min/wc2.1_2.5m_elev.tif")
names(alt)<-"elev"
#### Crop the climate layers and clip to background boundary
bio_fut1_crop <- crop(bio_fut1, extent(bkg))
alt_crop <- crop(alt, extent(bkg))
bio_fut1_proj <- stack(bio_fut1_crop, alt_crop)

bio_fut1_proj <- stack(mask(bio_fut1_proj, intersect_mask(bio_fut1_proj))) 

bio_fut1_proj <- bio_fut1_proj[[grep("bio|elev", paste(sel.var), value = T)]]

rm(list= ls()[(ls() %in% c("bio_fut1_list", "bio_fut1", "bio_fut1_crop", "alt", "alt_crop"))])

#projecting FUTURE
myBiomodProj_F1 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                     new.env = bio_fut1_proj,
                                     proj.name = 'future_2050_ssp585_cnrm-esm2',
                                     xy.new.env = NULL,
                                     selected.models = "all",
                                     binary.meth = 'TSS',
                                     compress = F,
                                     build.clamping.mask = F,
                                     do.stack=F,
                                     output.format = '.img')

#ensemble projecting FUTURE                                                                                                                             
myBiomodProj_EM_F1_all <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all,
                                                     projection.output = myBiomodProj_F1,
                                                     new.env = NULL,
                                                     xy.new.env = NULL,
                                                     selected.models = 'all',
                                                     proj.name = 'future_2050_ssp585_cnrm-esm2_ensemble',
                                                     binary.meth = 'TSS',
                                                     filtered.meth = NULL,
                                                     compress = NULL,
                                                     output.format = '.img',
                                                     total.consensus = TRUE)


rm(list= ls()[(ls() %in% c("bio_fut1_proj", "myBiomodProj_F1", "myBiomodProj_EM_F1_all"))])
#plot(raster("Hcandicans/proj_future_2050_ssp585_cnrm-esm2_ensemble/individual_projections/Hcandicans_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img"))
#


#### FUTURE:: 2050, SSP5-8.5, MIROC6  ####
bio_fut2 <- stack("../../../GIS/clima/2050_ssp585_miroc6/wc2.1_2.5m_bioc_MIROC6_ssp585_2041-2060.tif")
names(bio_fut2)<-gsub(pattern = "wc2.1_2.5m_bioc_MIROC6_ssp585_2041.2060.", replacement = "bio_", names(bio_fut2))
alt <- raster("../../../GIS/clima/current_2-5min/wc2.1_2.5m_elev.tif")
names(alt)<-"elev"
#### Crop the climate layers and clip to background boundary
bio_fut2_crop <- crop(bio_fut2, extent(bkg))
alt_crop <- crop(alt, extent(bkg))
bio_fut2_proj <- stack(bio_fut2_crop, alt_crop)

bio_fut2_proj <- stack(mask(bio_fut2_proj, intersect_mask(bio_fut2_proj))) 

bio_fut2_proj <- bio_fut2_proj[[grep("bio|elev", paste(sel.var), value = T)]]

rm(list= ls()[(ls() %in% c("bio_fut2_list", "bio_fut2", "bio_fut2_crop", "alt", "alt_crop"))])

#projecting FUTURE
myBiomodProj_F2 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                     new.env = bio_fut2_proj,
                                     proj.name = 'future_2050_ssp585_miroc6',
                                     xy.new.env = NULL,
                                     selected.models = "all",
                                     binary.meth = 'TSS',
                                     compress = F,
                                     build.clamping.mask = F,
                                     do.stack=F,
                                     output.format = '.img')

#ensemble projecting FUTURE                                                                                                                             
myBiomodProj_EM_F2_all <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all,
                                                     projection.output = myBiomodProj_F2,
                                                     new.env = NULL,
                                                     xy.new.env = NULL,
                                                     selected.models = 'all',
                                                     proj.name = 'future_2050_ssp585_miroc6_ensemble',
                                                     binary.meth = 'TSS',
                                                     filtered.meth = NULL,
                                                     compress = NULL,
                                                     output.format = '.img',
                                                     total.consensus = TRUE)


rm(list= ls()[(ls() %in% c("bio_fut2_proj", "myBiomodProj_F2", "myBiomodEM_F2_all", "myBiomodProj_EM_F2_all"))])
#plot(raster("Hcandicans/proj_future_2050_ssp585_miroc6_ensemble/individual_projections/Hcandicans_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img"))
#


#### uncertainty ####
fut1.list <- list.files(path = "Hcandicans/proj_future_2050_ssp585_cnrm-esm2/individual_projections", pattern = ".img$", full.names = T, recursive = T)
fut1.list <- grep(pattern = "_TSSbin", fut1.list, value = T, invert = T)
fut2.list <- list.files(path = "Hcandicans/proj_future_2050_ssp585_miroc6/individual_projections", pattern = ".img$", full.names = T, recursive = T)
fut2.list <- grep(pattern = "_TSSbin", fut2.list, value = T, invert = T)
models.list <- c(fut1.list, fut2.list)

individual_models<-stack(models.list)
#names(individual_models)
#plot(individual_models[[1:4]])
rm(list= ls()[(ls() %in% c("models.list", "fut1.list", "fut2.list"))])

Mtotal <- mean(individual_models)
Mcn <- mean(individual_models[[grep("cnrm.esm2", names(individual_models), value = T)]])
Mmr <- mean(individual_models[[grep("miroc6", names(individual_models), value = T)]])
Mcnglm <- mean(individual_models[[grep("cnrm.esm2.*GLM", names(individual_models), value = T)]])
Mcnmaxent <- mean(individual_models[[grep("cnrm.esm2.*MAXENT.Phillips", names(individual_models), value = T)]])
Mmrglm <- mean(individual_models[[grep("miroc6.*GLM", names(individual_models), value = T)]])
Mmrmaxent <- mean(individual_models[[grep("miroc6.*MAXENT.Phillips", names(individual_models), value = T)]])


SST <- sum((individual_models-Mtotal)^2)

SSGCM <- sum((individual_models[[grep("cnrm.esm2", names(individual_models), value = T)]]-Mcn)^2)+
         sum((individual_models[[grep("miroc6", names(individual_models), value = T)]]-Mmr)^2)
between.gcm <- SSGCM/SST


SSMET <- sum((individual_models[[grep("cnrm.esm2.*GLM", names(individual_models), value = T)]]-Mcnglm)^2)+
         sum((individual_models[[grep("cnrm.esm2.*MAXENT.Phillips", names(individual_models), value = T)]]-Mcnmaxent)^2)+
         sum((individual_models[[grep("miroc6.*GLM", names(individual_models), value = T)]]-Mmrglm)^2)+
         sum((individual_models[[grep("miroc6.*MAXENT.Phillips", names(individual_models), value = T)]]-Mmrmaxent)^2)
between.algo.within.gcm <- SSMET/SST


SScnglm <-sum((individual_models[[grep("cnrm.esm2.*GLM", names(individual_models), value = T)]]-Mcnglm)^2)
within.cnglm <- SScnglm/SST
SScnmaxent <- sum((individual_models[[grep("cnrm.esm2.*MAXENT.Phillips", names(individual_models), value = T)]]-Mcnmaxent)^2)
within.cnmaxent <- SScnmaxent/SST
SSmrglm <-sum((individual_models[[grep("miroc6.*GLM", names(individual_models), value = T)]]-Mmrglm)^2)
within.mrglm <- SSmrglm/SST
SSmrmaxent <- sum((individual_models[[grep("miroc6.*MAXENT.Phillips", names(individual_models), value = T)]]-Mmrmaxent)^2)
within.mrmaxent <- SSmrmaxent/SST

uncert.part<-stack(between.gcm, between.algo.within.gcm, within.cnglm, within.cnmaxent, within.mrglm, within.mrmaxent)
names(uncert.part)<-c("between_gcm", "between_algo_within_gcm", "within_cnglm", "within_cnmaxent", "within_mrglm", "within_mrmaxent")
plot(uncert.part)

writeRaster(uncert.part, "Hcandicans/future_uncertainty.grd", format="raster")



rm(list= ls()[!(ls() %in% c("bkg", "myBiomodData", "myBiomodEM_all", "myBiomodModelOut", "myBiomodOption", "occur2", "pres.bkg", "sel.var", "intersect_mask"))])


#



















