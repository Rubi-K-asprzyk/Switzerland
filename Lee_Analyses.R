

# Working directory
setwd("~/Documents/PhD-Thesis/Research/Switzerland")

#------------------------#
##### INITIALISATION #####
#------------------------#

  #### . Packages . ##### 

# Install/load pacman. 
suppressPackageStartupMessages(if(!require(pacman)){install.packages("pacman");library(pacman)})

# Install/load tons of packages.
p_load("vegan",
       "ecospat",
       "rdacca.hp",
       "ade4",
       "factoextra",
       "adehabitatHR",
       "raster",
       "rgdal",
       "maptools",
       "rgeos",
       "ecospat",
       "VennDiagram",
       "betapart",
       "spdep",
       "prioritizr",
       "RRF",
       "glmnet",
       "pls",
       "covsel",
       "tidyverse",
       "caret"
)


  #### . Data Preparation . ##### 

# -- Bryophytes -- # 

# Load the datasets
  # Occurence
bryoData <- read.csv("BryophyteData.csv",row.names=1) 
  # Traits
bryoTraits <- read.csv("Bryophytes_Traits_SpRichness.csv",row.names=1)
  # Environment Variables
bryoEnv <- read.csv("EnvDataForBryophytes.csv",row.names=1) # mnt_mean was changed directly in the dataset to "z"
alti <- bryoEnv$z # Modified to follow Flavien's script afterwards. 

# Load the phylotrees
Liver_Tree <- read.tree("timetree50mod-liverwortsV2.nwk")
Mosses_Tree <- read.tree("timetree50mod-mossesV2.nwk")

# -- Computation of the wanted Metrics -- #

  # Species richness
SR.bryo <- apply(bryoData[,5:ncol(bryoData)],1,sum)
  # PD

  # MPD

# -- Split between Mosses and Liverworts -- #

# Load the names of the mosses present in the dataframe of Occ_Data_Moss for further splitting between Mosses and Liverworts.
Bryo_Names <- read.csv(file = "Utilities/Bryophyte_Names.csv", row.names = 1)
Moss_Names <- dplyr::select(filter(Bryo_Names,Taxa == "Moss"),Species) %>% as.matrix()
Liver_Names <- dplyr::select(filter(Bryo_Names,Taxa == "Liverwort"),Species) %>% as.matrix()

# Extract the meta data names
Bryo_Meta <- colnames(bryoData)[1:4]

# Split the Occ_Data_Moss between Mosses and Liverworts
liverData <- dplyr::select(bryoData,any_of(c(Bryo_Meta,Liver_Names)))
mossData <- dplyr::select(bryoData,any_of(c(Bryo_Meta,Moss_Names)))

# Get the names of the liverworts and mosses present
Moss_Names <- colnames(mossData)[-c(1:4)]
Liver_Names <- colnames(liverData)[-c(1:4)]

# -- Tracheophytes -- # 

# Load the datasets
  # Occurence
TracheoOcc <- read.csv("TracheophyteData.csv",row.names=1)
  # Traits
tracheoTraits <- read.csv("Traits_Tracheophytes_Fraction.csv",row.names = 1)
  # Environment Variables
TracheoEnv <- read.csv("EnvDataForTracheophytes.csv",row.names = 1)

# Load the phylotrees
# OPTION 2: Zanne 2014 / Three keys to the radiation of angiosperms into freezing environments
ZanneTree <- read.tree(file = "PhyloTree/Zanne2014/Vascular_Plants_rooted.dated.tree")

# -- Computation of the wanted Metrics -- #

  # Species richness
SR.tracheo <- apply(TracheoData[,3:ncol(TracheoData)],1,sum)
  # PD

  # MPD

# --- Modify the data --- #

# Split between MetaData and OccurenceData
TracheoMeta <- TracheoOcc[,1:2]       # Verification: table(as.matrix(TracheoOcc))
TracheoOcc <- TracheoOcc[,-c(1:2)]

# Get the species names
Sp_names <- colnames(TracheoOcc)

# Replace "." by "_"
Sp_names <- gsub(".","_",Sp_names,fixed = T)
# Keep only the genuses names
Gn_names <- sub("_.*","",Sp_names) %>%
  unique() # Keep the genuses names


#-------------------------#
##### Alpha Diversity #####
#-------------------------#

cor(SR.bryo,SR.tracheo) # 0.320

test <- glm(SR.bryo~alti+I(alti^2),family = poisson())
with(summary(test), 1 - deviance/null.deviance) # Explained deviance = 0.143

summary(test) 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -4.690e-01  1.589e-01  -2.952  0.00315 ** 
#   alti         3.282e-03  1.959e-04  16.754  < 2e-16 ***
#   I(alti^2)   -9.041e-07  5.822e-08 -15.531  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 2744.2  on 412  degrees of freedom
# Residual deviance: 2353.0  on 410  degrees of freedom
# AIC: 3823.4
# 
# Number of Fisher Scoring iterations: 5



alti.pred <- as.data.frame(c(min(alti):max(alti)))
colnames(alti.pred)="alti"
SR.bryo.pred <- predict(test,newdata = alti.pred,type="response",se.fit=T)
min.SR.bryo.pred <- SR.bryo.pred$fit - SR.bryo.pred$se.fit
max.SR.bryo.pred <- SR.bryo.pred$fit + SR.bryo.pred$se.fit
SR.bryo.pred <- SR.bryo.pred$fit

test <- glm(SR.tracheo~alti+I(alti^2),family = poisson())
with(summary(test), 1 - deviance/null.deviance) # Explained deviance = 0.255
summary(test)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  2.893e+00  6.270e-02   46.14   <2e-16 ***
#   alti         1.649e-03  8.751e-05   18.84   <2e-16 ***
#   I(alti^2)   -6.759e-07  2.845e-08  -23.76   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 5754.6  on 412  degrees of freedom
# Residual deviance: 4289.2  on 410  degrees of freedom
# AIC: 6393.6
alti.pred <- as.data.frame(c(min(alti):max(alti)))
colnames(alti.pred)="alti"
SR.tracheo.pred <- predict(test,newdata = alti.pred,type="response",se.fit=T)
min.SR.tracheo.pred <- SR.tracheo.pred$fit - SR.tracheo.pred$se.fit
max.SR.tracheo.pred <- SR.tracheo.pred$fit + SR.tracheo.pred$se.fit
SR.tracheo.pred <- SR.tracheo.pred$fit

plot(SR.tracheo.pred~alti.pred[,1])
alti.pred[which.max(SR.tracheo.pred),1] # 1219.84m
alti.pred[which.max(SR.bryo.pred),1] # 1814.844m

fin <- rbind.data.frame(cbind.data.frame(SR=SR.tracheo.pred,Taxa="Tracheophytes",alti.pred,mini = min.SR.tracheo.pred, maxi = max.SR.tracheo.pred),
                        cbind.data.frame(SR=SR.bryo.pred,Taxa="Bryophytes",alti.pred,mini = min.SR.bryo.pred, maxi = max.SR.bryo.pred))

# png("CorrelationsBetweenSRTracheoAndBryo/Fig2.png",res=500,width=6000,height = 4500)

p <- ggplot(fin,aes(x=alti,y=SR,fill=Taxa,color=Taxa,ymin = mini, ymax = maxi))+geom_line(linewidth=2)+ 
  geom_segment(aes(x=alti.pred[which.max(SR.tracheo.pred),1],xend=alti.pred[which.max(SR.tracheo.pred),1],y=0,yend=max(SR.tracheo.pred)) , linetype = "dotted",color="#56876D",linewidth=2)+ 
  geom_segment(aes(x=alti.pred[which.max(SR.bryo.pred),1],xend=alti.pred[which.max(SR.bryo.pred),1],y=0,yend=max(SR.bryo.pred)), linetype = "dotted",color="#B68F68",linewidth=2)+
  geom_ribbon(alpha= 0.2,linetype="dotted",linewidth=0.5)+
  scale_color_discrete(type=c("#B68F68","#56876D"))+
  scale_fill_discrete(type=c("#B68F68","#56876D"))+
  ylab("Predicted Species richness")+xlab("Elevation (m)")+
  theme_classic()+
  theme(axis.text.x = element_text(size=15,angle = 30, hjust = 1),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title =element_text(size=15),
        legend.position = "bottom")
p
dev.off()


cor(SR.bryo,alti) #0.200
cor(SR.tracheo,alti) #-0.409

#-------------------------------#
##### Alpha Phylo Diversity #####
#-------------------------------#


cor(SR.bryo,SR.tracheo) # 0.320

test <- glm(SR.bryo~alti+I(alti^2),family = poisson())
with(summary(test), 1 - deviance/null.deviance) # Explained deviance = 0.143

summary(test) 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -4.690e-01  1.589e-01  -2.952  0.00315 ** 
#   alti         3.282e-03  1.959e-04  16.754  < 2e-16 ***
#   I(alti^2)   -9.041e-07  5.822e-08 -15.531  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 2744.2  on 412  degrees of freedom
# Residual deviance: 2353.0  on 410  degrees of freedom
# AIC: 3823.4
# 
# Number of Fisher Scoring iterations: 5



alti.pred <- as.data.frame(c(min(alti):max(alti)))
colnames(alti.pred)="alti"
SR.bryo.pred <- predict(test,newdata = alti.pred,type="response",se.fit=T)
min.SR.bryo.pred <- SR.bryo.pred$fit - SR.bryo.pred$se.fit
max.SR.bryo.pred <- SR.bryo.pred$fit + SR.bryo.pred$se.fit
SR.bryo.pred <- SR.bryo.pred$fit

test <- glm(SR.tracheo~alti+I(alti^2),family = poisson())
with(summary(test), 1 - deviance/null.deviance) # Explained deviance = 0.255
summary(test)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  2.893e+00  6.270e-02   46.14   <2e-16 ***
#   alti         1.649e-03  8.751e-05   18.84   <2e-16 ***
#   I(alti^2)   -6.759e-07  2.845e-08  -23.76   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 5754.6  on 412  degrees of freedom
# Residual deviance: 4289.2  on 410  degrees of freedom
# AIC: 6393.6
alti.pred <- as.data.frame(c(min(alti):max(alti)))
colnames(alti.pred)="alti"
SR.tracheo.pred <- predict(test,newdata = alti.pred,type="response",se.fit=T)
min.SR.tracheo.pred <- SR.tracheo.pred$fit - SR.tracheo.pred$se.fit
max.SR.tracheo.pred <- SR.tracheo.pred$fit + SR.tracheo.pred$se.fit
SR.tracheo.pred <- SR.tracheo.pred$fit

plot(SR.tracheo.pred~alti.pred[,1])
alti.pred[which.max(SR.tracheo.pred),1] # 1219.84m
alti.pred[which.max(SR.bryo.pred),1] # 1814.844m

fin <- rbind.data.frame(cbind.data.frame(SR=SR.tracheo.pred,Taxa="Tracheophytes",alti.pred,mini = min.SR.tracheo.pred, maxi = max.SR.tracheo.pred),
                        cbind.data.frame(SR=SR.bryo.pred,Taxa="Bryophytes",alti.pred,mini = min.SR.bryo.pred, maxi = max.SR.bryo.pred))

# png("CorrelationsBetweenSRTracheoAndBryo/Fig2.png",res=500,width=6000,height = 4500)

p <- ggplot(fin,aes(x=alti,y=SR,fill=Taxa,color=Taxa,ymin = mini, ymax = maxi))+geom_line(linewidth=2)+ 
  geom_segment(aes(x=alti.pred[which.max(SR.tracheo.pred),1],xend=alti.pred[which.max(SR.tracheo.pred),1],y=0,yend=max(SR.tracheo.pred)) , linetype = "dotted",color="#56876D",linewidth=2)+ 
  geom_segment(aes(x=alti.pred[which.max(SR.bryo.pred),1],xend=alti.pred[which.max(SR.bryo.pred),1],y=0,yend=max(SR.bryo.pred)), linetype = "dotted",color="#B68F68",linewidth=2)+
  geom_ribbon(alpha= 0.2,linetype="dotted",linewidth=0.5)+
  scale_color_discrete(type=c("#B68F68","#56876D"))+
  scale_fill_discrete(type=c("#B68F68","#56876D"))+
  ylab("Predicted Species richness")+xlab("Elevation (m)")+
  theme_classic()+
  theme(axis.text.x = element_text(size=15,angle = 30, hjust = 1),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title =element_text(size=15),
        legend.position = "bottom")
p
dev.off()


cor(SR.bryo,alti) #0.200
cor(SR.tracheo,alti) #-0.409

##################################

### Lee ----

nbg <- graph2nb(gabrielneigh(bryoData[,c("x","y")]), sym = TRUE) #Spatial links
listw <- nb2listw(nbg)
Li <- lee(SR.bryo,SR.tracheo,listw,n=413) #Global = 0.1275
lee.mc(SR.bryo,SR.tracheo,listw,nsim=1000, alternative="greater",spChk=T)
# Monte-Carlo simulation of Lee's L
# 
# data:  SR.bryo ,  SR.tracheo 
# weights: listw  
# number of simulations + 1: 1001 
# 
# statistic = 0.12757, observed rank = 909, p-value = 0.09191
# alternative hypothesis: greater


test <- glm(Li$localL~alti+I(alti^2))
with(summary(test), 1 - deviance/null.deviance) # Explained deviance = 0.269

summary(test) 

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.042e-01  1.513e-01   2.010    0.045 *  
#   alti        -8.171e-04  1.993e-04  -4.100 4.99e-05 ***
#   I(alti^2)    3.794e-07  6.125e-08   6.194 1.43e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 0.1885117)
# 
# Null deviance: 105.67  on 412  degrees of freedom
# Residual deviance:  77.29  on 410  degrees of freedom
# AIC: 487.9
# 
# Number of Fisher Scoring iterations: 2

alti.pred <- as.data.frame(c(min(alti):max(alti)))
colnames(alti.pred)="alti"
SR.bryo.pred <- predict(test,newdata = alti.pred,type="response")

plot(SR.bryo.pred~alti.pred[,1],type="l",xlab="Elevation (m)", ylab="Lee's L")
points(Li$localL~alti)

#### Li acrosss the altitudinal bands ----
sum(alti<1000) #64
sum(alti<1400 & alti>=1000)#67
sum(alti<1800 & alti>=1400)#72
sum(alti<2000 & alti>=1800)#62
sum(alti<2200 & alti>=2000)#77
sum(alti>=2200) #71

altibandMin <- c(0,1000,1400,1800,2000,2200)
altibandMax <- c(1000,1400,1800,2000,2200,Inf)
fin <- NULL



tracheo.Traits<- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/Traits_Tracheophytes_Fraction.csv",row.names = 1)


bryo.Traits <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/Bryophytes_Traits_SpRichness.csv",row.names = 1)
alti <- bryoEnv$mnt_mean

source("CorrelationsBetweenSRTracheoAndBryo/CovselModifGaussian.R")
##Test full

covdata <- cbind(bryoEnv,bryo.Traits,tracheo.Traits)
resp <- Li$localL

bb <- covsel.filterGau(covdata = covdata,pa=resp)
dd <- covsel.embedGau(covdata = bb,pa=resp,ntree=1000)

Rf <-(dd$ModelRF)
Rf$rsq[length(Rf$rsq)]
# Mean of squared residuals: 0.1488341
# % Var explained: 41.83

mod <- dd$ModelGLM
coef(mod,s=mod$lambda.1se)
mod$glmnet.fit$dev.ratio[mod$glmnet.fit$lambda == mod$lambda.1se] #0.34


## Selected Variables
write.table(dd$rf.beta,"CorrelationsBetweenSRTracheoAndBryo/covsel/FullDataRF_ImportanceVar.txt",sep="\t")
write.table(dd$glm.beta,"CorrelationsBetweenSRTracheoAndBryo/covsel/FullDataGLM_ImportanceVar.txt",sep="\t")
saveRDS(mod,"CorrelationsBetweenSRTracheoAndBryo/covsel/LM.rds")
saveRDS(Rf,"CorrelationsBetweenSRTracheoAndBryo/covsel/RRF.rds")

fin <- fin2 <- fin3<- NULL
for(i in 1:length(altibandMin)){
  ########
  print(cor(SR.bryo[alti>=altibandMin[i] & alti<altibandMax[i]],SR.tracheo[alti>=altibandMin[i] & alti<altibandMax[i]]))
  
  
  ## Preparation of the data
  petMat.bryo <- bryoData[alti>=altibandMin[i] & alti<altibandMax[i],5:ncol(bryoData)]
  petMat.bryo <- petMat.bryo[,apply(petMat.bryo,2,sum)>0]
  
  SR.bryo.alti <- apply(petMat.bryo,1,sum)
  
  petMat.bryo.env <- bryoEnv[alti>=altibandMin[i] & alti<altibandMax[i],]
  petMat.bryo.env <- petMat.bryo.env[,apply(petMat.bryo.env,2,sum)!=0]
  
  env <- dudi.pca(petMat.bryo.env,scannf = F,nf=ncol(petMat.bryo.env)-1,scale=T)
  NAxesToKeep <- which(cumsum(100*env$eig / sum(env$eig )) >90)[1]
  env <- env$li[,1:NAxesToKeep]
  colnames(env) = paste0("Axis_Env",1:ncol(env))
  
  petMat.bryo.traits <- bryo.Traits[alti>=altibandMin[i] & alti<altibandMax[i],]
  petMat.bryo.traits <- petMat.bryo.traits[,apply(petMat.bryo.traits,2,sum)>0]
  
  bryo.Traits.pca <- dudi.pca(petMat.bryo.traits,scannf = F,nf=ncol(petMat.bryo.traits)-1,scale=T)
  NAxesToKeep <- which(cumsum(100*bryo.Traits.pca$eig / sum(bryo.Traits.pca$eig )) >90)[1]
  bryo.Traits.pca <- bryo.Traits.pca$li[,1:NAxesToKeep]
  colnames(bryo.Traits.pca) = paste0("Axis_BryoTraits",1:ncol(bryo.Traits.pca))
  
  petMat.tracheo <- TracheoData[alti>=altibandMin[i] & alti<altibandMax[i],3:739]
  petMat.tracheo <- petMat.tracheo[,apply(petMat.tracheo,2,sum)>0]
  
  
  SR.tracheo.alti <- apply(petMat.tracheo,1,sum)
  
  
  petMat.tracheo.traits <- tracheo.Traits[alti>=altibandMin[i] & alti<altibandMax[i],]
  petMat.tracheo.traits <- petMat.tracheo.traits[,apply(petMat.tracheo.traits,2,sum)>0]
  
  tracheo.Traits.pca <- dudi.pca(petMat.tracheo.traits,scannf = F,nf=ncol(petMat.tracheo.traits)-1,scale=T)
  NAxesToKeep <- which(cumsum(100*tracheo.Traits.pca$eig / sum(tracheo.Traits.pca$eig )) >90)[1]
  tracheo.Traits.pca <- tracheo.Traits.pca$li[,1:NAxesToKeep]
  colnames(tracheo.Traits.pca) = paste0("Axis_TracheoTraits",1:ncol(tracheo.Traits.pca))
  
  print("Computation of LEE")
  nbg <- graph2nb(gabrielneigh(bryoData[alti>=altibandMin[i] & alti<altibandMax[i],c("x","y")]), sym = TRUE) #Spatial links
  listw <- nb2listw(nbg)
  Li.alti <- lee(SR.bryo.alti,SR.tracheo.alti,listw,n=length(SR.bryo.alti)) #Global = 0.1275
  Li.Global <- Li.alti$L
  Li.Local <- Li.alti$localL
  StatLi <- lee.mc(SR.bryo.alti,SR.tracheo.alti,listw,nsim=1000, alternative="greater",spChk=F)
  pval <- StatLi$p.value
  
  ToStore <- cbind.data.frame(elevationalBand = paste0("[",altibandMin[i],";", altibandMax[i],"["), Li = Li.Global, pval = pval)
  fin <- rbind.data.frame(fin, ToStore)
  
  covdata <- cbind(petMat.bryo.env,petMat.tracheo.traits,petMat.bryo.traits)
  ## Remove variables that do not vary among plots
  covdata <- covdata[,apply(covdata,2,function(x){
    if(length(unique(x))>2){
      return(T)
    }else{
      return(F)
    }
  })]
  covdata <- cbind(resp=Li.Local,covdata)
  
  pls.model <- plsr(formula=resp~.,data=covdata, scale=TRUE, validation="CV")
  cv = RMSEP(pls.model)
  best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1 #1
  summary(pls.model)
  write.table(summary(pls.model),paste0("CorrelationsBetweenSRTracheoAndBryo/plsr/PLS_",altibandMin[i],"_",altibandMax[i],".txt"),sep="\t") #did not run need to check if it runs
  # Data: 	X dimension: 71 98 
  # Y dimension: 71 1
  # Fit method: kernelpls
  # Number of components considered: 62
  # 
  # VALIDATION: RMSEP
  # Cross-validated using 10 random segments.
  # (Intercept)  1 comps  2 comps  3 comps  4 comps  5 comps  6 comps  7 comps  8 comps  9 comps  10 comps  11 comps  12 comps  13 comps  14 comps  15 comps  16 comps  17 comps
  # CV           0.872   0.7652    1.256    1.981    2.035    2.167    3.150    3.359    3.420    3.371     3.057     2.898     2.978     2.860     2.604     2.597     2.777     2.943
  # adjCV        0.872   0.7615    1.209    1.889    1.939    2.062    2.992    3.190    3.247    3.200     2.902     2.751     2.827     2.714     2.471     2.465     2.636     2.794
  # 18 comps  19 comps  20 comps  21 comps  22 comps  23 comps  24 comps  25 comps  26 comps  27 comps  28 comps  29 comps  30 comps  31 comps  32 comps  33 comps  34 comps  35 comps
  # CV        3.016     3.170     2.964     2.968     3.036     3.206     3.432     3.705     4.255     4.431     4.723     4.952     4.997     5.113     5.367     5.729     5.900     6.053
  # adjCV     2.863     3.008     2.813     2.817     2.881     3.043     3.257     3.516     4.039     4.206     4.483     4.701     4.744     4.853     5.094     5.438     5.601     5.746
  # 36 comps  37 comps  38 comps  39 comps  40 comps  41 comps  42 comps  43 comps  44 comps  45 comps  46 comps  47 comps  48 comps  49 comps  50 comps  51 comps  52 comps  53 comps
  # CV        6.147     6.243     6.341     6.450     6.535     6.680     6.795     7.017     7.013     7.033     7.107     7.143     7.052     6.989     6.963     6.861     6.854     6.818
  # adjCV     5.836     5.926     6.020     6.123     6.204     6.341     6.451     6.661     6.657     6.676     6.747     6.781     6.694     6.635     6.610     6.513     6.506     6.472
  # 54 comps  55 comps  56 comps  57 comps  58 comps  59 comps  60 comps  61 comps  62 comps
  # CV        6.866     6.867     6.938     7.023     7.024     7.027     7.024     7.023     7.023
  # adjCV     6.518     6.519     6.586     6.667     6.667     6.671     6.668     6.666     6.666
  # 
  # TRAINING: % variance explained
  # 1 comps  2 comps  3 comps  4 comps  5 comps  6 comps  7 comps  8 comps  9 comps  10 comps  11 comps  12 comps  13 comps  14 comps  15 comps  16 comps  17 comps  18 comps  19 comps
  # X       27.12    33.99    39.46    46.99    51.49    54.61    57.00    59.80    63.33     65.37     67.01     68.39     70.53     72.35     74.47     76.03     77.80     79.14     80.11
  # resp    37.25    63.41    75.97    81.16    85.95    90.13    93.11    94.41    95.22     96.00     96.62     97.18     97.52     97.81     98.03     98.36     98.59     98.73     98.88
  # 20 comps  21 comps  22 comps  23 comps  24 comps  25 comps  26 comps  27 comps  28 comps  29 comps  30 comps  31 comps  32 comps  33 comps  34 comps  35 comps  36 comps  37 comps
  # X        81.33     82.33     83.31     84.13     85.37     86.20     86.90     87.76     88.49     89.20      89.9     90.42     91.29     91.89     92.41     93.06     93.49     94.02
  # resp     98.96     99.04     99.12     99.23     99.29     99.35     99.39     99.41     99.45     99.48      99.5     99.56     99.61     99.67     99.71     99.73     99.77     99.80
  # 38 comps  39 comps  40 comps  41 comps  42 comps  43 comps  44 comps  45 comps  46 comps  47 comps  48 comps  49 comps  50 comps  51 comps  52 comps  53 comps  54 comps  55 comps
  # X        94.72     95.26     95.70     96.29     96.62     97.02     97.57     97.69     97.84     98.47     98.61     98.96     99.29     99.43     99.46     99.54     99.64     99.71
  # resp     99.82     99.83     99.85     99.86     99.88     99.91     99.92     99.95     99.96     99.96     99.97     99.97     99.98     99.98     99.99     99.99     99.99     99.99
  # 56 comps  57 comps  58 comps  59 comps  60 comps  61 comps  62 comps
  # X        99.76     99.84     99.87     99.89     99.91     99.94     99.96
  # resp    100.00    100.00    100.00    100.00    100.00    100.00    100.00
  
  ## Here only one component is useful
  coefficients = coef(pls.model)
  sum.coef = sum(sapply(coefficients, abs))
  coefficients = coefficients * 100 / sum.coef
  coefficients = sort(coefficients[, 1 , 1],decreasing = T)
  
  write.table(paste0("CorrelationsBetweenSRTracheoAndBryo/plsr/Coefficient_",altibandMin[i],"_",altibandMax[i],".txt"),sep="\t")
  
  # barplot(head(coefficients, 5))
  # barplot(tail(coefficients, 5))
  # moliniid                bio1_t_8110_LV95                  eastSlope_mean ch_edaphic_eivdescombes_pixel_h                     SCD40_final
  # 3.51269113                      2.77763480                      2.69900188                      2.30908750                      2.30883041
  # ch_edaphic_eivdescombes_pixel_d            bio9_tdryq_8110_LV95            bio14_pdry_8110_LV95                    loiseleureid                 northSlope_mean
  # 2.02527092                      2.01705175                      1.94809865                      1.88406745                      1.86968207
  # Amblystegiids                           sradY             CanopyVaud2m_median                        TRI_mean              bio15_ps_8110_LV95
  # 1.81980262                      1.70152015                      1.66531612                      1.66080699                      1.22735255
  # bio6_tminc_8110_LV95             bio3_tiso_8110_LV95                      asperulids    meanGCI_annualmean_2017_2021            bio13_pwet_8110_LV95
  # 1.15553860                      1.03224834                      1.02975976                      0.96706284                      0.95057264
  # empetrid                      fragariids                      myrtillids                    Polytrichids                    AI_8110_LV95
  # 0.92698564                      0.92345528                      0.89032282                      0.88115062                      0.86878262
  # ch_edaphic_eivdescombes_pixel_n     meanEVI_annualmin_2017_2021                          Bryids                       primulids                    Anomodontids
  # 0.83877603                      0.81409673                      0.73248002                      0.71106210                      0.69526327
  # arrhenaterid   meanEVI_annualrange_2017_2021                      epipactids             Rocks_Min80Per_prop                       aparinids
  # 0.64037759                      0.60964473                      0.48190819                      0.44090514                      0.41725592
  # bio16_pwetq_8110_LV95          bio18_pwarmq_8110_LV95                        trisetid                         thymids                     ranunculids
  # 0.41245048                      0.41245048                      0.40748674                      0.34578812                      0.31925321
  # bio7_tar_8110_LV95              bio2_tdr_8110_LV95                        linariid            bio5_tmaxw_8110_LV95 ch_edaphic_eivdescombes_pixel_k
  # 0.28911955                      0.28817135                      0.26866605                      0.12295830                      0.07419457
  # CanopyVaud2m_mean                     illecebrids                      slope_mean           bio17_pdryq_8110_LV95          bio11_tcoldq_8110_LV95
  # 0.05881744                     -0.03392920                     -0.06978115                     -0.15695729                     -0.16801941
  # arabids                      Bazzaniids               bio12_p_8110_LV95 ch_edaphic_eivdescombes_pixel_l                    hierochloids
  # -0.18456318                     -0.19438459                     -0.20970842                     -0.23552668                     -0.23983809
  # gdd0Y_8110_LV95            bio8_twetq_8110_LV95          bio10_twarmq_8110_LV95     meanGCI_annualmin_2017_2021                        mnt_mean
  # -0.33377295                     -0.34896414                     -0.34896414                     -0.43411248                     -0.43864462
  # ajugids                      verbascids     meanEVI_annualmax_2017_2021                          bellid                CanopyVaud2m_max
  # -0.46715483                     -0.49946727                     -0.51381593                     -0.52277175                     -0.56423322
  # gdd3Y_8110_LV95     meanEVI_annualstd_2017_2021 ch_edaphic_eivdescombes_pixel_f                     aegopodiids ch_edaphic_eivdescombes_pixel_r
  # -0.56638499                     -0.56775153                     -0.57935831                     -0.59586622                     -0.63512077
  # north_mean                       Thuidiids                 gdd5Y_8110_LV95     meanGCI_annualstd_2017_2021                        TPI_mean
  # -0.64793194                     -0.68610889                     -0.71138177                     -0.72566107                     -0.75691957
  # CanopyVaud2m_min                 CanopyVaud2m_sd                    gypsophilids                       east_mean   meanGCI_annualrange_2017_2021
  # -0.80622840                     -0.86116351                     -0.87226053                     -0.88802749                     -0.89232036
  # roughness_mean     meanGCI_annualmax_2017_2021    meanEVI_annualmean_2017_2021                     Pleuroziids                  etpY_8110_LV95
  # -0.89900615                     -0.93885936                     -1.03999362                     -1.20293046                     -1.21965067
  # Plagiotheciids          bio19_pcoldq_8110_LV95                          nardid                    Glacier_prop                         oxalids
  # -1.31720739                     -1.32000612                     -1.34296412                     -1.54244665                     -1.60469984
  # meanEVI_annualmedian_2017_2021                     Peltigerids  meanGCI_annualmedian_2017_2021                    androsaceids                     SCD15_final
  # -1.66481770                     -1.84514303                     -1.99921907                     -2.12181866                     -2.34266972
  # RocheMeubler_prop ch_edaphic_eivdescombes_pixel_w               bio4_ts_8110_LV95
  # -2.37461212                     -3.03972665                     -4.99391382
  
  covdata <- cbind(env,bryo.Traits.pca,tracheo.Traits.pca)
  resp <- Li.Local
  print("Variable Selection")
  bb <- covsel.filterGau(covdata = covdata,pa=resp)
  dd <- covsel.embedGau(covdata = bb,pa=resp,ntree=500)
  
  Rf <-(dd$ModelRF)
  rsqRF <- Rf$rsq[length(Rf$rsq)]
  
  
  mod <- dd$ModelGLM
  
  rsqlm <- mod$glmnet.fit$dev.ratio[mod$glmnet.fit$lambda == mod$lambda.1se] #0.34
  ToKeep <- cbind.data.frame(elevationalBand = paste0("[",altibandMin[i],";", altibandMax[i],"["), rsqRF,rsqlm)
  fin3 <- rbind.data.frame(fin3, ToKeep)
  
  ## Selected Variables
  # write.table(dd$rf.beta,paste0("CorrelationsBetweenSRTracheoAndBryo/covsel/FullDataRF_",altibandMin[i],"_",altibandMax[i],"_ImportanceVar.txt"),sep="\t")
  # write.table(dd$glm.beta,paste0("CorrelationsBetweenSRTracheoAndBryo/covsel/FullDataGLM_",altibandMin[i],"_",altibandMax[i],"_ImportanceVar.txt"),sep="\t")
  # saveRDS(mod,paste0("CorrelationsBetweenSRTracheoAndBryo/covsel/LM_",altibandMin[i],"_",altibandMax[i],".rds"))
  # saveRDS(Rf,paste0("CorrelationsBetweenSRTracheoAndBryo/covsel/RRF_",altibandMin[i],"_",altibandMax[i],".rds"))
  
  write.table(dd$rf.beta,paste0("CorrelationsBetweenSRTracheoAndBryo/covsel/FullDataRF_",altibandMin[i],"_",altibandMax[i],"_ImportanceVar_PCA.txt"),sep="\t")
  write.table(dd$glm.beta,paste0("CorrelationsBetweenSRTracheoAndBryo/covsel/FullDataGLM_",altibandMin[i],"_",altibandMax[i],"_ImportanceVar_PCA.txt"),sep="\t")
  saveRDS(mod,paste0("CorrelationsBetweenSRTracheoAndBryo/covsel/LM_",altibandMin[i],"_",altibandMax[i],"_PCA.rds"))
  saveRDS(Rf,paste0("CorrelationsBetweenSRTracheoAndBryo/covsel/RRF_",altibandMin[i],"_",altibandMax[i],"_PCA.rds"))

  
}

write.csv2(fin,"CorrelationsBetweenSRTracheoAndBryo/LiGlobal_altiBand.csv")
write.csv2(fin3,"CorrelationsBetweenSRTracheoAndBryo/covsel/AltiCovselModPerf.csv")

barplot(fin[,2],names.arg=fin[,1])


### Turn-Over -----
ToKeep <- bryoData[,3:324]
ToKeep <- apply(ToKeep,1,sum) > 0
sum(ToKeep) #375
BryoTO <- beta.pair(bryoData[ToKeep,3:324])$beta.sim
bryoEnv_dist <- dist(scale(bryoEnv[ToKeep,]))
tracheoTraits_dist <- dist(tracheoTraits[ToKeep,])
b <- dist(alti[ToKeep])
mantel(BryoTO,b)
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = BryoTO, ydis = b) 
# 
# Mantel statistic r: 0.3876 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0129 0.0166 0.0209 0.0242 
# Permutation: free
# Number of permutations: 999

mantel.partial(BryoTO,bryoEnv_dist,tracheoTraits_dist)
# Partial Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel.partial(xdis = BryoTO, ydis = bryoEnv_dist, zdis = tracheoTraits_dist) 
# 
# Mantel statistic r: 0.2774 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0174 0.0223 0.0257 0.0312 
# Permutation: free
# Number of permutations: 999
mantel.partial(BryoTO, tracheoTraits_dist, bryoEnv_dist) 
# Partial Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel.partial(xdis = BryoTO, ydis = tracheoTraits_dist, zdis = bryoEnv_dist) 
# 
# Mantel statistic r: 0.08862 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0175 0.0220 0.0251 0.0289 
# Permutation: free
# Number of permutations: 999



ToKeep <- TracheoData[,3:739]
ToKeep <- apply(ToKeep,1,sum) > 0
sum(ToKeep) #400
TracheoTO <- beta.pair(TracheoData[ToKeep,3:739])$beta.sim
BryoTraits_dist <- dist(bryoTraits[ToKeep,])
TracheoEnv_dist <- dist(scale(TracheoEnv[ToKeep,]))
b <- dist(alti[ToKeep])
# mantel(TracheoTO,b)
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = TracheoTO, ydis = b) 
# 
# Mantel statistic r: 0.5311 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0184 0.0237 0.0267 0.0351 
# Permutation: free
# Number of permutations: 999
mantel.partial(TracheoTO, TracheoEnv_dist, BryoTraits_dist) 
# Partial Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel.partial(xdis = TracheoTO, ydis = TracheoEnv_dist, zdis = BryoTraits_dist) 
# 
# Mantel statistic r: 0.4504 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0247 0.0336 0.0399 0.0445 
# Permutation: free
# Number of permutations: 999
mantel.partial(TracheoTO, BryoTraits_dist, TracheoEnv_dist) 
# Partial Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel.partial(xdis = TracheoTO, ydis = BryoTraits_dist, zdis = TracheoEnv_dist) 
# 
# Mantel statistic r: 0.1661 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0239 0.0321 0.0398 0.0503 
# Permutation: free
# Number of permutations: 999

fin2 <- NULL
fin1<- NULL
for(i in 1:length(altibandMin)){
  
  print(i)
  
  petMat.bryo <- bryoData[alti>=altibandMin[i] & alti<altibandMax[i],3:324]
  petMat.bryo <- petMat.bryo[,apply(petMat.bryo,2,sum)>0]
  ToKeep <- apply(petMat.bryo,1,sum)>0
  petMat.bryo <- petMat.bryo[ToKeep,]
  petMat.bryo <- beta.pair(petMat.bryo)$beta.sim
  
  TO <- petMat.bryo[lower.tri(as.matrix(petMat.bryo),diag=F)]
  fin1 <- rbind.data.frame(fin1,cbind.data.frame(Taxon="Bryophytes",
                                                 TO = TO,
                                                 elevationalBand = paste0("[",altibandMin[i],";", altibandMax[i],"[")))
  
  petMat.bryo.env <- bryoEnv[alti>=altibandMin[i] & alti<altibandMax[i],]
  petMat.bryo.env <- petMat.bryo.env[ToKeep,]
  petMat.bryo.env <- petMat.bryo.env[,apply(petMat.bryo.env,2,sum)!=0]
  petMat.bryo.env <- dist(scale(petMat.bryo.env))
  
  petMat.tracheo.traits <- tracheoTraits[alti>=altibandMin[i] & alti<altibandMax[i],]
  petMat.tracheo.traits <- petMat.tracheo.traits[ToKeep,]
  petMat.tracheo.traits <- petMat.tracheo.traits[,apply(petMat.tracheo.traits,2,sum)>0]
  petMat.tracheo.traits <- dist(petMat.tracheo.traits)
  bEnv <- mantel.partial(petMat.bryo,petMat.bryo.env,petMat.tracheo.traits)
  bTraits <- mantel.partial(petMat.bryo,petMat.tracheo.traits,petMat.bryo.env)
  
  inter <- cbind.data.frame(dataset="bryophytes", 
                            envR= bEnv$statistic,
                            envSigni = bEnv$signif,
                            traitsR = bTraits$statistic,
                            traitsSigni = bTraits$signif,
                            elevationalBand = paste0("[",altibandMin[i],";", altibandMax[i],"["))
  
  fin2 <- rbind.data.frame(fin2,inter)
  
  
  
  petMat.tracheo <- TracheoData[alti>=altibandMin[i] & alti<altibandMax[i],3:739]
  petMat.tracheo <- petMat.tracheo[,apply(petMat.tracheo,2,sum)>0]
  ToKeep <- apply(petMat.tracheo,1,sum)>0
  petMat.tracheo <- petMat.tracheo[ToKeep,]
  petMat.tracheo <- beta.pair(petMat.tracheo)$beta.sim
  
  TO <- petMat.tracheo[lower.tri(as.matrix(petMat.tracheo),diag=F)]
  fin1 <- rbind.data.frame(fin1,cbind.data.frame(Taxon="Tracheophytes",
                                                 TO = TO,
                                                 elevationalBand = paste0("[",altibandMin[i],";", altibandMax[i],"[")))
  
  
  petMat.tracheo.env <- TracheoEnv[alti>=altibandMin[i] & alti<altibandMax[i],]
  petMat.tracheo.env <- petMat.tracheo.env[ToKeep,]
  petMat.tracheo.env <- petMat.tracheo.env[,apply(petMat.tracheo.env,2,sum)!=0]
  petMat.tracheo.env <- dist(scale(petMat.tracheo.env))
  
  petMat.bryo.traits <- bryoTraits[alti>=altibandMin[i] & alti<altibandMax[i],]
  petMat.bryo.traits <- petMat.bryo.traits[ToKeep,]
  petMat.bryo.traits <- petMat.bryo.traits[,apply(petMat.bryo.traits,2,sum)>0]
  petMat.bryo.traits <- dist(petMat.bryo.traits)
  
  bEnv <- mantel.partial(petMat.tracheo,petMat.tracheo.env,petMat.bryo.traits)
  bTraits <- mantel.partial(petMat.tracheo,petMat.bryo.traits,petMat.tracheo.env)
  
  inter <- cbind.data.frame(dataset="tracheophytes", 
                            envR= bEnv$statistic,
                            envSigni = bEnv$signif,
                            traitsR = bTraits$statistic,
                            traitsSigni = bTraits$signif,
                            elevationalBand = paste0("[",altibandMin[i],";", altibandMax[i],"["))
  
  fin2 <- rbind.data.frame(fin2,inter)
  
  
}

write.csv2(fin2,"MantelPartial.csv")
write.csv2(fin1,"TO_alti.csv")
fin1$elevationalBand = as.factor(fin1$elevationalBand)
levels(fin1$elevationalBand)[6] = ">= 2200m"
levels(fin1$elevationalBand)[1]  = "< 1000 m"
png("TO.png",res=500,width=6000,height = 4500)
p <- ggplot(fin1,aes(y=TO,x=elevationalBand,fill=Taxon))+geom_boxplot()+scale_fill_discrete(type=c("#B68F68","#56876D"))+ylab("Turnover")+xlab("Elevation Band (m)")+
  theme_classic()+
  theme(axis.text.x = element_text(size=15,angle = 30, hjust = 1),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title =element_text(size=15),
        legend.position = "bottom")
p
dev.off()


fin2 <- NULL
for(i in 1:length(altibandMin)){
  
  print(i)
  
  petMat.bryo <- bryoData[alti>=altibandMin[i] & alti<altibandMax[i],3:324]
  petMat.bryo <- petMat.bryo[,apply(petMat.bryo,2,sum)>0]
  
  petMat.tracheo <- TracheoData[alti>=altibandMin[i] & alti<altibandMax[i],3:739]
  petMat.tracheo <- petMat.tracheo[,apply(petMat.tracheo,2,sum)>0]
  
  ToKeep <- apply(petMat.bryo,1,sum)>0 & apply(petMat.tracheo,1,sum)>0
  petMat.bryo <- petMat.bryo[ToKeep,]
  
  petMat.bryo <- beta.pair(petMat.bryo)$beta.sim
  
  
  petMat.tracheo <- petMat.tracheo[ToKeep,]
  petMat.tracheo <- beta.pair(petMat.tracheo)$beta.sim
  
  bEnv <- mantel(petMat.tracheo,petMat.bryo)
  
  inter <- cbind.data.frame(tracheoBryo.R= bEnv$statistic,
                            tracheoBryo.Signi = bEnv$signif,
                            elevationalBand = paste0("[",altibandMin[i],";", altibandMax[i],"["),
                            nPlots = sum(ToKeep))
  
  fin2 <- rbind.data.frame(fin2,inter)
  
  
}

write.csv2(fin2,"MantelPartial_TracheoBryo.csv")





#################
## Final plot for the figure ----
library(ggplot2)
library(viridis)
bryoData <- read.csv("BryophyteData.csv",row.names=1)
sum(apply(bryoData[,3:324],2,sum)>0)

bryoEnv <- read.csv("EnvDataForBryophytes.csv",row.names=1)
alti <- bryoEnv$mnt_mean

fin <- cbind.data.frame(ElevationalBand = factor(c("<1000m","[1000m ; 1400m[","[1400m ; 1800m[","[1800m ; 2000m[","[2000m ; 2200m[",">2200m"),
                                                 levels = c("<1000m","[1000m ; 1400m[","[1400m ; 1800m[","[1800m ; 2000m[","[2000m ; 2200m[",">2200m")),
                        NPlots = c(sum(alti<1000),sum(alti<1400 & alti>=1000),sum(alti<1800 & alti>=1400),sum(alti<2000 & alti>=1800),
                                   sum(alti<2200 & alti>=2000),sum(alti>=2200)))
png("Fig1_Part2.png",res=500,width=2500,height = 3000)
p <- ggplot(fin,aes(x=ElevationalBand, y=NPlots, fill= ElevationalBand, colour= ElevationalBand)) + geom_col() + 
  xlab("") + ylab("Number of plots") + 
  scale_fill_viridis_d(option="inferno")+ scale_color_discrete(type=rep("black",6))+
  theme_classic()+
  theme(axis.text.x = element_text(size=15,angle = 30, hjust = 1),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.position = "none")
p
dev.off()

#### Prioritization ----
bryoData <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/BryophyteData.csv",row.names=1)

bryoEnv <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/EnvDataForBryophytes.csv",row.names=1)
alti <- bryoEnv$mnt_mean

TracheoData <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/TracheophyteData.csv",row.names=1)


# Bryophytes ----
BryoData <- bryoData[,5:ncol(bryoData)]
PU <- cbind.data.frame(id = 1:413,cost = 1)
features = cbind.data.frame(id = 1:ncol(BryoData),name = colnames(BryoData))
#rij should be a super long matrix with sp identifier + 1 when present and 0 when absent in the plot
rij <- NULL
for(i in 1: 413){
  inter <- cbind.data.frame(pu=i,species = 1:ncol(BryoData), amount = as.numeric(BryoData[i,]))
  rij <- rbind(rij,inter)
}

p0 <- problem(PU, features = features, rij = rij,cost_column = "cost")%>%
  add_min_set_objective() %>%
  add_relative_targets(0.3)%>%
  add_binary_decisions() %>%
  add_default_solver(verbose=FALSE)

s0 <- solve(p0)
sum(s0$solution_1) #107


# Tracheo ----
TracheoData <- TracheoData[,3:ncol(TracheoData)]
PU <- cbind.data.frame(id = 1:413,cost = 1)
features = cbind.data.frame(id = 1:ncol(TracheoData),name = colnames(TracheoData))
#rij should be a super long matrix with sp identifier + 1 when present and 0 when absent in the plot
rij <- NULL
for(i in 1: 413){
  inter <- cbind.data.frame(pu=i,species = 1:ncol(TracheoData), amount = as.numeric(TracheoData[i,]))
  rij <- rbind(rij,inter)
}
p1 <- problem(PU, features = features, rij = rij,cost_column = "cost")%>%
  add_min_set_objective() %>%
  add_relative_targets(0.3)%>%
  add_binary_decisions() %>%
  add_default_solver(verbose=FALSE)

s1 <- solve(p1)
sum(s1$solution_1) #163


sum(s0$solution_1 ==  1 & s1$solution_1 == 1)/413
sum(s0$solution_1 ==  0 & s1$solution_1 == 1)/413
sum(s0$solution_1 ==  1 & s1$solution_1 == 0)/413

elevPrio.Bryo <- alti[s0$solution_1==1]
elevPrio.Tracheo <- alti[s1$solution_1==1]
altibandMin <- c(0,1000,1400,1800,2000,2200)
altibandMax <- c(1000,1400,1800,2000,2200,Inf)
test <-2*s0$solution_1 + s1$solution_1

fin <- NULL
for(i in 1:length(altibandMax)){
  ToKeep <- test == 3 & alti>=altibandMin[i] & alti<altibandMax[i]
  
  inter <- cbind.data.frame(elevationalBand = paste0("[",altibandMin[i],";", altibandMax[i],"["),
                            NPlot = sum(alti>=altibandMin[i] & alti<altibandMax[i]),
                            OnlyBryo = sum(test == 2 & alti>=altibandMin[i] & alti<altibandMax[i]),
                            OnlyTracheo = sum(test == 1 & alti>=altibandMin[i] & alti<altibandMax[i]),
                            Both = sum(test == 3 & alti>=altibandMin[i] & alti<altibandMax[i]))
  fin <- rbind(fin,inter)
}
fin$TotBryo = fin$OnlyBryo + fin$Both
fin$TotTracheo = fin$OnlyTracheo + fin$Both
fin$BothPerc = 100*fin$Both/fin$NPlot
write.csv(fin,"PrioritisationV1.csv")
write.csv(s0$solution_1,"PrioritisationV1_Bryo.csv")
write.csv(s1$solution_1,"PrioritisationV1_Tracheo.csv")


##### Prioritisation V2 ----

bryoData <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/BryophyteData.csv",row.names=1)

bryoEnv <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/EnvDataForBryophytes.csv",row.names=1)
alti <- bryoEnv$mnt_mean

TracheoData <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/TracheophyteData.csv",row.names=1)


# Bryophytes ----
BryoData <- bryoData[,5:ncol(bryoData)]
PU <- cbind.data.frame(id = 1:413,cost = 1)
features = cbind.data.frame(id = 1:ncol(BryoData),name = colnames(BryoData))
#rij should be a super long matrix with sp identifier + 1 when present and 0 when absent in the plot
rij <- NULL
for(i in 1: 413){
  inter <- cbind.data.frame(pu=i,species = 1:ncol(BryoData), amount = as.numeric(BryoData[i,]))
  rij <- rbind(rij,inter)
}
budget = 124
p0 <- problem(PU, features = features, rij = rij,cost_column = "cost")%>%
  add_max_features_objective(budget = budget) %>%
  add_relative_targets(0.3) %>%
  add_binary_decisions() %>%
  add_linear_constraints(threshold = 124, "=", PU$cost) %>%
  add_default_solver(verbose=FALSE)

s0 <- solve(p0)
sum(s0$solution_1) #107


# Tracheo ----
TracheoData <- TracheoData[,3:ncol(TracheoData)]
PU <- cbind.data.frame(id = 1:413,cost = 1)
features = cbind.data.frame(id = 1:ncol(TracheoData),name = colnames(TracheoData))
#rij should be a super long matrix with sp identifier + 1 when present and 0 when absent in the plot
rij <- NULL
for(i in 1: 413){
  inter <- cbind.data.frame(pu=i,species = 1:ncol(TracheoData), amount = as.numeric(TracheoData[i,]))
  rij <- rbind(rij,inter)
}
p1 <- problem(PU, features = features, rij = rij,cost_column = "cost")%>%
  add_max_features_objective(budget = budget) %>%
  add_relative_targets(0.3) %>%
  add_binary_decisions() %>%
  add_linear_constraints(threshold = 124, "=", PU$cost) %>%
  add_default_solver(verbose=FALSE)


s1 <- solve(p1)
sum(s1$solution_1) #163


sum(s0$solution_1 ==  1 & s1$solution_1 == 1)/413
sum(s0$solution_1 ==  0 & s1$solution_1 == 1)/413
sum(s0$solution_1 ==  1 & s1$solution_1 == 0)/413

elevPrio.Bryo <- alti[s0$solution_1==1]
elevPrio.Tracheo <- alti[s1$solution_1==1]
altibandMin <- c(0,1000,1400,1800,2000,2200)
altibandMax <- c(1000,1400,1800,2000,2200,Inf)
test <-2*s0$solution_1 + s1$solution_1

fin <- NULL
for(i in 1:length(altibandMax)){
  ToKeep <- test == 3 & alti>=altibandMin[i] & alti<altibandMax[i]
  
  inter <- cbind.data.frame(elevationalBand = paste0("[",altibandMin[i],";", altibandMax[i],"["),
                            NPlot = sum(alti>=altibandMin[i] & alti<altibandMax[i]),
                            OnlyBryo = sum(test == 2 & alti>=altibandMin[i] & alti<altibandMax[i]),
                            OnlyTracheo = sum(test == 1 & alti>=altibandMin[i] & alti<altibandMax[i]),
                            Both = sum(test == 3 & alti>=altibandMin[i] & alti<altibandMax[i]))
  fin <- rbind(fin,inter)
}
fin$TotBryo = fin$OnlyBryo + fin$Both
fin$TotTracheo = fin$OnlyTracheo + fin$Both
fin$BothPerc = 100*fin$Both/fin$NPlot
write.csv(fin,"PrioritisationV2.csv")
write.csv(s0$solution_1,"PrioritisationV2_Bryo.csv")
write.csv(s1$solution_1,"PrioritisationV2_Tracheo.csv")


## Check if Prioritisation Match Li Values ----

sol.Bryo <- read.csv("CorrelationsBetweenSRTracheoAndBryo/PrioritisationV2_Bryo.csv",row.names = 1)
sol.Tracheo <- read.csv("CorrelationsBetweenSRTracheoAndBryo/PrioritisationV2_Tracheo.csv",row.names = 1)

bryoData <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/BryophyteData.csv",row.names=1)
SR.bryo <- apply(bryoData[,5:ncol(bryoData)],1,sum)
bryoEnv <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/EnvDataForBryophytes.csv",row.names=1)
alti <- bryoEnv$mnt_mean
TracheoData <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/TracheophyteData.csv",row.names=1)
SR.tracheo <- apply(TracheoData[,3:ncol(TracheoData)],1,sum)

nbg <- graph2nb(gabrielneigh(bryoData[,c("x","y")]), sym = TRUE) #Spatial links
listw <- nb2listw(nbg)
Li <- lee(SR.bryo,SR.tracheo,listw,n=413) #Global = 0.1275

Common <- sol.Tracheo+sol.Bryo

boxplot(Li$localL[Common==2], Li$localL[Common==1], Li$localL[Common==0], names=c("Both","Prio1","None"),ylab = "Li")

sol.Bryo <- read.csv("CorrelationsBetweenSRTracheoAndBryo/PrioritisationV1_Bryo.csv",row.names = 1)
sol.Tracheo <- read.csv("CorrelationsBetweenSRTracheoAndBryo/PrioritisationV1_Tracheo.csv",row.names = 1)
Common <- sol.Tracheo+sol.Bryo
boxplot(Li$localL[Common==2], Li$localL[Common==1], Li$localL[Common==0], names=c("Both","Prio1","None"),ylab = "Li")


## Beta Diversity analyses with GDM ----
bryoData <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/BryophyteData.csv",row.names=1)
bryoSp <- bryoData[,5:ncol(bryoData)]


TracheoData <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/TracheophyteData.csv",row.names=1)
TracheoSp <- TracheoData[,3:ncol(TracheoData)]

ToKeep <- apply(bryoSp, 1, sum)>0 & apply(TracheoSp, 1, sum)>0 
bryoSp <- bryoSp[ToKeep ,]
TO.bryo <- as.dist(beta.pair(bryoSp)$beta.sim	)
TracheoSp <- TracheoSp[ToKeep ,]
TO.tracheo <- as.dist(beta.pair(TracheoSp)$beta.sim)

mantel(TO.bryo,TO.tracheo)
nbg <- graph2nb(gabrielneigh(bryoData[ToKeep,c("x","y")]), sym = TRUE) #Spatial links
listw <- nb2listw(nbg)
library(adespatial )


# Mantel <- mantel.randtest(TO.tracheo,TO.bryo)
# rSEnviTrop <- Mantel$obs
# signiSEnviTrop <- Mantel$pvalue
# result<- msr(x=Mantel,listw,999)
# 
# 
# rEnviCorTrop <- result$obs-result$expvar["Expectation"]
# signiEnviCorTrop <- result$pvalue

library(gdm)

TO.bryo.gdm <- cbind(Sites = 1:sum(ToKeep), as.matrix(TO.bryo))
bryoEnv <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/EnvDataForBryophytes.csv",row.names=1)
bryoEnv <- bryoEnv[ToKeep,]

env <- dudi.pca(bryoEnv,scannf = F,nf=ncol(bryoEnv)-1,scale=T)
NAxesToKeep <- which(cumsum(100*env$eig / sum(env$eig )) >90)[1]
env <- env$li[,1:NAxesToKeep]
env <- dist(env)
env <- cbind(Sites = 1:sum(ToKeep), as.matrix(env))

tracheoTraits <- read.csv("DatasetsToUse_Cleaned/paper_STOTEN/Traits_Tracheophytes_Fraction.csv",row.names = 1)


tracheoTraits <- dudi.pca(tracheoTraits[ToKeep,],scannf = F,nf=ncol(tracheoTraits)-1,scale=T)
NAxesToKeep <- which(cumsum(100*tracheoTraits$eig / sum(tracheoTraits$eig )) >90)[1]
tracheoTraits <- tracheoTraits$li[,1:NAxesToKeep]
tracheoTraits <- dist(tracheoTraits)
tracheoTraits <- cbind(Sites = 1:sum(ToKeep), as.matrix(tracheoTraits))

coord <- cbind(Sites = 1:sum(ToKeep),bryoData[ToKeep,3:4])

TO.tracheo.gdm <- cbind(Sites = 1:sum(ToKeep), as.matrix(TO.tracheo))
FormatGDM <- formatsitepair(bioData = TO.bryo.gdm,
                            bioFormat = 3,
                            XColumn="x",
                            YColumn="y",
                            predData=coord,
                            siteColumn="Sites",
                            distPreds = list(env = env,
                                             Tracheo = TO.tracheo.gdm,
                                             Traits = tracheoTraits))


gdm.1 <- gdm.varImp(FormatGDM, geo=T,
                    nPerm = 100,
                    parallel = T,
                    cores = 5,
                    predSelect = T) #geo = T if you want to take spatial autocorrelation into account
summary(gdm.1)
gdm.1$`Predictor Importance`
# Geographic   matrix_1   matrix_2 
# 2.085      7.530     37.092 
barplot(sort(gdm.1$`Predictor Importance`, decreasing=T))
