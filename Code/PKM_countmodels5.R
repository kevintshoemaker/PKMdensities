
# PKM analysis- count models 

# QUESTIONS  ----------------------


# TODO -------------------------


# TIMELINE  --------------------


# submit summer 2024


# clear workspace --------------------

rm(list=ls())

# load packages ----------------------

library(glmmTMB)
library(DHARMa)
library(car)
library(lme4)
library(buildmer)
#library(broom.mixed)
library(MuMIn)
library(rsq)
library(caret)
#library(glmmLasso)
#library(glmnet)
library(MASS)
library(MuMIn)

# library(sp)
# library(rgeos)
# library(raster)

library(sf)
library(terra)

library(lunar)


# load functions ------------------------------

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

VisualizeRelation <- function(data=master_df2,model=bestmod3,predvar=topvars[3],allvars=allvars){
  len <- 100
  
  predvar2 <- gsub("_std","",predvar)
  
  predvar3 <- gsub("_[[:digit:]][[:digit:]]nn", "",  predvar2)
  predvar4 <- gsub("_[[:digit:]]nn", "",  predvar2)
  predvar5 <- gsub("_[[:digit:]][[:digit:]][[:digit:]]nn", "",  predvar2)
  
  xscal <- as.numeric(gsub("[^[:digit:]]","",predvar))
  
  dataclasses <- sapply(data,class)
  
  dim <- data[,predvar]
  range <- seq(min(dim,na.rm=T),max(dim,na.rm=T),length=len)
  
  realmean <- mean(data[[predvar2]])
  realsd <- sd(data[[predvar2]])
  
  newdata <- data.frame(temp=range)
  names(newdata) <- c(predvar)
  
  othervars <- allvars[!allvars%in%c(predvar,"totcaps")]
  
  
  var = othervars[5]
  for(var in othervars){
    thisvar <- data[[var]]
    if(is.factor(thisvar)){
      tab <- table(thisvar)
      vals <- names(tab)
      levs <- levels(thisvar)
      mostcom <- vals[which.max(tab)]
      newvec <- factor(rep(mostcom,times=nrow(newdata)),levels=levs)
      newdata[,var] <- newvec
    }else{
      newdata[,var] <- 0 #mean(thisvar)
    }
  }
  
  
  pred <- predict(model,newdata,type="response",se.fit=T,allow.new.levels = T)
  
  plot(range,pred$fit,xlab=sprintf("%s (%s m)",covar_df$clearname[covar_df$origname%in%c(predvar3,predvar4,predvar5)],xscal),
       ylab="Exp Count",type="l",lwd=2,xaxt="n",ylim=c(0,3))
  points(range,pred$fit+pred$se.fit,type="l",lty=2)
  points(range,pred$fit-pred$se.fit,type="l",lty=2)
  ats <- seq(min(range),max(range),length=6)
  axis(1,ats,labels = round(realmean+ats*realsd,1))
  rug(jitter(data[seq(1,nrow(data),5),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  
  # svg(sprintf("pdplots_%s.svg",predvar2),width=4.5,height = 4)
  # 
  # plot(range,pred$fit,xlab=covar_df$clearname[covar_df$origname==predvar2],ylab="Exp Count",type="l",lwd=2,xaxt="n",ylim=c(0,3))
  # points(range,pred$fit+pred$se.fit,type="l",lty=2)
  # points(range,pred$fit-pred$se.fit,type="l",lty=2)
  # ats <- seq(min(range),max(range),length=6)
  # axis(1,ats,labels = round(realmean+ats*realsd,2))
  # rug(jitter(data[seq(1,nrow(data),50),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  # 
  # dev.off()
}

VisualizeInteraction <- function(data=master_df,model=bestmod3,var1=ints[[1]][1],var2=ints[[1]][2],allvars=allvars){
  len <- 30     # increase this for higher-res figures
  
  dataclasses <- sapply(data,class)
  
  var1_2 <- gsub("_std","",var1)
  var2_2 <- gsub("_std","",var2)
  
  temp3 <- gsub("_[[:digit:]]nn", "",  var1_2)
  temp4 <- gsub("_[[:digit:]][[:digit:]]nn", "",  var1_2)
  temp5 <- gsub("_[[:digit:]][[:digit:]][[:digit:]]nn", "",  var1_2)
  
  temp6 <- gsub("_[[:digit:]]nn", "",  var2_2)
  temp7 <- gsub("_[[:digit:]][[:digit:]]nn", "",  var2_2)
  temp8 <- gsub("_[[:digit:]][[:digit:]][[:digit:]]nn", "",  var2_2)
  
  var1_2 <- c(temp3,temp4,temp5)[which.min(c(nchar(temp3),nchar(temp4),nchar(temp5)))]
  var2_2 <- c(temp6,temp7,temp8)[which.min(c(nchar(temp6),nchar(temp7),nchar(temp8)))]
  
  xscal <- as.numeric(gsub("[^[:digit:]]","",var1))
  yscal <- as.numeric(gsub("[^[:digit:]]","",var2))
  
  realmean1 <- mean(data[[var1_2]],na.rm=T)
  realsd1 <- sd(data[[var1_2]],na.rm=T)
  realmean2 <- mean(data[[var2_2]],na.rm=T)
  realsd2 <- sd(data[[var2_2]],na.rm=T)
  standvar1 <- var1
  standvar2 <- var2
  
  
  firstdim <- data[,standvar1]
  seconddim <- data[,standvar2]
  range1 <- seq(min(firstdim,na.rm=T),max(firstdim,na.rm=T),length=len)
  range2 <- seq(min(seconddim,na.rm=T),max(seconddim,na.rm=T),length=len)
  newdata <- expand.grid(range1,range2)
  # head(newdata,50)
  names(newdata) <- c(standvar1,standvar2)
  
  othervars <- allvars[!allvars%in%c(standvar1,standvar2,"totcaps")]
  
  var = othervars[2]
  for(var in othervars){
    thisvar <- data[[var]]
    if(is.factor(thisvar)){
      tab <- table(thisvar)
      vals <- names(tab)
      levs <- levels(thisvar)
      mostcom <- vals[which.max(tab)]
      newvec <- factor(rep(mostcom,times=nrow(newdata)),levels=levs)
      newdata[[var]] <- newvec
    }else{
      newdata[[var]] <- mean(thisvar,na.rm=T)
    }
  }
  
  newdata$Gd2<-1
  pred <- predict(model,newdata,type="response")
  
  predmat <-  matrix(pred,nrow=len,ncol=len)
  
  predmat[which(predmat>4,arr.ind = T)] = 4
  
  
  
  #svg(sprintf("intplots_%s_%s.svg",var1_2,var2_2),width=4.5,height = 4)
  
  # layout(matrix(c(1:1),nrow=1))
  par(mai=c(.3,.5,.1,.1))
  
  persp(realmean1+realsd1*range1,realmean2+realsd2*range2,predmat,theta = 55, phi = 40, r = sqrt(10), d = 3, 
        ticktype = "detailed", mgp = c(4, 1, 0),
        zlab="\nExp count",
        xlab=sprintf("\n%s (%s m)",covar_df$clearname[covar_df$origname==var1_2],xscal),
        ylab=sprintf("\n%s (%s m)",covar_df$clearname[covar_df$origname==var2_2],yscal))-> res
  
  points(trans3d(data[[var1_2]], 
                 data[[var2_2]], max(predmat)+0, pmat = res), col = gray(0.8), pch = 19, cex=0.5)
  # points(trans3d(jitter(data[[var1_2]],range(data[[var1_2]])), 
  #                jitter(data[[var2_2]],range(data[[var2_2]])), max(predmat)+0, pmat = res), col = gray(0.8), pch = 19, cex=0.5)
  # 
  
  #dev.off() 
}

# load global vars --------------------

global <- list()

global$CRS <- 26911 #CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

global$BaseDir <- getwd()
global$NAIP_Raster_Dir <- "M:\\GIS\\NAIP\\predictions"


# load and process data ----------------

dat <- list()

## trap locations -----------------

dat$traplocs <- read.csv("trap_locations.csv",header=T)

dat$traplocs$Det <- dat$traplocs$X.Det
dat$traplocs$X.Det <- NULL
dat$traplocs$X <- NULL
dat$traplocs$Gd <- dat$traplocs$GdNo
dat$traplocs$GdNo <- NULL

dat$traplocs <- st_as_sf(dat$traplocs,coords = c("Easting","Northing"),crs=global$CRS)

rownames(dat$traplocs) <- dat$traplocs$Det

plot(dat$traplocs$geometry)

dat$distmat <- as.matrix(dist(st_coordinates(dat$traplocs)))
nrow(dat$distmat)
ncol(dat$distmat)
nrow(dat$traplocs)
rownames(dat$distmat) <- as.character(dat$traplocs$Det)
colnames(dat$distmat) <- as.character(dat$traplocs$Det)

## Trap covariates ----------------

dat$trapcovars <- read.csv("TrapCovars.csv",header=T)
summary(dat$trapcovars)
head(dat$trapcovars)
names(dat$trapcovars)

dat$trapcovars$X.Det

names(dat$trapcovars) <- gsub("X.","",names(dat$trapcovars))

rownames(dat$trapcovars) <- dat$trapcovars$Det

length(unique(dat$trapcovars$Det)) == nrow(dat$trapcovars)   # true

global$alldets <- sort(unique(dat$traplocs$X.Det))
length(global$alldets)  # 1152 unique detector locations

global$allgds <- sort(unique(dat$trapcovars$Gd))
# global$allgds

dat$trapcovars$Q_DnCnt <- dat$trapcovars$Q_SDnCnt+dat$trapcovars$Q_MDnCnt

names(dat$trapcovars)

covar_df <- data.frame(
  origname = names(dat$trapcovars)[-c(1,2)],
  std = sprintf("%s_std",names(dat$trapcovars)[-c(1,2)]),
  clearname = c("BunchGrass","BareGround (%)","Cobble","Litter","CoarseGravel","Gravel","SmallDunes","MiniDunes",
                "Grass","Forbs","ShrubCover","Dist_NearestShrub","ShrubCanopy","ShrubHeight (cm)","ShrubCount",
                "ShrubSR","GrassForb","Dunes"),
  stringsAsFactors = F
)
covar_df

t=covar_df$origname[3]
covar_df$mean <- sapply(covar_df$origname,function(t) mean(dat$trapcovars[[t]],na.rm=T) )
covar_df$min <- sapply(covar_df$origname,function(t) min(dat$trapcovars[[t]],na.rm=T) )
covar_df$max <- sapply(covar_df$origname,function(t) max(dat$trapcovars[[t]],na.rm=T) )
covar_df$sd <- sapply(covar_df$origname,function(t) sd(dat$trapcovars[[t]],na.rm=T) )

temp <- sapply(1:length(covar_df$origname), function(t)   # standardize all covars
  dat$trapcovars[[covar_df$std[t]]] <<- (dat$trapcovars[[covar_df$origname[t]]]-covar_df$mean[t])/covar_df$sd[t] )


## Dates ----------------------

dat$dates <- read.csv("dates.csv",header=T)

## Soils ---------------------

# TODO: Revise with new categories- determined with Marjorie and Sarah 11/16/23
#   Marjorie sent new datafile for soils... 
#   Sarah sent updated datafile for soils...
#   if zero sum constraint, maybe leave out silt/clay...

# dat$soils_orig <- read.csv("Soils.csv",header=T)   # original summarized soils data
# summary(dat$soils_orig)

# alldet_orig <- sort(unique(dat$soils_orig$Det)) 
# length(alldet_orig)
# setdiff(alldet_orig,global$alldets) #same

dat$soilsraw <- read.csv("Soils2_Wblanks_15Jan24.csv",header=T)  # raw soils data
summary(dat$soilsraw)

# nrow(dat$soils_orig)  # 1152 observations  
nrow(dat$soilsraw)       # also 1152 rows

# dat$soils_orig

# dat$soils_orig$Det
# setdiff(dat$soils_orig$Det,dat$soils$UnqID)  # these are not the same

temp <- dat$soilsraw
dat$soilsraw$UnqID <- tolower(dat$soilsraw$UnqID)
uncertain <- which(grepl("xx",dat$soilsraw$UnqID))
dat$soilsraw$Det <- gsub("-","_",dat$soilsraw$Gd.St)
dat$soilsraw$Det <- as.numeric(gsub("_","",dat$soilsraw$Det))
alldets <- sort(unique(dat$soilsraw$Det))
length(alldets)  #1150 (2 fewer than total detectors)
duplicated <- which(duplicated(dat$soilsraw$Det))  # two duplicates

dup_test <- dat$soils[dat$soilsraw$Det%in%dat$soilsraw$Det[duplicated],]
dup_test <- dup_test[order(dup_test$Det),]
# write.csv(dup_test,"duplicated_soils.csv",row.names = F)

missing <- setdiff(dat$trapcovars$Det,alldets)  # a couple missing...
length(missing)
names(dat$soilsraw)

# nrow(dat$traplocs) # same as trap locations

# ndx <- match(as.numeric(row.names(dat$trapcovars)),dat$soils$Det)  # closer?? 

# cbind(row.names(dat$trapcovars),dat$soils$Det[ndx])    # check

### Sarah: combine silt and clay into a "fine" category-- so coarse, sand and fine... 

soilcats <- read.csv("soilcats.csv")
soil_varnames <- setdiff(names(soilcats),c("X","UnqID","Grid","Station","Gd.St","Det","Total","UnqID.save"))

names(dat$soilsraw)

dat$soilsraw[uncertain,soil_varnames][] <- NA    # set uncertain values to NA

hasdata <- apply(dat$soilsraw[,soil_varnames],1,function(t) !all(is.na(t)) )
dat$soilsraw <- dat$soilsraw[hasdata,]
nrow(dat$soilsraw)  #1104 data points vs 1152 for full dataset.  Consider interpolating by location?

# dat$soils$Fine <- dat$soils$Clay+dat$soils$Silt
# hist(dat$soils$Fine)


### explore whether nearby stations have similar soil characteristics (if so, use interpolation for missing data)

# dat$distmat[1,]
nearest3 <- t(sapply(dat$traplocs$Det,function(t) dat$traplocs$Det[order(dat$distmat[as.character(t),])[1:3]]  ))
dim(nearest3)
rownames(nearest3) <- dat$traplocs$Det

near3_raw <- nearest3[match(dat$soilsraw$Det,rownames(nearest3)),]
nrow(near3_raw)
v=1
for(v in 1:length(soil_varnames)){
  thisvar <- soil_varnames[v]
  dat$soilsraw[[paste0(thisvar,"_near3")]] <- sapply(1:nrow(near3_raw),function(t) mean(dat$soilsraw[[thisvar]][dat$soilsraw$Det%in%near3_raw[t,1:3]],na.rm=T)   )
}

soil_varnames
lm(cobble~cobble_near3,dat$soilsraw)  # near3 is very good predictor of cobble
plot(cobble~cobble_near3,dat$soilsraw)

lm(gravel~gravel_near3,dat$soilsraw)  # near3 is very good predictor of cobble
plot(gravel~gravel_near3,dat$soilsraw)

lm(coarsesand~coarsesand_near3,dat$soilsraw)  # near3 is very good predictor of cobble
plot(coarsesand~coarsesand_near3,dat$soilsraw)

lm(sand3~sand3_near3,dat$soilsraw)  # near3 is very good predictor of cobble
plot(sand3~sand3_near3,dat$soilsraw)

summary(lm(sand2~sand2_near3,dat$soilsraw))  # near3 is very good predictor of cobble
plot(sand2~sand2_near3,dat$soilsraw)

# do interpolation of soils for missing data

nrow(dat$soilsraw)
length(unique(dat$soilsraw$Det))  # no duplicated data left

missing <- setdiff(dat$traplocs$Det,dat$soilsraw$Det)
length(missing)   # 48 data points to interpolate...

temp <- data.frame(Det=missing)
near3_miss <- nearest3[match(temp$Det,rownames(nearest3)),]
nrow(near3_miss)
v=1
for(v in 1:length(soil_varnames)){
  thisvar <- soil_varnames[v]
  temp[[thisvar]] <- sapply(1:nrow(near3_miss),function(t) mean(dat$soilsraw[[thisvar]][dat$soilsraw$Det%in%near3_miss[t,1:3]],na.rm=T)   )
}

thesecols <- match(names(temp),names(dat$soilsraw))

temp2 <- data.frame(UnqID = rep(NA,times=length(missing)))
temp3 <- lapply(names(dat$soilsraw)[-1],function(t) temp2[[t]] <<- NA )
temp2[,thesecols] <- temp

ncol(temp2)
ncol(dat$soilsraw)

if(nrow(dat$soilsraw)<nrow(dat$traplocs)){
  dat$soilsraw <- rbind(dat$soilsraw,temp2)
}
nrow(dat$soilsraw)

## make new soils dataset with new categories

# library(tidyverse)

soilcats$X

soilcats2 <- soilcats[soilcats$X%in%c("medfinevf","allusda","coarsemed","orig","medfine","medium","nottoofine","finevf","fine"),]
rownames(soilcats2) <- c("medfinevf","allusda","coarsemed","orig","medfine","medium","nottoofine","finevf","fine")
soilcats2$X <- NULL
all(colnames(soilcats2)==soil_varnames)

temp <- data.frame(Det=dat$traplocs$Det)
nrow(temp)  

# soilcats2
i=1
for(i in 1:nrow(soilcats2)){
  thisvar <- rownames(soilcats2)[i]
  thesecols <- soil_varnames[which((soilcats2[rownames(soilcats2)==rownames(soilcats2)[i],]==1))]
  # thesecols_ndx <- match(thesecols,names(dat$soilsraw))
  # t=dat$soilsraw[1,]
  temp[[thisvar]] <- sapply(1:nrow(dat$soilsraw),function(t) sum(dat$soilsraw[t,thesecols]) )
}

temp2 <- dat$soilsraw[,soil_varnames]

dat$soils <- cbind(temp,temp2)

covar_df2 <- data.frame(
  origname = names(dat$soils)[-c(1)],
  std = sprintf("%s_std",names(dat$soils)[-c(1)]),
  clearname = c(rownames(soilcats2),soil_varnames),
  stringsAsFactors = F
)
covar_df2

# t=covar_df2$origname[3]
covar_df2$mean <- sapply(covar_df2$origname,function(t) mean(dat$soils[[t]],na.rm=T) )
covar_df2$min <- sapply(covar_df2$origname,function(t) min(dat$soils[[t]],na.rm=T) )
covar_df2$max <- sapply(covar_df2$origname,function(t) max(dat$soils[[t]],na.rm=T) )
covar_df2$sd <- sapply(covar_df2$origname,function(t) sd(dat$soils[[t]],na.rm=T) )

temp <- sapply(1:length(covar_df2$origname), function(t)   # standardize all covars
  dat$soils[[covar_df2$std[t]]] <<- (dat$soils[[covar_df2$origname[t]]]-covar_df2$mean[t])/covar_df2$sd[t] )

ndx <- match(dat$trapcovars$Det,dat$soils$Det)
# cbind(dat$trapcovars$Det, dat$soils$Det[ndx])
master_df <- cbind(dat$trapcovars,dat$soils[ndx,])

covar_df <- rbind(covar_df,covar_df2)

rm(covar_df2,temp,dup_test,near3_miss,near3_raw,nearest3,soilcats,soilcats2,
   t,temp2,temp3,alldets,duplicated,hasdata,i,missing,ndx, soil_varnames,
   thesecols,thisvar,uncertain,v)

head(master_df)
nrow(master_df)

# Weather covariates (temporal covars) -----------------

dat$temporalcovars <- read.csv("TemporalCovars2.csv",header=T)
summary(dat$temporalcovars)
head(dat$temporalcovars)
names(dat$temporalcovars)

temp <- strsplit(dat$temporalcovars$Sess.Gd.Occ.ID,split="_")
dat$temporalcovars$Session2 <- as.numeric(sapply(temp,function(t) t[1]))
dat$temporalcovars$Grid <- as.numeric(sapply(temp,function(t) t[2]))
dat$temporalcovars$Occasion <- as.numeric(sapply(temp,function(t) t[3]))

alloccasions <- sort(unique(dat$temporalcovars$Occasion))
allgrids <- sort(unique(dat$temporalcovars$Grid))
allsessions <- sort(unique(dat$temporalcovars$Session2))

mintemp <- list()
maxtemp <- list()
precip <- list()
wind <- list()
moon <- list()

g=1
for(g in 1:length(allgrids)){
  thisgrid <- allgrids[g]
  thisgrid2 <- sprintf("Grid_%s",thisgrid)
  temp <- subset(dat$temporalcovars,Grid==thisgrid)
  thissess <- sort(unique(temp$Session2))
  wind[[thisgrid2]] <- matrix(0,nrow=length(thissess),ncol=length(alloccasions))
  mintemp[[thisgrid2]] <- matrix(0,nrow=length(thissess),ncol=length(alloccasions))
  maxtemp[[thisgrid2]] <- matrix(0,nrow=length(thissess),ncol=length(alloccasions)) 
  precip[[thisgrid2]] <- matrix(0,nrow=length(thissess),ncol=length(alloccasions))
  moon[[thisgrid2]] <- matrix(0,nrow=length(thissess),ncol=length(alloccasions))
  s = 1
  for(s in 1:length(thissess)){
    thissess2 <- thissess[s]
    temp2 <- subset(temp,Session==thissess2)
    temp2 <- temp2[order(temp2$Occasion),]
    wind[[thisgrid2]][s,] <- temp2$Avg.Wind[1:length(alloccasions)]
    mintemp[[thisgrid2]][s,] <- temp2$Min.Temp[1:length(alloccasions)]
    maxtemp[[thisgrid2]][s,] <- temp2$Max.Temp[1:length(alloccasions)]
    precip[[thisgrid2]][s,]  <- ifelse(temp2$Precip..mm.[1:length(alloccasions)]>0,1,0)
    moon[[thisgrid2]][s,] <- temp2$Max.Temp[1:length(alloccasions)]
  }
  
}

mintemp$Grid_1

rm(temp, temp2, g, s, thisgrid, thisgrid2,thissess, thissess2)

# PKM captures -----------------

dat$PKMcaps <- read.csv("PKM_caps2.csv",header=T)
summary(dat$PKMcaps)

head(dat$PKMcaps)

dets <- as.numeric(row.names(dat$trapcovars))
# cbind(dets,master_df$Det)  # ok

t=dets[1]
master_df$totcaps <- sapply(dets,function(t) nrow(subset(dat$PKMcaps,Det==t)))

# # Other species captures ------------------ 
 
dat$othercaps <- read.csv("Other_caps2.csv",header=T)
names(dat$othercaps)[names(dat$othercaps)=="Gd.St"] <- "Det"
#dat$othercaps$Gd.St

dat$othercaps[is.na(dat$othercaps)] <- 0

alldets <- sort(unique(master_df$Det))
alldets2 <- sort(unique(dat$othercaps$Det))
alldets3 <- sort(unique(dat$PKMcaps$Det))

all(alldets2%in%alldets)    # okay
all(alldets%in%alldets2)    # okay
all(alldets3%in%alldets)    # okay
all(alldets%in%alldets3)    # okay - some traps never had PKM

summary(dat$othercaps)
head(dat$othercaps)

setdiff(master_df$Det,dat$othercaps$Det)   # check
setdiff(dat$othercaps$Det,master_df$Det)   # okay

setdiff(dat$PKMcaps$Det,dat$othercaps$Det)  # okay
nopkm <- setdiff(dat$othercaps$Det,dat$PKMcaps$Det)   # traps with no PKM captures
nopkm_but <- setdiff(nopkm,dat$othercaps$Det[dat$othercaps$MIPA==0])


#unique(dat$othercaps$X)   # not useful  

dat$othercaps$PKMCAPS2 <- master_df$totcaps[match(dat$othercaps$Det,master_df$Det)]
check <- dat$othercaps[,c("MIPA","PKMCAPS2")]
check

## debug- look for discrepancies and try to explain
discrep <- dat$othercaps$Det[which((dat$othercaps$MIPA==0&dat$othercaps$PKMCAPS2>0)|(dat$othercaps$MIPA>1&dat$othercaps$PKMCAPS2==0))]   
#discrep2 <- dat$othercaps$Det[which((dat$othercaps$M_MIPA>0)&(dat$othercaps$PKMCAPS2>0)&(dat$othercaps$M_MIPA!=dat$othercaps$PKMCAPS2))]   

nopkm_butpresent <- intersect(discrep,nopkm_but)
pkm_butabsent <- setdiff(discrep,nopkm_but)

trap <-discrep[1] 
subset(dat$PKMcaps,Det==trap )
dat$othercaps[dat$othercaps$Det==trap ,]

dat$othercaps[dat$othercaps$Det%in%discrep ,]

length(discrep)

names(dat$othercaps)
otherspec <- c("DIDE","DIME","PELO","DIOR","ONsp","PEMA","PETR")
#otherspec <- c("DIDE","DIME","DIOR","PELO")    # NOTE: separate desert K from the others
nonseed <- c("PEMA", "ONsp", "PETR")
seedeaters <- setdiff(otherspec,nonseed)
desertk <- c("DIDE")
krats <- c("DIDE","DIME","DIOR")
krats2 <- c("DIME","DIOR")

dets <- as.numeric(row.names(dat$trapcovars))

ndx <- match(master_df$Det,dat$othercaps$Det)
temp <- sapply(otherspec,function(t) master_df[[t]] <<- dat$othercaps[[t]][ndx]  )  

head(master_df)


master_df$KratCaps <- apply(master_df[,krats],1,sum)
#master_df$KratCaps2 <- apply(master_df[,krats2],1,sum)
master_df$seedeaters <- apply(master_df[,seedeaters],1,sum)

rm(check,temp,alldets,alldets2,alldets3,dets,discrep,ndx,nopkm,nopkm_but,nopkm_butpresent,pkm_butabsent,t,trap)

covar_df2 <- data.frame(
  origname = c(otherspec,"KratCaps","seedeaters"),
  std = sprintf("%s_std",c(otherspec,"KratCaps","seedeaters")),
  clearname = c("DesertKrat","MerriamKrat","LittlePocket","OrdsKrat","ON_sp","DeerMouse","PinyonMouse","K-rats","Seed eaters"),
  stringsAsFactors = F
)


# Standardize covariates

covar_df2$mean <- sapply(covar_df2$origname,function(t) mean(master_df[[t]],na.rm=T) )
covar_df2$min <- sapply(covar_df2$origname,function(t) min(master_df[[t]],na.rm=T) )
covar_df2$max <- sapply(covar_df2$origname,function(t) max(master_df[[t]],na.rm=T) )
covar_df2$sd <- sapply(covar_df2$origname,function(t) sd(master_df[[t]],na.rm=T) )

covar_df <- rbind(covar_df,covar_df2)

temp <- sapply(1:length(covar_df2$origname), function(t)   # standardize all covars
  master_df[[covar_df2$std[t]]] <<- (master_df[[covar_df2$origname[t]]]-covar_df2$mean[t])/covar_df2$sd[t] )


rm(covar_df2)


# Visualize covars  -----------------------

# covar_df <- covar_df[-which(covar_df$origname=='Total'),]

i=1
for(i in 1:nrow(covar_df)){
  hist(master_df[[covar_df$std[i]]],main=covar_df$clearname[i])
}

# severe outlier in soil-gravel

# max(master_df$Gvl_std)
# max(master_df$Gvl)
# master_df$Gvl_std[master_df$Gvl_std>4] <- 4


# Visualize response -------------------

hist(master_df$totcaps)


# Pick a sand variable ---------------

covar_df$origname
sandvars <- c("medfinevf","allusda","coarsemed","orig","medfine","medium","nottoofine","finevf","fine")
sandvars <- c(sandvars,"cobble","gravel","coarsesand","sand3","sand2","sand1","finesand","veryfine","siltclay")


## loop through all sand vars- which ones are most correlated with PKM?

sandtest = data.frame(var=sandvars)
sandtest$effsize <- NA 
sandtest$rsquared <- NA
sandtest$pval <- NA
i=13
for(i in 1:length(sandvars)){
  thisvar <- sandvars[i]
  thisvar_std <- paste0(thisvar,"_std")
  # master_df[[thisvar_std]]
  # master_df$totcaps
  form <- formula( sprintf("totcaps~%s",thisvar_std))
  temp100 <- subset(master_df,sand3_std<5.5)  # remove outliers...
  thismod <- glm.nb(form,data=temp100)  #glm.nb(form,data=master_df)
  nullmod <- glm.nb(totcaps~1,data=temp100)
  thissum <- summary(thismod)
  # plot(form,data=master_df)
  sandtest$effsize[i] <- coef(thismod)[2]
  sandtest$rsquared[i] <- MuMIn::r.squaredLR(thismod,nullmod)
  sandtest$pval[i] <- thissum$coefficients[2,4]
}

sandtest  # none perform very well!

sandvar <-  "sand3"   # "medfinevf" # "allusda"

covar_df_orig <- covar_df
covar_df <- covar_df[!covar_df$origname%in%setdiff(sandvars,sandvar),]

# cor(master_df[,sandvars])  # correlation of sand vars...


# Correlation analyses ----------------

temp <- cor(master_df[,covar_df$std],use="complete.obs")

ndx <- which(abs(temp)>0.7&temp<1,arr.ind = T)

temp[rownames(ndx),rownames(ndx)]

covar_df <- covar_df[covar_df$std!="Sh_Cnt_std",]   # remove shrub count- keep canopy pct
covar_df <- covar_df[covar_df$std!="Sh_SR_std",]    # remove shrubSR, not sure what this is anyway!

covar_df <- covar_df[covar_df$std!="Q_FbsAvg_std",]   # remove forbs AND grasses separately (model together)
covar_df <- covar_df[covar_df$std!="Q_GrsAvg_std",]   # a priori.. 
covar_df <- covar_df[covar_df$std!="Q_BGrsCnt_std",]   # a priori.

# covar_df <- covar_df[covar_df$std!="Cbl_std",]    
# covar_df <- covar_df[covar_df$std!="Gvl_std",]    # a priori reasoning  -- use "coarse" instead

# covar_df <- covar_df[covar_df$std!="Silt_std",]    
# covar_df <- covar_df[covar_df$std!="Clay_std",]   
# covar_df <- covar_df[covar_df$std!="Total_std",]   
covar_df



# delete quadrat variables related to soil: coarse gravel and gravel, cobble, rock and boulder. 

covar_df <- covar_df[covar_df$std!="Q_CGvlAvg_std",]   
covar_df <- covar_df[covar_df$std!="Q_CRBAvg_std",]   
covar_df <- covar_df[covar_df$std!="Q_GvlAvg_std",]   

# delete quadrat shrub average!

covar_df <- covar_df[covar_df$std!="Q_ShbAvg_std",] 

# delete small and mini dunes [Q: did we remove these due to high correlation??]

covar_df <- covar_df[covar_df$std!="Q_SDnCnt_std",]   
covar_df <- covar_df[covar_df$std!="Q_MDnCnt_std",]   

covar_df


# finally, delete "Fine", due to zero-sum constraint
# covar_df <- covar_df[covar_df$std!="Fine_std",] 
# covar_df

nrow(master_df)

# finally, delete other species except for total krats
#covar_df <- covar_df[covar_df$std!="DIDE_std",] 
covar_df <- covar_df[covar_df$std!="DIOR_std",] 
covar_df <- covar_df[covar_df$std!="DIME_std",] 
covar_df <- covar_df[covar_df$std!="PELO_std",]
covar_df <- covar_df[covar_df$std!="DIDE_std",]
covar_df <- covar_df[covar_df$std!="ONsp_std",]
covar_df <- covar_df[covar_df$std!="PEMA_std",]
covar_df <- covar_df[covar_df$std!="PETR_std",]
#covar_df <- covar_df[covar_df$std!="KratCaps_std",]
covar_df <- covar_df[covar_df$std!="seedeaters_std",]
covar_df

temp <- cor(master_df[,covar_df$std],use="complete.obs")
temp

covar_df <- covar_df[covar_df$std!="Sh_NNAvg_std",]   # high correlation with two others
# covar_df <- covar_df[covar_df$std!="Crs_std",]   # high negative correlation with sand

# hist(master_df$Sh_SR)

temp <- cor(master_df[,covar_df$std],use="complete.obs")
temp
ndx <- which(abs(temp)>0.7&temp<1,arr.ind = T)
temp[rownames(ndx),rownames(ndx)]    

covar_df$clearname[covar_df$origname==sandvar] <- "Coarse-Medium Sand"
covar_df

# Visualize response vs covars at local trap scale ----------------------

par(mfrow=c(1,1))
for(i in 1:nrow(covar_df)){
  plot(master_df[[covar_df$std[i]]],master_df$totcaps,main=covar_df$clearname[i])    # not much really comes out here
}


# Perform initial diagnoses of appropriate count distributions ----------------------

# test Poisson
fit <- vcd::goodfit(master_df$totcaps,type="poisson")               # test fails for Poisson count model every time
summary(fit)
vcd::rootogram(fit)                                                 # bad fit
vcd::Ord_plot(master_df$totcaps)                                    # this diagnoses that neg binom distribution is appropriate!

# test NegBin
fit <- vcd::goodfit(master_df$totcaps, type="nbinom",method="ML")   # goodness of fit test for neg binom
summary(fit)                                                  
vcd::rootogram(fit)                                           
vcd::distplot(master_df$totcaps, type="nbinom")               # Neg binom fits great!


# determine which distribution to use ---------------- 

covar_df

form <- as.formula(sprintf("totcaps ~ %s + (1|Gd)",paste(covar_df$std,collapse="+")))

form2 <- as.formula(sprintf("totcaps ~ (%s)^2 + (1|Gd)",paste(covar_df$std,collapse="+")))

### try poisson regression with no zero component

test1.mod <- glmmTMB(form,
                     master_df,
                     ziformula = ~0 ,
                     family= poisson(link = "log"))

summary(test1.mod)


test1.res <- DHARMa::simulateResiduals(test1.mod,n=300)
plot(test1.res)                        # okay fit
DHARMa::testUniformity(test1.res)      # fail
DHARMa::testResiduals(test1.res)       # fails outlier test, dispersion test, 
DHARMa::testOutliers(test1.res,type="bootstrap")  # fails outlier test

DHARMa::testZeroInflation(test1.res)   # barely passes


### try negative binomial regression with no zero component

test2.mod <- glmmTMB(form,
                     master_df,
                     ziformula = ~0 ,
                     dispformula = ~1,            
                     family= nbinom1(link = "log"))

summary(test2.mod)

test2.res <- DHARMa::simulateResiduals(test2.mod,n=300)
plot(test2.res)      
DHARMa::testUniformity(test2.res)         # pass
DHARMa::testResiduals(test2.res)          # passes all tests
DHARMa::testOutliers(test1.res,type="bootstrap")  # passes outlier test with bootstrap
DHARMa::testZeroInflation(test2.res)      # passes zero inflation test


# AIC model selection table
bbmle::AICtab(test1.mod,test2.mod,weights=TRUE,mnames=c("Poisson","NegBin") )   # negbin wins!


# Try grid-level analysis -------------------------

# TODO: use weighted regression

meandensities2 <- read.csv("meandens_jags.csv",row.names = 1)

meandensities2$Grid <- as.numeric(gsub("grid_","",rownames(meandensities2)))
#master_df$Gd

meandensities2 <- meandensities2[order(meandensities2$Grid),]

ngrids <- length(unique(master_df$Gd))
allgrids <- unique(master_df$Gd)
gd_df <- data.frame(temp=numeric(ngrids))

temp <- sapply(covar_df$std,function(t) gd_df[[t]] <<- tapply(master_df[[t]],master_df$Gd,mean,na.rm=T ) )

gd_df[["meancaps"]] <- tapply(master_df$totcaps,master_df$Gd,mean )
gd_df <- gd_df[,-1]
gd_df

rownames(gd_df) <- sort(allgrids)

gd_df$meddens <- meandensities2$X50.[match(rownames(gd_df),meandensities2$Grid)]
gd_df$meddensu <- meandensities2$X97.5.[match(rownames(gd_df),meandensities2$Grid)]
gd_df$meddensl <- meandensities2$X2.5.[match(rownames(gd_df),meandensities2$Grid)]

#response <- "meancaps"
response <- "meddens"

ndx <- 1:nrow(covar_df)
temp <- cor(gd_df[,covar_df$std[ndx]])
ndx <- which(abs(temp)>0.7&temp<1,arr.ind = T)
temp[rownames(ndx),rownames(ndx)]    

covar_df$std
# ndx <- c(2,3,4,5,6,7,8,9)
# covar_df$std[ndx]

gd_df$weights <-  1/(gd_df$meddensu-gd_df$meddens)^2

form <- as.formula(sprintf("%s ~ %s",response,paste(covar_df$std,collapse="+")))
grd_mod1 <- lm(form,data=gd_df,weights=weights)
summary(grd_mod1)

grd_mod2 <- stepAIC(grd_mod1,k=log(nrow(gd_df)))
summary(grd_mod2)      # not bad!    Rsquared of .5
#plot(grd_mod2)

temp <- as.character(grd_mod2$call$formula)[3]
ndx2 <- which(covar_df$std%in% strsplit(temp,split = " \\+ ")[[1]])

form2 <- as.formula(sprintf("%s ~ (%s)^2",response,temp))
grd_mod_int <- lm(form2,data=gd_df,weights=weights)
grd_mod_int2 <- stepAIC(grd_mod_int,k=log(nrow(gd_df)))
summary(grd_mod_int2)     # interaction between sand and seed eaters?


 # visualize the relationships (grid level)
graphics.off()

# svg("sitelevelfig4.svg",5.5,7)
png("sitelevelfig4.png",5.5,7,units="in",res=600)
par(mfrow=c(3,1))
par(mai=c(0.8,0.9,0.2,0.1))

ordr <- order(meandensities2$X50.,decreasing = T)
Hmisc::errbar(1:nrow(meandensities2),meandensities2$X50.[ordr],meandensities2$X97.5.[ordr],meandensities2$X2.5.[ordr], ylim=c(0,20),ylab="PKM Density (per ha)",xlab="",xaxt="n")
axis(1,at=1:nrow(meandensities2),labels = rownames(meandensities2[ordr,]),las=2)
#title("Mean PKM Densities")

fig_label("a)",region="figure","topleft",cex=2)
#legend("topleft",legend="a)",cex=1.5,bty="n")

covar_df$clearname[covar_df$origname=="Sh_CpyPct"] <- "Shrub Canopy (%)"
covar_df$clearname[covar_df$origname=="Q_BGrdAvg"] <- "Bare Ground (%)"

counter=1
i=covar_df$std[ndx2][2]
for(i in covar_df$std[ndx2[1:2]]){
  # plot(gd_df[[i]],gd_df[[response]],main=covar_df$clearname[covar_df$std==i],
  #      xlab=covar_df$clearname[covar_df$std==i],ylab=response) 
  Hmisc::errbar(gd_df[[i]],gd_df[[response]],gd_df$meddensu,gd_df$meddensl,main="", # covar_df$clearname[covar_df$std==i]
                xlab=covar_df$clearname[covar_df$std==i],ylab="PKM density (per ha)",xaxt="n")
  newdata <- data.frame(
    this = seq(min(gd_df[[i]]),max(gd_df[[i]]),0.01)
  )
  others <- setdiff(covar_df$std[ndx2],i)
  for(j in others){
    newdata[[j]] = mean(gd_df[[j]]) 
  }
  names(newdata)[1] <- i
  ys <- predict(grd_mod_int2,newdata=newdata,interval="confidence")
  lines(newdata[[i]],ys[,1],lwd=2)
  lines(newdata[[i]],ys[,2],lty=2)
  lines(newdata[[i]],ys[,3],lty=2)
  fig_label(c("b)","c)")[counter],region="figure","topleft",cex=2)
  axis(1,at=seq(-2,2,0.5),labels = round(covar_df$mean[covar_df$std==i]+covar_df$sd[covar_df$std==i]*seq(-2,2,0.5)) )
  counter=counter+1
}

dev.off()

temp <- cor(gd_df[,covar_df$std[ndx2]])

temp

# NOTE: old model had sand as a top variable at the grid level. This effect was not recovered after correcting the sand categories...

# Neighborhood covariate values #1: use raw observed data --------------------

plot(dat$traplocs$geometry[dat$traplocs$Gd==18])  # look at grids individually- not square!

trapdists <- st_distance(dat$traplocs,dat$traplocs)
trapdists <- as.matrix(trapdists)
class(trapdists) <- "numeric"
colnames(trapdists) <- dat$traplocs$Det
rownames(trapdists) <- dat$traplocs$Det
# trapdists[1,]

mean(apply(trapdists,1,function(t) min(t[which(t>0)])))    # mean dist between traps is 8 m


# nrow(subset(master_df,Gd==1))   # number of traps per grid

#nns <- c(1,9,25,64)
maxdists <- c(1,8,16,24,90)
maxdists2 <- sapply(maxdists,function(t) sqrt(t^2*2))
maxdists3 <- round(maxdists2)


nn <- maxdists2[1]
for(nn in maxdists2){
  covar_df[[sprintf("std_%s",round(nn))]] <- sapply(covar_df$std,function(t) sprintf("%s_%snn",t,round(nn)))
  covar_df[[sprintf("orig_%s",round(nn))]] <- sapply(covar_df$origname,function(t) sprintf("%s_%snn",t,round(nn)))
  temp <- sapply(covar_df[[sprintf("std_%s",round(nn))]],function(t) master_df[[t]] <<- numeric(nrow(master_df)))
  temp <- sapply(covar_df[[sprintf("orig_%s",round(nn))]],function(t) master_df[[t]] <<- numeric(nrow(master_df)))
  i=1
  for(i in 1:length(master_df$Det)){
    thisdet <- as.character(master_df$Det[i])
    distsfrom <- sort(trapdists[,thisdet],decreasing = F)
    nearest <- as.numeric(names(distsfrom[which(distsfrom<=nn)]))
    
    tmp2 <-  subset(master_df,Det%in%nearest)
    j=1
    for(j in 1:nrow(covar_df)){
      master_df[[covar_df[[sprintf("std_%s",round(nn))]][j]]][i] <- mean(tmp2[[covar_df$std[j]]],na.rm=T)
      master_df[[covar_df[[sprintf("orig_%s",round(nn))]][j]]][i] <- mean(tmp2[[covar_df$origname[j]]],na.rm=T)
    }
  }
}

summary(master_df)
hist(master_df$Q_FbGsAvg_std_127nn)


# Find the best scale for each variable (observed env data) --------------------------------

# 
# graphics.off()
# 

master_df2 <- master_df[complete.cases(master_df[,!names(master_df)=="Total_std"]),]

master_df2$Gd2 <- as.numeric(as.factor(master_df2$Gd))
master_df2$Gd3 <- as.factor(master_df2$Gd2)
# 
# 
# var=1
# bestscales <- numeric(nrow(covar_df))
# for(var in 1:nrow(covar_df)){
#   nn=maxdists3[2]
#   counter=1
#   rsqs <- numeric(length(maxdists3))
#   for(nn in maxdists3){
#     form <- as.formula(sprintf("totcaps ~ %s + (1|Gd2)",covar_df[[sprintf("std_%s",nn)]][var]))
#     temp1 <- glmmTMB(form,
#             master_df2,
#             ziformula = ~0 ,
#             dispformula = ~1,                     
#             family= nbinom1(link = "log"))
#     rsqs[counter] <- MuMIn::r.squaredGLMM(temp1)["trigamma",2]
#     counter=counter+1
#   }
#   bestscales[var] <- maxdists3[which.max(rsqs)]
#   plot(maxdists3,rsqs,main=covar_df$clearname[var],type="b")
# }
# 
# names(bestscales) <- covar_df$clearname
# bestscales

#?glmmLasso::glmmLasso (no nb dist)

#?mpath::glmregNB (no mixed model)


# Find the best scale for each variable using random forest?

# note: using importance probably doesn't make sense, since correlations with other vars may cause the importance to go up!
#    solution: use R-squared instead. this is the scale that explains the most variation, considering interactions etc..   


# Try predicting using OOB=T?  

#library(party)
#library(randomForest)
library(ranger)

graphics.off()

## run initial random forest model

thesecovs <- sprintf("%s_%snn",covar_df$std,maxdists3[3])
thisform <-  as.formula(sprintf("totcaps ~ %s",paste(thesecovs,collapse = "+") ))
thismod <- ranger::ranger(thisform,data=master_df2,importance = "permutation")
varimp <- sort(importance(thismod),decreasing = T)
varimp[1]

var=1
#bestscales2 <- numeric(nrow(covar_df))
bestscales <- numeric(nrow(covar_df))
names(bestscales) <- covar_df$origname

#rsqs <- matrix(NA,nrow=nrow(covar_df),ncol= length(maxdists3)+2)

prev <- NULL
svg("OptScales_all4.svg",8,8)
# png("OptScales_all4.png",8,8,units="in",res=600)
par(mfrow=c(3,3))
var=3
for(var in 1:nrow(covar_df)){    # loop through variables
  counter=1
  
  thisvar <- paste(strsplit(names(varimp)[var],"_")[[1]][-length(strsplit(names(varimp)[var],"_")[[1]])],collapse="_")
  thisvarndx <- which(covar_df$std==thisvar)
  
  rsqs <- numeric(length(maxdists3))
  
  #vars_less1 <- setdiff(names(varimp), thisvar)
  
  nn=maxdists3[1]
  for(nn in maxdists3){
    
    thisform <- as.formula(sprintf("totcaps ~ %s", paste(c(prev,sprintf("%s_%snn",thisvar,nn)),collapse = "+") ))
    

    thismod <- ranger::ranger(thisform,data=master_df2)
                            
    pred <- predict(thismod,data=master_df2)
    rsqs[counter] <- 1-sum((master_df2$totcaps - pred$predictions)^2)/sum((master_df2$totcaps - mean(master_df2$totcaps))^2)
    
    counter=counter+1
  }
  
  bestscales[thisvarndx] <- maxdists3[which.max(rsqs)]
  
  prev <- c(prev,sprintf("%s_%snn",thisvar,bestscales[thisvarndx]))
  
  plot(c(maxdists3),rsqs,main=covar_df$clearname[thisvarndx],type="b", xlab=covar_df$clearname[thisvarndx],
       ylab="R-squared")
  abline(v=bestscales[thisvarndx],lwd=2,lty=2)
  
}
dev.off()

## note: for final 'scales' figure, use inkscape to grab the plots needed for final figure. 

# randomForest::varImpPlot(thismod)
# partialPlot(thismod,x.var="Crs_std_11nn",pred.data = master_df2)

# colnames(rsqs) <- c(0,0,maxdists3)
# rownames(rsqs) <- covar_df$clearname
#  #  rownames(rsqs) <- c(covar_df$clearname[1:4],"ShrubHeight",covar_df$clearname[5:9])
# rsqs
# 
# bestscales3 <- apply(rsqs[,-c(1,2)],1,function(t) as.numeric(names(t)[which.max(t)]) )

#bestscales3["SoilSand"] <- 23

# names(bestscales3) <- covar_df$clearname
# bestscales2

# oldbest <- bestscales
# oldbest <- rbind(oldbest,bestscales3)
# 
# bestscales <- bestscales3
# bestscales["SoilSand"] <- oldbest["oldbest","SoilSand"]

# ### try RF cross validation
# trycv <- rfcv(trainx=master_df2[,paste(sprintf("%s_%snn",
#                           rep(covar_df$std,each=length(maxdists3)),   #
#                           rep(maxdists3,times=(nrow(covar_df)))))],
#      trainy=master_df2$totcaps, cv.fold=5, scale="log", step=0.25,
#      recursive=T)
# with(trycv, plot(n.var, error.cv, log="x", type="o", lwd=2))


# note: remove shrub height- models with it were not better than models without it

# covar_df2 <- covar_df[-which(covar_df$clearname=="Dist_NearestShrub"),]
# bestscales3 <- bestscales3[-5]


# Try Boruta?

# library(Boruta)
# vars_all <- paste(sprintf("%s_%snn",
#                             rep(covar_df$std,each=length(maxdists3)),   #
#                             rep(maxdists3,times=(nrow(covar_df)))))
# 
# form_all <- as.formula(sprintf("totcaps ~ (%s)",vars_all ))
# 
# result <- Boruta::Boruta(form_all,data=master_df2)
# 
# # trycv <- rfcv(trainx=master_df2[,paste(sprintf("%s_%snn",
# #                                                rep(covar_df$std,each=length(maxdists3)),   #
# #                                                rep(maxdists3,times=(nrow(covar_df)))))],
# #               trainy=master_df2$totcaps, cv.fold=5, scale="log", step=0.25,
# #               recursive=T)


# Try Caret recursive feature elimination -------------------

varstotry <- sprintf("%s_%snn",covar_df$std,bestscales)

master_df2$random1 <- runif(nrow(master_df2))
master_df2$random2 <- runif(nrow(master_df2))

# train a preliminary random forest algorithm using caret/ranger (takes a while- thinned dataset to allow faster computation since this is a preliminary step )

varstotry <- c(varstotry,"random1","random2")
response <- "totcaps"

formula = as.formula(paste0(response, "~", paste(varstotry,collapse = "+")))

# folds <- list()
# for(i in 1:5){
folds <- groupKFold(master_df2$Gd,k=5)   # changed to leave out individual grids
# }

control <- caret::trainControl(method="cv", index=folds, allowParallel = TRUE)   # changed from repeatedcv  number=10, 

model <- caret::train(formula, data=master_df2, method = 'ranger',
                      trControl=control, importance = 'permutation', preProcess="scale")   # method="ranger" ,


# estimate variable importance using ranger/caret

importance <- caret::varImp(model, scale=FALSE)


# summarize importance

print(importance)

# plot variable importance

plot(importance)   

## KTS: we can remove any variables less important than random. Let's also remove curvature, precip_warm, etc cuz they are near random

randimp <- max(importance$importance[c("random1","random2"),])*1.5
covars_new <- rownames(importance$importance)[importance$importance>randimp]    # change to whatever the importance of the random var is...

covars_new <- c(covars_new,c("random1","random2"))
formula <- as.formula(paste0(response,"~",paste(covars_new,collapse = "+")))   # formula with (potentially) reduced covars


# define the settings for recursive feature elimination [kts- might want to change number to 10, below]
folds <- groupKFold(master_df2$Gd,k=10)   # changed to leave out individual grids
control <- caret::rfeControl(functions=rfFuncs, method="cv", index = folds, rerank=F)   # number=5,

# run the RFE algorithm (takes a while to run)

results <- caret::rfe(master_df2[,covars_new], master_df2[,response], sizes=c(1:length(covars_new)), rfeControl=control)

# remember to save the workspace after running this!

# summarize the results

print(results)

# list the chosen features

predictors(results)      

# plot the results

png("rfefig4.png",4,3,units="in",res=600)
plot(results, type=c("g", "o"))    # visualize number of important features
dev.off()

bestvars <- results$optVariables[1:5]
bestvars

covars_new2 <- bestvars   #results$optVariables[1:7]

formula <- as.formula(paste0(response,"~",paste(covars_new2,collapse = "+")))   # formula with further reduced covars



# # Find the best scale for each variable (classified NAIP imagery)
# 
# var=1
# r_bestscales <- numeric(nrow(covar_df2))
# for(var in 1:nrow(covar_df2)){
#   nn=maxdists3[1]
#   counter=1
#   r_rsqs <- numeric(length(maxdists3))
#   for(nn in maxdists3){
#     form <- as.formula(sprintf("totcaps ~ %s + (1|Gd2)",covar_df2[[sprintf("orig_%s",nn)]][var]))
#     temp1 <- glmmTMB(form,
#                      master_df2,
#                      ziformula = ~0 ,
#                      dispformula = ~1,                     
#                      family= nbinom1(link = "log"))
#     r_rsqs[counter] <- MuMIn::r.squaredGLMM(temp1)["trigamma",2]
#     counter=counter+1
#   }
#   r_bestscales[var] <- maxdists3[which.max(r_rsqs)]
#   plot(maxdists3,r_rsqs,main=covar_df2$orig[var],type="b")
# }
# 
# names(r_bestscales) <- covar_df2$orig
# r_bestscales


# 
# # Correlation (field measured vars)
#
# 
# ndx <- c(1,2,3,4,5,6,7,8,9,10)
# temp <- cor(master_df2[,sprintf("%s_%snn",covar_df$std,bestscales3)][ndx])    # note: maybe try not eliminating shrub height?
# temp
# 
# ndx <- c(1,2,3,4,7,8,9,10)
# temp <- cor(master_df2[,sprintf("%s_%snn",covar_df$std,bestscales3)][ndx])
# temp
# 
# covar_df2 <- covar_df[ndx,]     # reduce covar_df to uncorrelated vars.



# Correlation analyses --------------------------

bestscales
temp <- cor(master_df[,paste0(covar_df$std,"_",bestscales,"nn")],use="complete.obs")

ndx <- which(abs(temp)>0.7&temp<1,arr.ind = T)

temp[rownames(ndx),rownames(ndx)]   # shrub height and canopy are correlated

bestvars <- setdiff(bestvars,"Sh_HtAvg_std_23nn")  # remove shrub height if included in best vars



# Determine non-linearities and interactions using the best variable set! ------------------

# covar_df
# bestvars <- paste0(covar_df$std,"_",bestscales,"nn")[c(1,3,8,9)]
# covars_new2 <- bestvars   #results$optVariables[1:7]
# response <- "totcaps"
# formula <- as.formula(paste0(response,"~",paste(covars_new2,collapse = "+")))   # formula with further reduced covars


temp <- strsplit(bestvars,"_")
ndx <- sapply(temp,function(t) which(t=="std"))
temp2 <- sapply(1:length(temp),function(t) temp[[t]][ndx[t]-1])
temp3 <- sapply(1:length(temp),function(t) if(temp2[t]=="CpyPct") temp2[t] <<- paste0("Sh_",temp2[t])  )
temp3 <- sapply(1:length(temp),function(t) if(temp2[t]=="SR") temp2[t] <<- paste0("Sh_",temp2[t])  )
temp3 <- sapply(1:length(ndx),function(t) if(ndx[t]>2&!temp2[t]%in%c("Sh_CpyPct","Sh_SR")) temp2[t] <<- paste0("Q_",temp2[t])  )
covar_df2 <- covar_df[covar_df$origname%in%temp2,]
covar_df2$best <- bestvars[match(covar_df2$origname,temp2)]


# Non-linearities ----------------------

library(ranger)
thismod <- ranger(formula,data=master_df2,importance = 'permutation')

varimp <- importance(thismod)
#plot(varimp)
#?importance.ranger

options(scipen=20)

varimp <- varimp[order(varimp,decreasing=F)]

png("varimp_finalvars4.png",4,4,units="in",res=600)
par(mai=c(1,2,0.1,0.1))
par(las=2)
plot(1:length(varimp)~varimp,yaxt="n",xlab="importance",ylab="",xlim=c(0,1))
axis(2,at=1:length(varimp),labels = covar_df2$clearname[match(names(varimp),covar_df2$best)]  )
dev.off()

# library(abcrf)
# variableImpPlot(varimp)

## partial dependence plots

png("effects_topvars4_rf.png",6,4,units="in",res=600)
par(mai=c(.8,.8,.4,.1))
par(mfrow=c(2,3))
p=1
for(p in 1:length(bestvars)){
  thisvar <- bestvars[p]
  
  nd <- data.frame(x=seq(min(master_df2[[thisvar]]),max(master_df2[[thisvar]]),length=50))
  names(nd) <- thisvar
  
  othervars <- setdiff(bestvars,thisvar)
  temp <- sapply(othervars,function(t) nd[[t]] <<- mean(master_df2[[t]])  )
  #nd
  
  pred = predict(thismod,data=nd,type="response")$predictions
  
  plot(pred~nd[,1],type="l",xlab=thisvar,main=covar_df2$clearname[covar_df2$best==thisvar])
  rug(master_df2[[thisvar]])
  
}
dev.off()


# Interactions ------------------------


allcomb <- as.data.frame(t(combn(bestvars,2)))
names(allcomb) <- c("var1","var2")

allcomb$int1 <- NA
allcomb$int2 <- NA

p=1
for(p in 1:nrow(allcomb)){
  var1 = allcomb$var1[p]
  var2 = allcomb$var2[p]
  
  all1 = seq(min(master_df2[[var1]]),max(master_df2[[var1]]),length=10)
  all2 = seq(min(master_df2[[var2]]),max(master_df2[[var2]]),length=10)
  
  nd <- expand.grid(all1,all2)
  names(nd) <- c(var1,var2)
  
  othervars <- setdiff(bestvars,c(var1,var2))
  temp <- sapply(othervars,function(t) nd[[t]] <<- mean(master_df2[[t]])  )
  
  pred = predict(thismod,data=nd,type="response")$predictions
  
  additive_model <- lm(pred~as.factor(nd[[var1]])+as.factor(nd[[var2]]))
  
  pred_add = predict(additive_model)
  
  allcomb$int1[p] <- sqrt(mean((pred-pred_add)^2))
  
  maximp <- mean(varimp[c(var1,var2)])
  
  allcomb$int2[p] <- allcomb$int1[p]/maximp
  
}

allcomb[order(allcomb$int1,decreasing = T),]
allcomb <- allcomb[order(allcomb$int2,decreasing = T),]
allcomb   

### visualize interaction

ints.torun <- 1:4
int=3
for(int in 1:length(ints.torun)){
  png(sprintf("rf_int4%s.png",int),5,5,units = "in",res = 600)
  par(mfrow=c(1,1))
  thisint <- ints.torun[int]
  var1 = allcomb$var2[thisint]
  var2 = allcomb$var1[thisint]
  
  all1 = seq(min(master_df2[[var1]]),max(master_df2[[var1]]),length=20)
  all2 = seq(min(master_df2[[var2]]),max(master_df2[[var2]]),length=20)
  
  nd <- expand.grid(all1,all2)
  names(nd) <- c(var1,var2)
  
  othervars <- setdiff(bestvars,c(var1,var2))
  temp <- sapply(othervars,function(t) nd[[t]] <<- mean(master_df2[[t]])  )
  
  pred = predict(thismod,data=nd,type="response")$predictions
  
  predmat = matrix(pred,nrow=length(all1),ncol=length(all2))
  
  persp(all1,all2,predmat,theta=45,phi=35,xlab=var1,ylab=var2,zlab="sel int")
  dev.off()
}


### use top two interactions?


# # Correlation (NAIP derived vars)
# 
# ndx <- c(1,2,3,4,5,6,7,8)
# temp <- cor(master_df[,covar_df2[[sprintf("orig_%s",maxdists3[3])]]][ndx])
# temp
# 
# ndx <- c(3,4,8)
# temp <- cor(master_df[,covar_df2[[sprintf("orig_%s",maxdists3[3])]]][ndx])
# temp
# 
# covar_df2 <- covar_df2[ndx,]     # reduce covar_df to uncorrelated vars.
# 

# Build a multiple GLMM with vars at their best scales... (observed vars) -----------------


#### First try using a lasso to find the best variable set

# form <- as.formula(sprintf("totcaps ~ %s + 0",paste(sprintf("%s_%snn",covar_df$std,bestscales2),collapse="+") ))
# 
# form2 <- as.formula(sprintf("totcaps ~ (%s)^2 + 0",paste(sprintf("%s_%snn",covar_df$std,bestscales2),collapse="+") ))
# 
# x <- model.matrix(form2,master_df2)
# 
# #temp <- glmnet(x,master_df$totcaps,family="poisson")
# temp <- cv.glmnet(x,master_df2$totcaps,family="poisson")   # foldid = master_df2$Gd2
# 
# #plot(temp)
# #opt = temp$lambda.min
# coefs <- coef(temp, s = "lambda.min")
# 
# bestcoefs <- rownames(coefs)[abs(coefs[,1])>0.000000001][-1]
# 
# bestform <- as.formula(sprintf("totcaps ~ %s + (1|Gd2)",paste(bestcoefs,collapse="+")))
# 
# bestmod_alt <-  glmmTMB(bestform,
#                  master_df2,
#                  ziformula = ~0 ,
#                  dispformula = ~1,
#                  family= nbinom1(link = "log"))
# 
# summary(bestmod_alt)    # okay, not great!
# 
# rsq <- MuMIn::r.squaredGLMM(bestmod_alt)
# rsq


# ####   Use backward stepwise selection to find best model (takes a while!)
# 
# form <- as.formula(sprintf("totcaps ~ %s + (1|Gd2)",paste(sprintf("%s_%snn",covar_df2$std,bestscales3[ndx]),collapse="+") ))
# 
# form2 <- as.formula(sprintf("totcaps ~ (%s)^2 + (1|Gd2)",paste(sprintf("%s_%snn",covar_df2$std,bestscales3[ndx]),collapse="+") ))
# 
# #bestmod <- buildmer::buildglmmTMB(form2, ziformula = ~0,data=master_df2, dispformula = ~1, direction="backward",crit="AIC",
 #                                                      reduce.random=F,family=nbinom1(link = "log"))

# Alternative: use the top interactions from RF (conditional forest) (this is better I think!)
   # todo: try quadratic terms on dist nearest shrub and sand? or not...

# mainterms <- paste(covar_df2$best,collapse="+")
# 
# ints <- paste(sapply(1:2,function(t) paste(c(allcomb$var1[t],allcomb$var2[t]),collapse=':') ),collapse="+")
# 
# polyterms <- sprintf("poly(%s,2)+poly(%s,2)+poly(%s,2)+poly(%s,2)",bestvars[1],bestvars[2],bestvars[3],bestvars[4])
# 
# polyterms_int <- sprintf("(poly(%s,2)+poly(%s,2)+poly(%s,2)+poly(%s,2))^2",bestvars[1],bestvars[2],bestvars[3],bestvars[4])
# 
# mainterms_int <- sprintf("(%s+%s+%s+%s)^2",bestvars[1],bestvars[2],bestvars[3],bestvars[4])
# 
# covar_df2$best
# 
# form <- as.formula(sprintf("totcaps ~ %s + (1|Gd2)",mainterms ))
# 
# form2 <- as.formula(sprintf("totcaps ~ %s + %s + (1|Gd2)",mainterms,ints ))
# 
# form2_2 <- as.formula(sprintf("totcaps ~ %s + %s",mainterms,ints ))
# 
# form3 <- as.formula(sprintf("totcaps ~ %s + (1|Gd2)",polyterms ))
# 
# form4 <- as.formula(sprintf("totcaps ~ %s + (1|Gd2)",polyterms_int ))
# 
# form5 <- as.formula(sprintf("totcaps ~ %s + (1|Gd2)",mainterms_int ))
# 
# 
# bestmod <- buildmer::buildglmmTMB(form3, ziformula = ~0,data=master_df2, dispformula = ~1, direction="backward",#crit="AIC",
#                         family=nbinom1(link = "log"),REML=F)
# bestmod3 <- bestmod@model
# bestform3 <- as.formula(bestmod@model$modelInfo$allForm$formula)
# 
# summary(bestmod3)  

# mod2 <- glmmTMB(form2,ziformula = ~0,data=master_df2, dispformula = ~1, family=nbinom1(link = "log"),REML = F)
# summary(mod2)
# 
# mod3 <- MASS::glm.nb(form2_2,data=master_df2)
# summary(bestmod2)

# mod_full <- glmmTMB(form2, ziformula = ~0,data=master_df2, dispformula = ~1, family=nbinom1(link = "log"),REML = F)
# summary(mod_full)
# AIC(mod_full)
# 
# mod_noint <- glmmTMB(form, ziformula = ~0,data=master_df2, dispformula = ~1, family=nbinom1(link = "log"),REML = F)
# summary(mod_noint)
# AIC(mod_noint,mod_full)

# TODO: use spline fit within model using splines::ns or something similar

bestvars
response
#library(splines)
# form <- as.formula("totcaps~Sh_CpyPct_std_34nn+(poly(Snd_std_23nn,2)+KratCaps_std_23nn+poly(Q_BGrdAvg_std_34nn),2)^2 + (1|Gd2)"  )
                   
# mod_full1 <- glmmTMB(form, ziformula = ~0,data=master_df2, dispformula = ~1, family=nbinom2(link = "log"),REML = F)
# summary(mod_full1)


form <- as.formula("totcaps~Q_BGrdAvg_std_34nn + 
                   # Q_FbGsAvg_std_23nn +
                   KratCaps_std_23nn + Sh_CpyPct_std_34nn + 
                   # sand3_std_23nn +   # sand not coming out as important
                   # Q_FbGsAvg_std_23nn:KratCaps_std_23nn + 
                   (1|Gd2)") 

mod_full2 <- glmmTMB(form, ziformula = ~0,data=master_df2, dispformula = ~1, family=nbinom2(link = "log"),REML = F)
summary(mod_full2)


# bestmod <- glmmTMB(bestform3,ziformula = ~0,data=master_df2, dispformula = ~1, family=nbinom1(link = "log"),REML = F)
# summary(bestmod)


form <- as.formula("totcaps~Sh_CpyPct_std_34nn + Q_FbGsAvg_std_23nn+
                   Q_BGrdAvg_std_34nn+KratCaps_std_23nn+
                   Q_FbGsAvg_std_23nn:Sh_CpyPct_std_34nn + (1|Gd2)")  
mod_full3 <- glmmTMB(form, ziformula = ~0,data=master_df2, dispformula = ~1, family=nbinom2(link = "log"),REML = F)
summary(mod_full3)


form <- as.formula("totcaps~Sh_CpyPct_std_34nn +
                   Q_BGrdAvg_std_34nn+KratCaps_std_23nn+ (1|Gd2)")  
mod_full4 <- glmmTMB(form, ziformula = ~0,data=master_df2, dispformula = ~1, family=nbinom2(link = "log"),REML = F)
summary(mod_full4)


bestmod <- mod_full4

nullmod <- glmmTMB(totcaps~(1|Gd2),ziformula = ~0,data=master_df2, dispformula = ~1, family=nbinom2(link = "log"),REML = F)
 
MuMIn::r.squaredGLMM(bestmod,nullmod)    # still not that great- R2 of 0.21

bbmle::AICtab(nullmod,bestmod,weights=TRUE,mnames=c("Null","Best") )   # best model wins!


tmp <- rownames(coef(summary(bestmod))$cond)[-1]

#topvars <- tmp[!grepl(":",tmp)]
topvars2 <- tmp #  bestvars        #c("Snd_std_23nn","Q_BGrdAvg_std_34nn")
# ints <- strsplit(tmp[grep(":",tmp)],":")

allvars <- c("totcaps","Gd2",topvars2)  #unique(c(topvars2,unlist(ints))))

i=2

## partial dependence plots --------------------

svg("PKM_main_effects5.svg",7.5,3)
# png("PKM_main_effects5.png",7.5,3,units="in",res=600)
layout(matrix(1:3,nrow=1,byrow = T))
par(mai=c(0.8,0.8,0.1,0.1))
if(length(topvars2)>0) tmp <- sapply(1:length(topvars2), function(i) VisualizeRelation(data=master_df2,
                                                                                     model=bestmod,
                                                                                     predvar=topvars2[i],
                                                                                     allvars = allvars) )     # run and save partial dependence plots for all top variables
dev.off()








#graphics.off()

svg("PKM_interactions.svg",8,8)
layout(matrix(1:1,nrow=1,byrow = T))
if(length(ints)>0 ) tmp <- sapply(1:length(ints), function(i) VisualizeInteraction(data=master_df2,
                                                                                   model=bestmod,
                                                                                   var1=ints[[i]][1],
                                                                                   var2=ints[[i]][2],
                                                                                   allvars=allvars) )
dev.off()


summary(bestmod3)








# OLD CODE ----------------------------





# ################
# # Cross validation?
# ################
# 
# library(caret)
# folds <- groupKFold(master_df2$Gd2)
# 
# pred_cv <- numeric(nrow(master_df2))
# pred_fit <- numeric(nrow(master_df2))
# 
# t = folds[[5]]
# temp <- lapply(folds,function(t){
#   thismd <- glmmTMB(bestform3, ziformula = ~0,data=master_df2[t,], dispformula = ~1,
#           family=nbinom1(link = "log"))
#   pred_cv[setdiff((1:nrow(master_df2)),t)] <<- predict(thismd,newdata=master_df2[setdiff((1:nrow(master_df2)),t),],type="response")
#   pred_fit[setdiff((1:nrow(master_df2)),t)] <<- predict(bestmod3,newdata=master_df2[setdiff((1:nrow(master_df2)),t),],type="response")
# })
# 
# 
# rsq_cv <- 1-sum((master_df2$totcaps - pred_cv)^2)/sum((master_df2$totcaps - mean(master_df2$totcaps))^2)
# rsq_fit <- 1-sum((master_df2$totcaps - pred_fit)^2)/sum((master_df2$totcaps - mean(master_df2$totcaps))^2)



















# ########################
# # Run with Tom's layers
# 
# ##########
# # Build a multiple GLMM with vars at their best scales... (classified NAIP vars)
# ###########
# 
# ####   Use backward stepwise selection to find best model (takes a while!)
# 
# form <- as.formula(sprintf("totcaps ~ %s + (1|Gd2)",paste(sprintf("%s_%snn",covar_df2$orig,r_bestscales),collapse="+") ))
# 
# form2 <- as.formula(sprintf("totcaps ~ (%s)^2 + (1|Gd2)",paste(sprintf("%s_%snn",covar_df2$orig,r_bestscales),collapse="+") ))
# 
# r_bestmod <- buildmer::buildglmmTMB(form2, ziformula = ~0,data=master_df2, dispformula = ~1, direction="backward",crit="AIC",
#                                   reduce.random=F,family=nbinom1(link = "log"))
# r_bestmod2 <- r_bestmod@model
# 
# r_bestmod2 <- glmmTMB(totcaps ~ r_Snd_small_127nn+ (1|Gd2), ziformula = ~0,data=master_df2, dispformula = ~1,    # rerun best model
#                                      family=nbinom1(link = "log"))
# 
# summary(r_bestmod2)   
# 
# MuMIn::r.squaredGLMM(r_bestmod2)    # r-squared not nearly as good...
# 
# tmp <- rownames(coef(summary(r_bestmod2))$cond)[-1]
# 
# topvars <- tmp[!grepl(":",tmp)]
# ints <- strsplit(tmp[grep(":",tmp)],":")
# 
# allvars <- c("totcaps","Gd2",unique(c(topvars,unlist(ints))))
# 
# i=1
# 
# #svg("PKM_main_effects3.svg",8,7)
# #layout(matrix(1:6,nrow=2,byrow = T))
# if(length(topvars)>0) tmp <- sapply(1:length(topvars), function(i) VisualizeRelation(data=master_df2,
#                                                                                      model=r_bestmod2,
#                                                                                      predvar=topvars[i],
#                                                                                      allvars = allvars) )     # run and save partial dependence plots for all top variables
# dev.off()
# 
# 
# 
# summary(bestmod2)
# 
# 





# ############
# # overlay predictions on orthophotos?
# 
# names(pred_fit) <- rownames(master_df2)
# 
# 
# dat$traplocs@data$pred <- NA
# ndx <- match(dat$traplocs@data$X.Det,as.numeric(names(pred_fit)))
# dat$traplocs@data$pred <- pred_fit[ndx]
# dat$traplocs@data$pred_rsc <- (dat$traplocs@data$pred-(min(dat$traplocs@data$pred,na.rm = T))+0.1)/(max(dat$traplocs@data$pred,na.rm=T)+0.2)
# 
# 
# okay <- !is.na(dat$traplocs@data$pred)
# graphics.off()
# cols <- numeric(length(dat$traplocs@data$X.Det))
# cols[] <- NA
# cols[okay] <- gray(1-dat$traplocs@data$pred_rsc[okay])
# plot(dat$traplocs,pch=20,cex=0.3,col=cols)
# 
# plot(1,1,col=gray(.9),pch=20,cex=3)
# 
# 
# ## add orthoimagery
# 
# library(ggmap)
# library(leaflet)
# library(RgoogleMaps)
# 
# ?RgoogleMaps::plotmap
# ?ggmap::get_map
# ?ggmap::ggmap
# 
# register_google("AIzaSyCsZjkVbEnK3p_DZPrzQ0OASi7TmtZwIi0",)
# 
# 
# tfm <- spTransform(dat$traplocs, CRS("+proj=longlat +datum=WGS84"))
# bbox <- c(left = min(tfm@coords[,1]), bottom =  min(tfm@coords[,2]), right =  max(tfm@coords[,1]), top =  max(tfm@coords[,2]))
# cord <- c(mean(c(bbox["left"],bbox["right"])),mean(c(bbox["top"],bbox["bottom"])))
# map <- ggmap::get_map(cord,maptype = "satellite",zoom=14,source="google")
# ggmap(map)+
#   geom_point(aes(x = tfm@coords[,1], y = tfm@coords[,2],color=pred_rsc),data=tfm@data) +  #color="darkred",
#              scale_colour_gradientn(colours = c("red","green"))
#   
# 
# 
# tfm <- spTransform(dat$traplocs[dat$traplocs@data$GdNo==3,], CRS("+proj=longlat +datum=WGS84"))
# bbox <- c(left = min(tfm@coords[,1]), bottom =  min(tfm@coords[,2]), right =  max(tfm@coords[,1]), top =  max(tfm@coords[,2]))
# cord <- c(mean(c(bbox["left"],bbox["right"])),mean(c(bbox["top"],bbox["bottom"])))
# map <- ggmap::get_map(cord,maptype = "satellite",source="google",zoom = 18)
# ggmap(map) +
#   geom_point(aes(x = tfm@coords[,1], y = tfm@coords[,2],color=pred_rsc),data=tfm@data) +  #color="darkred",
#   scale_colour_gradientn(colours = c("red","green"))
# 
# 
# 
# 
# 
# 
# plot(tfm)
# 
# 
# ############
# # Old code
# ############
# 
# 
# nn=maxdists3[2]
# 
# allmods <- list()
# 
# rsqs <- list()
# 
# for(nn in maxdists3){
#   
#   form <- as.formula(sprintf("totcaps ~ %s + 0",paste(covar_df2[[sprintf("std_%s",nn)]][ndx],collapse="+")))
#   
#   form2 <- as.formula(sprintf("totcaps ~ (%s)^2 + 0",paste(covar_df2[[sprintf("std_%s",nn)]][ndx],collapse="+")))
#   
#   x <- model.matrix(form2,master_df2)
#   
#   #temp <- glmnet(x,master_df$totcaps,family="poisson")
#   temp <- cv.glmnet(x,master_df2$totcaps,family="poisson",foldid = master_df2$Gd2)
#   
#   #plot(temp)
#   #opt = temp$lambda.min
#   coefs <- coef(temp, s = "lambda.min")
#   
#   bestcoefs <- rownames(coefs)[abs(coefs[,1])>0.000000001][-1]
#   
#   bestform <- as.formula(sprintf("totcaps ~ %s + (1|Gd2)",paste(bestcoefs,collapse="+")))
#   
#   allmods[[sprintf("mod_%snn",nn)]] <- glmmTMB(bestform,
#                                                master_df2,
#                                                ziformula = ~0 ,
#                                                dispformula = ~1,                     
#                                                family= nbinom1(link = "log"))
#   
#   rsqs[[sprintf("mod_%snn",nn)]] <- MuMIn::r.squaredGLMM(allmods[[sprintf("mod_%snn",nn)]])
#   
#   
# }
# 
# rsqs2 <- sapply(names(rsqs),function(t) rsqs[[t]]["trigamma",2] )   # top model has 6 nearest neighbors
# rsqs2
# aics <- sapply(names(allmods),function(t) AIC(allmods[[t]]))   # here, the top model has 1 nearest neighbor
# aics
# 
# bestmod2 <- allmods[[sprintf("mod_%snn",6)]]
# summary(bestmod2)
# 
# bestmod3 <- allmods[[sprintf("mod_%snn",1)]]
# summary(bestmod3)
# 




#master_df$Gd2 <- as.factor(master_df$Gd)
# temp <- glmmLasso(fix=form2,rnd=list(Gd2=~1), data=master_df, lambda = 100)    # note:: can add a smooth term if we want!

# form <- as.formula(sprintf("totcaps ~ %s",paste(covar_df2[[sprintf("std_%s",nn)]][ndx],collapse="+")))
# 
# glmnet()


# form2 <- as.formula(sprintf("totcaps ~ (%s)^2",paste(topvars,collapse = "+")))
# 
# bestmod2 <- buildmer::buildglmmTMB(form2, ziformula = ~0,data=master_df, dispformula = ~1, direction="backward",crit="AIC",
#                                   reduce.random=T,family=nbinom1(link = "log"))

### this one takes a while!

# summary(bestmod_6nn2)
# 
# topvars <- rownames(coef(summary(bestmod2))$cond)[-1]
# 
# formb <- bestmod_6nn2$call$formula
# 
# bestmod_6nn2 <- glmmTMB(formb,
#                      master_df,
#                      ziformula = ~0 ,
#                      dispformula = ~1,                    
#                      family= nbinom1(link = "log"))
# 
# 
# summary(bestmod_6nn2)

# nullmod <- glmmTMB(totcaps ~ 1,
#                    master_df,
#                    ziformula = ~0,
#                    dispformula = ~1 ,                                        # should we look at alternative dispersion formulas? Why just month?
#                    family= nbinom1(link = "log"))    # nbinom1(link = "log")  # dispformula = ~1,
# 
# 
# glance(nullmod)
# glance(bestmod_6nn2)
# 
# my_rsq(bestmod_6nn2)     # r-squared is 0.27
# 



# l1 <- as.numeric(logLik(bestmod_6nn2))
# l2 <- as.numeric(logLik(nullmod))
# 1-l1/l2                 # r-squared is 0.05. Null model is nearly as good as the best model!    Doing slightly better...

# AIC(nullmod,bestmod_6nn2)
# nullmod$obj$fn()

# 
# tmp <- rownames(coef(summary(bestmod3))$cond)[-1]
# 
# topvars <- tmp[!grepl(":",tmp)]
# ints <- strsplit(tmp[grep(":",tmp)],":")
# 
# allvars <- c("totcaps","Gd",unique(c(topvars,unlist(ints))))
# 
# i=1
# 
# svg("PKM_main_effects2.svg",5,7)
# layout(matrix(1:3,nrow=3,byrow = T))
# if(length(topvars)>0) tmp <- sapply(1:length(topvars), function(i) VisualizeRelation(master_df2,bestmod3,topvars[i],allvars = allvars) )     # run and save partial dependence plots for all top variables
# dev.off()
# 
# svg("PKM_interactions.svg",7,7)
# layout(matrix(1:2,nrow=1,byrow = T))
# if(length(ints)>0 ) tmp <- sapply(1:length(ints), function(i) VisualizeInteraction(master_df2,bestmod3,ints[[i]][2],ints[[i]][1],allvars=allvars) )
# dev.off()






###############
# Bolker r-2 code

# model<- allmods[[sprintf("mod_%snn",nn)]]
# 
# ## examples
# library(lme4)
# library(glmmTMB)
# 
# ## utility functions
# 
# ##' extract the 'conditional-model' term from a glmmTMB object;
# ##' otherwise, return x unchanged
# collapse_cond <- function(x)
#   if (is.list(x) && "cond" %in% names(x)) x[["cond"]] else x
# 
# ##' Cleaned-up/adapted version of Jon Lefcheck's code from SEMfit;
# ##' also incorporates some stuff from MuMIn::rsquaredGLMM.
# ##' Computes Nakagawa/Schielzeth/Johnson analogue of R^2 for
# ##' GLMMs. Should work for [g]lmer(.nb), glmmTMB models ...
# ##'
# ##' @param model a fitted model
# ##' @return a list composed of elements "family", "link", "marginal", "conditional"
# my_rsq <- function(model) {
#   
#   ## get basics from model (as generally as possible)
#   vals <- list(
#     beta=fixef(model),
#     X=getME(model,"X"),
#     vc=VarCorr(model),
#     re=ranef(model))
#   
#   ## glmmTMB-safety
#   if (is(model,"glmmTMB")) {
#     vals <- lapply(vals,collapse_cond)
#     nullEnv <- function(x) {
#       environment(x) <- NULL
#       return(x)
#     }
#     if (!identical(nullEnv(model$modelInfo$allForm$ziformula),nullEnv(~0)))
#       warning("R2 ignores effects of zero-inflation")
#     dform <- nullEnv(model$modelInfo$allForm$dispformula)
#     if (!identical(dform,nullEnv(~1)) &&
#         (!identical(dform,nullEnv(~0))))
#       warning("R2 ignores effects of dispersion model")
#   }
#   
#   # Test for non-zero random effects
#   if (any(sapply(vals$vc, function(x) any(diag(x)==0)))) {
#     ## FIXME: test more generally for singularity, via theta?
#     stop("Some variance components equal zero. Respecify random structure!")
#   }
#   
#   ## set family/link info
#   ret <- list()
#   if (is(model,"glmmTMB") || is(model,"glmerMod")) {
#     ret$family <- family(model)$family
#     ret$link <- family(model)$link
#   } else {
#     ret$family <- "gaussian"; ret$link <- "identity"
#   }
#   
#   ## Get variance of fixed effects: multiply coefs by design matrix
#   varF <- with(vals,var(as.vector(beta %*% t(X))))
#   
#   ## Are random slopes present as fixed effects? Warn.
#   random.slopes <- if("list" %in% class(vals$re)) {
#     ## multiple RE
#     unique(c(sapply(vals$re,colnames)))
#   } else {
#     colnames(vals$re)
#   }
#   if (!all(random.slopes %in% names(vals$beta))) 
#     warning("Random slopes not present as fixed effects. This artificially inflates the conditional R2. Respecify fixed structure!")
#   
#   ## Separate observation variance from variance of random effects
#   nr <- sapply(vals$re, nrow)
#   not.obs.terms <- names(nr[nr != nobs(model)])
#   obs.terms <- names(nr[nr==nobs(model)])
#   
#   ## Compute variance associated with a random-effects term
#   ## (Johnson 2014)
#   getVarRand <- function(terms) {
#     sum(
#       sapply(vals$vc[terms],
#              function(Sigma) {
#                Z <- vals$X[, rownames(Sigma), drop = FALSE]
#                Z.m <- Z %*% Sigma
#                return(sum(diag(crossprod(Z.m, Z))) / nobs(model))
#              } )
#     )
#   }
#   
#   ## Variance of random effects 
#   varRand <- getVarRand(not.obs.terms)
#   
#   if (is(model,"lmerMod") ||
#       (ret$family=="gaussian" && ret$link=="identity")) {
#     ## Get residual variance
#     varDist <- sigma(model)^2
#     varDisp <- 0
#   } else {
#     varDisp <- if (length(obs.terms)==0) 0 else getVarRand(obs.terms)
#     
#     badlink <- function(link,family) {
#       warning(sprintf("Model link '%s' is not yet supported for the %s distribution",link,family))
#       return(NA)
#     }
#     
#     if(ret$family == "binomial") {
#       varDist <- switch(ret$link,
#                         logit=pi^2/3,
#                         probit=1,
#                         badlink(ret$link,ret$family))
#     } else if (ret$family == "poisson" ||
#                grepl("nbinom",ret$family) ||
#                grepl("Negative Binomial", ret$family)) {
#       ## Generate null model (intercept and random effects only, no fixed effects)
#       
#       ## https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q4/023013.html
#       ## FIXME: deparse is a *little* dangerous
#       rterms <- paste0("(",sapply(findbars(formula(model)),deparse),")")
#       nullform <- reformulate(rterms,response=".")
#       null.model <- update(model,nullform)
#       
#       ## from MuMIn::rsquaredGLMM
#       
#       ## Get the fixed effects of the null model
#       null.fixef <- unname(collapse_cond(fixef(null.model)))
#       
#       ## in general want log(1+var(x)/mu^2)
#       logVarDist <- function(null.fixef) {
#         mu <- exp(null.fixef)
#         if (mu < 6)
#           warning(sprintf("mu of %0.1f is too close to zero, estimate may be unreliable \n",mu))
#         vv <- switch(ret$family,
#                      poisson=mu,
#                      nbinom1=,
#                      nbinom2=family(model)$variance(mu,sigma(model)),
#                      if (is(model,"merMod"))
#                        mu*(1+mu/getME(model,"glmer.nb.theta"))
#                      else mu*(1+mu/model$theta))
#         cvsquared <- vv/mu^2
#         return(log1p(cvsquared))
#       }
#       
#       varDist <- switch(ret$link,
#                         log=logVarDist(null.fixef),
#                         sqrt=0.25,
#                         badlink(ret$link,ret$family))
#     }
#   }
#   ## Calculate R2 values
#   ret$Marginal = varF / (varF + varRand + varDisp + varDist)
#   ret$Conditional = (varF + varRand) / (varF + varRand + varDisp + varDist)
#   return(ret)
# }
# 
# if (FALSE) {
#   fm1 <- lmer(Reaction~Days+(Days|Subject),data=sleepstudy)
#   fm2 <- glmmTMB(Reaction~Days+(Days|Subject),data=sleepstudy)
#   my_rsq(fm1)
#   my_rsq(fm2)
#   
#   ## devtools::install_github("jslefche/piecewiseSEM")
#   library(piecewiseSEM)
#   sem.model.fits(fm1)  ## same answer
#   
#   fm3 <- glmer(incidence/size~period+(1|herd),cbpp,
#                family=binomial,weights=size)
#   fm4 <- glmmTMB(incidence/size~period+(1|herd),cbpp,
#                  family=binomial,weights=size)
#   my_rsq(fm3)
#   my_rsq(fm4)
#   
#   fm5 <- glmer.nb(TICKS~YEAR+scale(HEIGHT)+(1|BROOD),grouseticks)
#   fm6 <- glmmTMB(TICKS~YEAR+scale(HEIGHT)+(1|BROOD),grouseticks,family=nbinom2)
#   my_rsq(fm5)
#   my_rsq(fm6)
# }
# 
# ## Tjur's coeff of determination, from sjmisc ... ????
# ## only does Bernoulli responses ???
# cod <- function(x) {
#   y <- model.response(model.frame(x))
#   pred <- predict(fm,type="response")
#   if (anyNA(rr <- residuals(x)))
#     pred <- pred[!is.na(rr)]
#   categories <- unique(y)
#   m1 <- mean(pred[which(y == categories[1])], na.rm = TRUE)
#   m2 <- mean(pred[which(y == categories[2])], na.rm = TRUE)
#   cod <- abs(m2 - m1)
#   return(cod)
# }
# 


