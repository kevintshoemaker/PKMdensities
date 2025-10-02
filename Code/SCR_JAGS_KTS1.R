

# NOTE code needs to be updated with new spatal packages (no more sp and rgeos)

###############
# TODO
###############

# Add covariate for time of year on the detection probability. [don't do this]

# Look at density across time.
    # look at 2, 4, 6 vs 1, 3, 5
    # look at 1,2, vs 3,4 vs 5,6 



###############
# Questions
###############

## Q: were there any grids with zero PKM captures for one or more sessions? If so, I need to know which grids were surveyed for all sessions
### SEE ABOVE- THIS IS IMPORTANT!   [use DATE sheet to extract info on when grids were surveyed]

################
# Clear workspace
################

# rm(list=ls())

################
# Load packages
################

library(secr)
library(raster)
library(sp)
# library(rgdal)   # no longer available
library(abind)
# library(rgeos)  # no longer available
library(lubridate)
library(lunar)

library(sf)
library(terra)


################
# Load functions
################

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


##############
# Set User
##############

KEVIN <- TRUE
SARAH <- FALSE

##############
# Global vars
##############

N_GRIDS <- 18

N_SESSIONS <- 6

N_DENSITIES <- N_GRIDS*N_SESSIONS

N_DENSITIES



###################
# KTS: use JAGS to get density estimates

global <- list()

# global$CRS <- CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
global$CRS <- 26911

#########
# trap locations

traplocs <- read.csv("trap_locations.csv",header=T)

traplocs$Det <- traplocs$X.Det
traplocs$X.Det <- NULL
traplocs$X <- NULL
# traplocs$Gd <- dat$traplocs$GdNo

traplocs <- st_as_sf(traplocs,coords = c("Easting","Northing"),crs=global$CRS,remove=F)
plot(traplocs$geometry)


# names(traplocs)[1] <- "Det"
# 
# traplocs <- SpatialPointsDataFrame(traplocs[,c("Easting","Northing")],data=traplocs,proj4string = global$CRS)
# rownames(traplocs@data) <- traplocs@data$Det
# plot(traplocs)

traplocs$GdNo <- paste0("grid_",traplocs$GdNo)


##########
# read in capture data

PKMcaps <- read.csv("PKM_caps2.csv",header=T)
summary(PKMcaps)

names(PKMcaps)[1] <- "Session"
PKMcaps$Session <- paste0("session_",PKMcaps$Session)

PKMcaps$GdNo <- traplocs$GdNo[match(PKMcaps$Det,traplocs$Det)] 

head(PKMcaps)

dets <- sort(unique(traplocs$Det))
ndets <- length(dets)
ndets     # 1152 detectors

t=dets[1]


##########
# Dates

PKMdates <- read.csv("dates.csv",header=T)

PKMdates$Date <- mdy(PKMdates$Date)  # convert to date format

PKMdates$moon <- lunar::lunar.illumination(PKMdates$Date)   # add moon illumination


head(PKMdates)


#########
# make capture history: grid, year, individual, detector 

length(dets)
dets

allgrids <- sort(unique(PKMdates$Grid))
allgrids2 <- paste0("grid_",allgrids)
ngrids <- length(allgrids)    # 18 grids

allinds<-sort(unique(PKMcaps$ID)) #individuals trapped 
ninds <- length(allinds)   # 700 individuals trapped

allinds_pergrid <- lapply(1:ngrids,function(t){ temp<-subset(PKMcaps,GdNo==allgrids2[t]); sort(unique(temp$ID))}  )
names(allinds_pergrid) <- allgrids2
ninds_pergrid <- sapply(1:ngrids,function(t){ temp<-subset(PKMcaps,GdNo==allgrids2[t]); length(unique(temp$ID))}  )
names(ninds_pergrid) <- allgrids2
 #max_indspergrid <- max(ninds_pergrid)

allsessions <- sort(unique(PKMdates$Session))
allsessions2 <- paste0("session_",allsessions)
nsessions <- length(allsessions)

alloccasions <- sort(unique(PKMdates$Occ))

allsessions_pergrid <- lapply(1:ngrids,function(t){ temp<-subset(PKMdates,Grid==allgrids[t]); allsessions2[sort(unique(temp$Session))] }  )
names(allsessions_pergrid) <- allgrids2
nsessions_pergrid <- sapply(1:ngrids,function(t){ temp<-subset(PKMdates,Grid==allgrids[t]); length(unique(temp$Session)) }  )
names(nsessions_pergrid) <- allgrids2

allinds_pergridses <- lapply(1:ngrids,function(t){lapply(1:nsessions_pergrid[t], function(i){temp<-subset(PKMcaps,GdNo==allgrids2[t]&Session==allsessions_pergrid[[t]][i]); sort(unique(temp$ID)) })    }  )
temp <- lapply(1:ngrids,function(t) names(allinds_pergridses[[t]])<<-allsessions_pergrid[[t]] )
names(allinds_pergridses) <- allgrids2
ninds_pergridses <- lapply(1:ngrids,function(t){sapply(1:nsessions_pergrid[t], function(i){temp<-subset(PKMcaps,GdNo==allgrids2[t]&Session==allsessions_pergrid[[t]][i]); length(unique(temp$ID)) })    }  )
temp <- lapply(1:ngrids,function(t) names(ninds_pergridses[[t]])<<-allsessions_pergrid[[t]] )
names(ninds_pergridses) <- allgrids2

alldets_pergrid <- lapply(1:ngrids,function(t){ temp<-subset(traplocs,GdNo==allgrids2[t]); sort(unique(temp$Det))}  )
names(alldets_pergrid) <- allgrids2
ndets_pergrid <- sapply(1:ngrids,function(t){ temp<-subset(traplocs,GdNo==allgrids2[t]); length(unique(temp$Det))}  )
names(ndets_pergrid) <- allgrids2

Capt <- list()
temp <- sapply(allgrids2,function(t) Capt[[as.character(t)]] <<- list())
temp <- sapply(allgrids2,function(t) sapply(1:nsessions_pergrid[t], function(i) Capt[[as.character(t)]][[allsessions_pergrid[[t]][i]]] <<- matrix(0,nrow=ninds_pergridses[[t]][i],ncol=ndets_pergrid[t]) ) )

 ## Q: were there any grids with zero captures for one or more sessions? If so, I need to know which grids were surveyed for all sessions
   ### SEE ABOVE- THIS IS IMPORTANT!

g=11
for(g in 1:ngrids){
  thisgrid <- allgrids2[g]
  sessions <- allsessions_pergrid[[thisgrid]]
  s=1
  for(s in 1:nsessions_pergrid[g]){
    thisses <- sessions[s]
    rownames(Capt[[thisgrid]][[thisses]]) <- allinds_pergridses[[thisgrid]][[thisses]]
    colnames(Capt[[thisgrid]][[thisses]]) <-  alldets_pergrid[[thisgrid]]
    temp <- subset(PKMcaps,GdNo==thisgrid&Session==thisses)
    i=1
    for(i in 1:ninds_pergridses[[g]][s]){
      temp2 <- subset(temp,ID==allinds_pergridses[[thisgrid]][[thisses]][i])
      tabb <- table(temp2$Det) 
      tabnums <- as.numeric(tabb)
      inn <- names(tabb)
      Capt[[thisgrid]][[thisses]][i,colnames(Capt[[thisgrid]][[thisses]])%in%inn] <- tabnums
    }
  }
}

sapply(sapply(allgrids2,function(t) sapply(allsessions_pergrid[[t]],function(i)  max(Capt[[t]][[i]])  ) ),max)   # data check: max number of caps of the same individual per detector


###########
# spatial data processing

head(traplocs)
traplocs_bygrid <- list()
g=1
for(g in 1:ngrids){
  thisgrid <- allgrids2[g]
  temp <- subset(traplocs,GdNo==thisgrid)
  traplocs_bygrid[[thisgrid]] <- temp[,c("Easting","Northing")]
  
  # traplocs_bygrid[[thisgrid]] <- SpatialPointsDataFrame(temp[,c("Easting","Northing")],data=temp@data,proj4string = global$CRS)
}
  
plot(traplocs_bygrid$grid_15)

lapply(traplocs_bygrid,plot)

buf <- 40

# a <- traplocs_bygrid$grid_15 
# sf::st_buffer(a,dist=buf)

traplocs_buf <-  lapply(traplocs_bygrid,function(t) st_buffer(t,dist=buf)) #gBuffer(t,width=buf) )  # buffer the points
lapply(traplocs_buf,plot)

t <- traplocs_buf$grid_1
st_coordinates(t)[,"X"]

minx <- sapply(traplocs_buf,function(t)  min(st_coordinates(t)[,"X"])  ) #min(t@polygons[[1]]@Polygons[[1]]@coords[,1]))
maxx <- sapply(traplocs_buf,function(t)  max(st_coordinates(t)[,"X"])  ) #max(t@polygons[[1]]@Polygons[[1]]@coords[,1]))
miny <- sapply(traplocs_buf,function(t)  min(st_coordinates(t)[,"Y"])  ) #min(t@polygons[[1]]@Polygons[[1]]@coords[,2]))
maxy <- sapply(traplocs_buf,function(t)  max(st_coordinates(t)[,"Y"])  ) #max(t@polygons[[1]]@Polygons[[1]]@coords[,2]))

areafun <- function(x){
  (maxx[x]-minx[x])*(maxy[x]-miny[x])*0.0001
}

#t <- traplocs_bygrid$grid_1
traplist <- lapply(traplocs_bygrid,function(t) st_coordinates(t) )   # t@coords )

area_inha <- sapply(allgrids2,areafun)
names(area_inha) <- allgrids2



#########
# covariates?

temporalcovars <- read.csv("TemporalCovars2.csv",header=T)
summary(temporalcovars)
head(temporalcovars)
names(temporalcovars)

temp <- strsplit(temporalcovars$Sess.Gd.Occ.ID,split="_")
temporalcovars$Session2 <- as.numeric(sapply(temp,function(t) t[1]))
temporalcovars$Grid <- as.numeric(sapply(temp,function(t) t[2]))
temporalcovars$Occasion <- as.numeric(sapply(temp,function(t) t[3]))

### link date and covariates dataframes


temporalcovars$Sess.Gd.Occ.ID
PKMdates$SGOID <- paste0(PKMdates$Session,"_", PKMdates$Grid,"_",PKMdates$Occ) 

ndx <- match(PKMdates$SGOID,temporalcovars$Sess.Gd.Occ.ID)

cbind(PKMdates$Grid,temporalcovars$Grid[ndx])   # check

mintemp <- list()
maxtemp <- list()
precip <- list()
wind <- list()

PKMdates$mintemp <- temporalcovars$Min.Temp[ndx]
PKMdates$maxtemp <- temporalcovars$Max.Temp[ndx]
PKMdates$precip <- temporalcovars$Precip..mm.[ndx]
PKMdates$wind <- temporalcovars$Avg.Wind[ndx]


cor(PKMdates[,c("mintemp","maxtemp","moon","wind","precip")])  # check correlations

cor(PKMdates[,c("mintemp","wind")])

## standardize covars


allcovars <- c("mintemp","moon","wind","precip")   # all covars to try

temp <- lapply(allcovars,function(t) PKMdates[[paste0(t,"_std")]]<<-scale(PKMdates[[t]])   )

mintemp <- list()
#maxtemp <- list()
precip <- list()
wind <- list()
moon <- list()

g=1
for(g in 1:length(allgrids)){
  thisgrid <- allgrids[g]
  thisgrid2 <- sprintf("Grid_%s",thisgrid)
  temp <- subset(PKMdates,Grid==thisgrid)
  thissess <- sort(unique(temp$Session))
  wind[[thisgrid2]] <- matrix(0,nrow=length(thissess),ncol=length(alloccasions))
  mintemp[[thisgrid2]] <- matrix(0,nrow=length(thissess),ncol=length(alloccasions))
  #maxtemp[[thisgrid2]] <- matrix(0,nrow=length(thissess),ncol=length(alloccasions)) 
  precip[[thisgrid2]] <- matrix(0,nrow=length(thissess),ncol=length(alloccasions))
  moon[[thisgrid2]] <- matrix(0,nrow=length(thissess),ncol=length(alloccasions))
  s = 1
  for(s in 1:length(thissess)){
    thissess2 <- thissess[s]
    temp2 <- subset(temp,Session==thissess2)
    temp2 <- temp2[order(temp2$Occ),]
    wind[[thisgrid2]][s,] <- temp2$wind_std[1:length(alloccasions)]
    mintemp[[thisgrid2]][s,] <- temp2$mintemp_std[1:length(alloccasions)]
    #maxtemp[[thisgrid2]][s,] <- temp2$maxtemp_std[1:length(alloccasions)]
    precip[[thisgrid2]][s,]  <- ifelse(temp2$precip_std[1:length(alloccasions)]>0,1,0)
    moon[[thisgrid2]][s,] <- temp2$moon_std[1:length(alloccasions)]
  }
  
}




#########
# JAGS code

BUGSfilename <- "pkm_jagscode2.txt"

cat("
    
model{

  p0    ~ dunif(0,1)
  p0.l <- log(p0/(1-p0))
  sigma ~ dunif(1,100)   # in units of meters
  beta.wind ~ dnorm(0,1)
  beta.moon ~ dnorm(0,1)
  beta.precip ~ dnorm(0,1)
  beta.temp ~ dnorm(0,1)
  
  for(grid in 1:ngrids){
    for(session in 1:nsessions[grid]){
      psi[grid,session] ~ dunif(0,1)
      for(ind in 1:ninds[grid,session]){
        z[grid,session,ind] ~ dbern(psi[grid,session])                   # determine if this is a real individual
        hrx[grid,session,ind] ~ dunif(minx[grid],maxx[grid])             # activity center
        hry[grid,session,ind] ~ dunif(miny[grid],maxy[grid])                 
        for(trap in 1:ntraps[grid]){
          D[grid,session,ind,trap] <- pow(pow(hrx[grid,session,ind]-trapmat[grid,trap,1],2) + pow(hry[grid,session,ind]-trapmat[grid,trap,2],2),1/2)  #euclidean distance
          logit(basecap[grid,session,ind,trap]) <- p0.l + beta.moon*moon[grid,session] + beta.temp*temp[grid,session] + beta.precip*precip[grid,session] + beta.wind*wind[grid,session] 
          halfnorm[grid,session,ind,trap]  <- basecap[grid,session,ind,trap] * exp(-1*(pow(D[grid,session,ind,trap],2))/(2*pow(sigma,2) ))  
          p[grid,session,ind,trap]<-halfnorm[grid,session,ind,trap]  * z[grid,session,ind] #probability of catching an individual at a given burrow system
          y[grid,session,ind,trap]~dbin(p[grid,session,ind,trap], 3)    # 3 is max number of times an individual could be captured at a particular trap system (detector)
        }
      }
    }
  }
  
  for(grid in 1:ngrids){
    for(session in 1:nsessions[grid]){
      N[grid,session]<- sum(z[grid,session,1:ninds[grid,session]])
      Density[grid,session]<-N[grid,session]/A[grid] #derive density
    }
  }

}

",file = BUGSfilename)




##########
# Initial data prep for JAGS

  # trap coordinates
trapmat <- do.call(abind,c(traplist,along=0))
trapmat[5,,]

maxnsessions <- max(nsessions_pergrid)   # 6 is the max number of sessions
nsessions_forjags <- nsessions_pergrid

temp <- lapply(ninds_pergridses,function(t){ if(length(t)<maxnsessions){ c(t,rep(0,times=(maxnsessions-length(t))))}else{t}    } )
ninds_forjags <- do.call(abind,c(temp,along=0))
colnames(ninds_forjags) <- NULL
ninds_forjags_aug <- ninds_forjags
naug <- 50
temp <- lapply(1:length(allgrids),function(t) ninds_forjags_aug[t,1:nsessions_pergrid[t]] <<- ninds_forjags_aug[t,1:nsessions_pergrid[t]] + naug  )

#quick data check

unique(subset(PKMcaps,Session=="session_2")$GdNo)
unique(subset(PKMdates,Session==2)$Grid)   # lots of surveys done but very few grids with recorded captures for session 2

maxninds <- max(sapply(ninds_pergridses,max))   # 24 is the max number of individuals

maxndets <- max(ndets_pergrid)

caphist_forjags <- array(0,dim=c(ngrids,maxnsessions,maxninds+naug,maxndets))

g=4
for(g in 1:ngrids){
  thisgrid <- allgrids2[g]
  sessions <- allsessions_pergrid[[thisgrid]]
  s=5
  for(s in 1:nsessions_pergrid[g]){
    thisses <- sessions[s]
    this <- Capt[[thisgrid]][[thisses]]
    if(nrow(this)>0){
      caphist_forjags[g,s,1:ninds_pergridses[[thisgrid]][thisses],] <- this
    }
  }
}


### covariates

temp <- lapply(c("Grid_19","Grid_20"),function(t) if(nrow(wind[[t]])!=6) wind[[t]] <<- rbind(wind[[t]],matrix(0,nrow=2,ncol=3)) )
windmat <- do.call(abind,c(wind,along=0))

temp <- lapply(c("Grid_19","Grid_20"),function(t) if(nrow(precip[[t]])!=6) precip[[t]] <<- rbind(precip[[t]],matrix(0,nrow=2,ncol=3)) )
precipmat <- do.call(abind,c(precip,along=0))

temp <- lapply(c("Grid_19","Grid_20"),function(t) if(nrow(moon[[t]])!=6) moon[[t]] <<- rbind(moon[[t]],matrix(0,nrow=2,ncol=3)) )
moonmat <- do.call(abind,c(moon,along=0))

temp <- lapply(c("Grid_19","Grid_20"),function(t) if(nrow(mintemp[[t]])!=6) mintemp[[t]] <<- rbind(mintemp[[t]],matrix(0,nrow=2,ncol=3)) )
mintempmat <- do.call(abind,c(mintemp,along=0))


moonmat2 <- apply(moonmat,c(1,2),mean)
windmat2 <- apply(windmat,c(1,2),mean)
mintempmat2 <- apply(mintempmat,c(1,2),mean)
precipmat2 <- apply(precipmat,c(1,2),mean)

#########
# Initialization stuff

initx <- array(NA,dim=c(ngrids,maxnsessions,maxninds+naug))
inity <- array(NA,dim=c(ngrids,maxnsessions,maxninds+naug))

g=1
for(g in 1:ngrids){
  thisgrid <- allgrids2[g]
  sessions <- allsessions_pergrid[[thisgrid]]
  s=1
  for(s in 1:nsessions_pergrid[g]){
    thisses <- sessions[s]
    thisch <- Capt[[thisgrid]][[thisses]]
    i=1
    
    if(ninds_pergridses[[thisgrid]][thisses]>0){
      for(i in 1:ninds_pergridses[[thisgrid]][thisses] ){
        traploc1 <- which(thisch[i,]==1)[1]
        initx[g,s,i] <- trapmat[g,traploc1,1]
        inity[g,s,i] <- trapmat[g,traploc1,2]
      }
    }
    
    i=8
    for(i in (ninds_pergridses[[thisgrid]][thisses]+1):(ninds_pergridses[[thisgrid]][thisses]+naug) ){
      xrange <- c(min(trapmat[g,,1]),max(trapmat[g,,1])) 
      yrange <- c(min(trapmat[g,,2]),max(trapmat[g,,2])) 
      initx[g,s,i] <- runif(1,xrange[1],xrange[2])
      inity[g,s,i] <- runif(1,yrange[1],yrange[2])
    }
     
  }
}

zinit <- array(NA,dim=c(ngrids,maxnsessions,maxninds+naug))

g=11
for(g in 1:ngrids){
  thisgrid <- allgrids2[g]
  sessions <- allsessions_pergrid[[thisgrid]]
  s=1
  for(s in 1:nsessions_pergrid[g]){
    thisses <- sessions[s]
    zinit[g,s,1:(ninds_pergridses[[thisgrid]][thisses]+naug)] <- 1
  }
}

zinit[11,1,]

#########
# Final data prep for JAGS

data.for.bugs<-list(
  ngrids = ngrids,
  nsessions = nsessions_forjags,
  ninds = ninds_forjags_aug,
  ntraps = ndets_pergrid,
  minx = minx,
  miny = miny,
  maxx = maxx,
  maxy = maxy,
  trapmat = trapmat,
  y = caphist_forjags,
  moon = moonmat2,
  wind = windmat2,
  temp = mintempmat2,
  precip = precipmat2,
  A = area_inha
)

initz.bugs<-function(){
  list(
    p0=runif(1,0.2,0.5), 
    sigma=runif(1,20,30),
    hrx = initx,
    hry = inity,
    z=zinit
  )
}

# initz.bugs()#we has numbers!

#########
# Run JAGS  (takes a while to run with all grids and all sessions)

win<-jagsUI::jags(data=data.for.bugs,
                  inits=initz.bugs,
                  parameters.to.save=c("psi","p0","Density","sigma",
                                       "N","hrx","hry","z","beta.wind","beta.temp","beta.moon","beta.precip"),   
                  n.iter=2000, #times each chain is run
                  model.file=BUGSfilename, 
                  n.chains = 1,
                  n.adapt = 1000,
                  parallel = TRUE,
                  n.cores=3,
                  n.burnin = 500, #discards the first 'this many' iterations
                  n.thin= 1 # saves every nth iteration
)




##########
#  Store results


### Save results to disk

filename <- sprintf("JAGS_SCR_PKM1.RData")

# save(win, file = filename, envir = .GlobalEnv)


load(filename)


##############
# Visualize results

# c("psi","p0","Density","sigma",
# "N","hrx","hry","z","beta.wind","beta.temp","beta.moon","beta.precip")

results <- win$samples[[1]]

nMCMC <- length(results[,"p0"])

plot(results[,"psi[2,1]"])
plot(results[,"p0"])
mean(results[,"p0"]);quantile(results[,"p0"],c(0.025,0.975))
plot(results[,"sigma"])
mean(results[,"sigma"]);quantile(results[,"sigma"],c(0.025,0.975))
plot(results[,"Density[1,1]"])

plot(results[,"beta.wind"])    
plot(results[,"beta.temp"])
plot(results[,"beta.moon"])
plot(results[,"beta.precip"])

meandensities = NULL
meandens_grdses = list()     # TODO: supplemental figure of abundance over time for each grid and session

sit = 1
for(sit in 1:ngrids){
  thissit = allgrids2[sit]
  dens = NULL
  temp = NULL
  ssn = 2
  for(ssn in 1:nsessions_pergrid[thissit]){
    mcm = results[,sprintf("Density[%s,%s]",sit,ssn)]
    dens = rbind(dens,quantile(mcm,c(0.025,0.5,0.975)))
    temp = cbind(temp,mcm)
  }
  rownames(dens) = allsessions_pergrid[[thissit]]
  dens <- as.data.frame(dens)
  dens$Session <- as.numeric(gsub("session_","",rownames(dens)))
  
  Hmisc::errbar(dens$Session,dens$`50%`,dens$`97.5%`,dens$`2.5%`, ylim=c(0,30),ylab="PKM Density (per ha)",xlab="Session")
  title(thissit)
  
  temp2 = apply(temp,1,mean) 
  
  meandensities = rbind(meandensities,quantile(temp2,c(0.025,0.5,0.975)) )
  
}

rownames(meandensities) <- allgrids2

meandensities <- as.data.frame(meandensities)
meandensities <- meandensities[order(meandensities$`50%`,decreasing = T),]

mean(meandensities[,2])

graphics.off()

svg("meandensities.svg",6,4)
Hmisc::errbar(1:nrow(meandensities),meandensities$`50%`,meandensities$`97.5%`,meandensities$`2.5%`, ylim=c(0,20),ylab="PKM Density (per ha)",xlab="",xaxt="n")
axis(1,at=1:nrow(meandensities),labels = rownames(meandensities),las=2)
title("Mean PKM Densities")
dev.off()

write.csv(meandensities,"meandens_jags.csv",row.names = T)



##############
#  Visualize capture prob as function of covariates

mintemps <- seq(min(PKMdates$mintemp_std),max(PKMdates$mintemp_std),length.out = 100)
winds <- seq(min(PKMdates$wind_std),max(PKMdates$wind_std),length.out = 100) 

predp_wind <- matrix(0,ncol=length(winds),nrow=nMCMC)
predp_temp <- matrix(0,ncol=length(winds),nrow=nMCMC)

allp0 <- qlogis(results[,"p0"])
allbwind <- results[,"beta.wind"]
allbtemp <- results[,"beta.temp"]
i=1
for(i in 1:length(winds)){
  predp_temp[,i] <- plogis(allp0 + allbtemp*mintemps[i] )
  predp_wind[,i] <- plogis(allp0 + allbwind*winds[i])
}


qtemp <- as.data.frame(t(apply(predp_temp,2,function(t) quantile(t, c(0.025,0.5,0.975))  )))
qwind <- as.data.frame(t(apply(predp_wind,2,function(t) quantile(t, c(0.025,0.5,0.975))  )))

png("p0bywindandtemp2.png",5,6,units="in",res=600)
# svg("p0bywindandtemp2.svg",5,6)
par(mfrow=c(2,1))
par(mai=c(0.8,1,0.5,0.1))
Hmisc::errbar(mintemps,qtemp$`50%`,qtemp$`97.5%`,qtemp$`2.5%`, type="l", ylim=c(0,.2),
              ylab="Detection prob (p0)",xlab="Min Temperature (C)",xaxt="n",cap=0)  # 
lines(mintemps,qtemp$`50%`,lwd=2)
lines(mintemps,qtemp$`97.5%`,lwd=1)
lines(mintemps,qtemp$`2.5%`,lwd=1)
#axis(1,at=,labels = rownames(meandensities),las=2)
#title("PKM detection by temperature")
fig_label("a)",region="figure","topleft",cex=1.7)
rug(PKMdates$mintemp_std)
axis(1,at=seq(-3,3,0.5),labels = round(mean(PKMdates$mintemp)+seq(-3,3,0.5)*sd(PKMdates$mintemp)) )

Hmisc::errbar(winds,qwind$`50%`,qwind$`97.5%`,qwind$`2.5%`, type="l", ylim=c(0,0.2),
              ylab="Detection prob (p0)",xlab="Wind (m/s)",xaxt="n",cap=0)  # ,xaxt="n"
lines(winds,qwind$`50%`,lwd=2)
lines(winds,qwind$`97.5%`,lwd=1)
lines(winds,qwind$`2.5%`,lwd=1)
#axis(1,at=,labels = rownames(meandensities),las=2)
#title("PKM detection by wind")
fig_label("b)",region="figure","topleft",cex=1.7)
rug(PKMdates$wind_std)
axis(1,at=seq(-3,3,0.5),labels = round(mean(PKMdates$wind)+seq(-3,3,0.5)*sd(PKMdates$wind),1) )

dev.off()


#############
# illustrate dropoff in detection probability with distance

distances <- seq(0.1,40,length=100)

predp_dist <- matrix(0,ncol=length(winds),nrow=nMCMC)

allp0 <- qlogis(results[,"p0"])
allbwind <- results[,"beta.wind"]
allbtemp <- results[,"beta.temp"]
allsigma <- results[,"sigma"]


i=1
for(i in 1:length(winds)){
  basecap <- plogis(allp0)
  predp_dist[,i]  <- basecap * exp(-distances[i]^2/(2*allsigma^2 ))
}


qdist <- as.data.frame(t(apply(predp_dist,2,function(t) quantile(t, c(0.025,0.5,0.975))  )))

svg("pbydist2.svg",3,5)

par(mfrow=c(1,1))
par(mai=c(0.8,0.6,0.5,0.1))

Hmisc::errbar(distances,qdist$`50%`,qdist$`97.5%`,qdist$`2.5%`, type="l", ylim=c(0,0.1),
              ylab="Detection prob",xlab="Distance (m)",cap=0)  # ,xaxt="n"
lines(distances,qdist$`50%`,lwd=2)
lines(distances,qdist$`97.5%`,lwd=1)
lines(distances,qdist$`2.5%`,lwd=1)
#axis(1,at=,labels = rownames(meandensities),las=2)
#title("PKM detection by distance from activity center")
#rug(PKMdates$wind_std)

dev.off()









#############
# OLD CODE









##############
# Set Working Directory
##############

#if(SARAH) setwd('d:/density users/sarah hegg/March 2015')
#if(KEVIN) setwd("E:\\Dropbox\\PKM_Data\\Kevin\\secr files")
if(KEVIN) setwd("C:\\Users\\Kevin\\Dropbox\\My collaborations\\Marjorie Matocq\\Sarah Hegg\\PKM_Data\\Kevin\\secr files")

##############
# Read in data
##############

####################
# Detector data

temppkmdets <- read.table("AllDets_TestSpCovs.txt", header=FALSE)
head(temppkmdets)

# Trap arrangement
## modified code from Sarah H   [Q: do you get multiple PKMs in a trap?]
pkmdets <- read.traps(file="AllDets_TestSpCovs.txt", detector="multi", covnames=c("Gd","Trp",
                                                                                  "S_Snd","S_Crs","M_TCaps"))

head(pkmdets)  # note there is lot more in this object than just location though...

# clusterID(pkmdets) <- covariates(pkmdets)$Gd       
# clustertrap(pkmdets) <- covariates(pkmdets)$Trp
# table(clustertrap(pkmdets), clusterID(pkmdets))

plot(pkmdets)
summary(pkmdets)
##summary(mash(pkmdets))   #**not sure if this does anything? I get an error msg when I run it

##usage(pkmdets)   
# NULL


################
# CH data

pkmCHdata <- read.table("AllGridsPKM.txt", header=FALSE)

head(pkmCHdata)   # session, ID, occasion, trapID

tail(pkmCHdata)

# to avoid cross-grid movements, append grid name to the individual id

gridids <- as.character(attributes(pkmdets)$covariates$Gd)
trapids <- as.numeric(rownames(attributes(pkmdets)$covariates))

ndx <- match(pkmCHdata$V4,trapids)

pkmCHdata$V2 <- paste(gridids[ndx],pkmCHdata$V2,sep="")  # alter the individual IDs


# make the capture history
pkmCH <- make.capthist(capture=pkmCHdata, traps=pkmdets, fmt="trapID")

## doesn't work (secr bug?) but would not be used if it did because we care about grid-level covariates
# mash(pkmCH)    ## KTS: yeah I agree!

summary(pkmCH, terse = T)      # KTS: why does the number of detections go up so much with each occasion?

#               1    2    3    4    5    6
# Occasions     3    3    3    3    3    3
# Detections   99   39  149  393  245  430
# Animals      64   29  121  222  139  228
# Detectors  1152 1152 1152 1152 1152 1152

################
## habitat mask
################

cdtmask <- make.mask(pkmdets, type="trapbuffer", spacing = 10, buffer=50)   #  "clusterrect" ,spacing=10
nrow(cdtmask)   # not sure what the spacing should be...
plot(cdtmask)
spacing(cdtmask)

# make points layer with habitat mask
nrow(as.data.frame(cdtmask))

cdtmask_sp <- SpatialPoints(coords=as.data.frame(cdtmask))
plot(cdtmask_sp)

ext <- extent(cdtmask_sp)

ext <- raster::extend(ext,50)

cdtmask_raster <- raster(ext,resolution=50)
cdtmask_raster <- rasterize(cdtmask_sp,cdtmask_raster,field=1)
cdtmask_raster <- clump(cdtmask_raster, directions=8, gaps=FALSE)

grids <- extract(cdtmask_raster,cdtmask_sp)

plot(cdtmask_raster)
plot(cdtmask_sp,add=T)


covariates(cdtmask) <- data.frame(GrdName = as.factor(grids))

#?clump

## TODO: add habitat mask here... 



# too large      
# [1] 24.92188
#***"too large" is MGE's comment. Not sure why it is too large. I understand it to 
#the average distance between all traps. Given this is all traps over all grids, its not
#surprising that its larger than ~10 because some of the grids are 800m apart. Is this 
#a problem? Am I understanding it correctly?  


###########
# First test [not separated by session/grid combinations... ]
###########

initialsigma <- RPSV(pkmCH,CC=TRUE)   # something weird about session 3...

# # dist(pkmdets[pkmCH$`3`["M007",],])   # moved long distance
# # 
# # dist(pkmdets[pkmCH$`3`["C51568",],])  # even longer distance!   This seems to be why the global algorithm doesn't converge... 
# # 
# # dist(pkmdets[pkmCH$`3`["C51338",],])
# # 
# # dist(pkmdets[pkmCH$`3`["C51345",],])
# 
#   # NOTE: much slower with trap-level covariates!
# 
# test1 <- secr.fit(pkmCH, mask = cdtmask, model = list(D~GrdName*session, g0~S_Snd, sigma~1), start = list(g0=exp(-1.9),sigma=exp(2.5)), #mask = cdtmask, buffer=4*initialsigma$`6`, 
#                   verify = FALSE, trace=TRUE,detectfn = "HN",method="Nelder-Mead")    # sessioncov = sesscovdf, , GrdName*S_Snd
# 
# 
# test1
# 
# #pkm.covs.test <- secr.fit (pkmCH, mask = cdtmask, model = list(D~1,g0~1, sigma~1), verify = FALSE)
# 	#***This seems to work just fine, I've been skipping it to save time
# ##############################################################################
# ## after some trial and error, this arcane R code lets us split the data into 
# ## a separate 'session' for each of the 108 time x grid combinations    # KTS: why do we need this again?
# 
# grid <- lapply(X=covariates(traps(pkmCH)), FUN='[[', 'Gd')   # extract the grid names for each primary sampling period
# 
# CH <- mapply(split, pkmCH, grid, bytrap = TRUE, SIMPLIFY = FALSE)	# now 108 sessions... 
# 	#***I get warnings here: all of them because no PKM captures on some occasions at some grids
# 
#   # But, each new session is embedded within a list, one element per primary sampling period.
# 
# CH[[1]]$B
# 
# plot(CH[[1]]$B)   # visualize the first grid, first session
# 
#   #CH <- unlist(CH, recursive = FALSE)  # now put this back together in a single object...  # KTS: this didn't seem to work so I changed it
# 
# CH <- c(CH[[1]],CH[[2]],CH[[3]],CH[[4]],CH[[5]],CH[[6]])     # KTS: not sure this is right...
#   
# class(CH) <- c('list', 'capthist')  ## restore class of collapsed list of capthist objects
# 
# 
# 
# summary(CH, terse=T)
# # 1.A 1.B 1.C 1.D 1.E 1.F 1.G 1.H 1.I 1.J 1.K 1.L 1.M 1.N 1.O 1.P 1.Q 1.R 2.A 2.B 2.C 2.D 2.E 2.F 2.G 2.H
# # Occasions    3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3
# # Detections   6  14   7   1   5   9  10  19   3  15   0   1   2   0   0   7   0   0  14  18   0   0   0   0   0   0
# # Animals      6   7   6   1   3   6   6  11   3   8   0   1   2   0   0   4   0   0  10  16   0   0   0   0   0   0
# # Detectors   64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64
# # 2.I 2.J 2.K 2.L 2.M 2.N 2.O 2.P 2.Q 2.R 3.A 3.B 3.C 3.D 3.E 3.F 3.G 3.H 3.I 3.J 3.K 3.L 3.M 3.N 3.O 3.P
# # Occasions    3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3
# # Detections   0   0   0   0   0   1   0   6   0   0  17   7  18   7   9  11  16   2  12   1   1   0  12   7   6   8
# # Animals      0   0   0   0   0   1   0   2   0   0  14   6  16   4   8  11  15   2   8   1   1   0   7   3   5   8
# # Detectors   64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64
# # 3.Q 3.R 4.A 4.B 4.C 4.D 4.E 4.F 4.G 4.H 4.I 4.J 4.K 4.L 4.M 4.N 4.O 4.P 4.Q 4.R 5.A 5.B 5.C 5.D 5.E 5.F
# # Occasions    3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3
# # Detections   5  10  29  17  51  19  33  20  41  14  13   5  11   0  14  19  25  28  31  23   9  25  19   8  15  15
# # Animals      5  10  18  11  22  12  22  12  22  10  10   5   7   0   7   9  12  12  17  15   6  13  12   6  10   8
# # Detectors   64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64
# # 5.G 5.H 5.I 5.J 5.K 5.L 5.M 5.N 5.O 5.P 5.Q 5.R 6.A 6.B 6.C 6.D 6.E 6.F 6.G 6.H 6.I 6.J 6.K 6.L 6.M 6.N
# # Occasions    3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3
# # Detections  27  32  19   4   9   7   7  13   9   4  18   5  34  30  22  19  23  17  35  35  11  35  11  25  21  28
# # Animals     13  17   9   3   6   3   4   7   5   3  12   2  24  18   9  11  12   7  16  18   7  19   6  14  10  15
# # Detectors   64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64  64
# # 6.O 6.P 6.Q 6.R
# # Occasions    3   3   3   3
# # Detections  21  13  31  19
# # Animals      9   7  18  11
# # Detectors   64  64  64  64
# 
# 
# ## build dataframe of session covariates
# ## add columns for habitat as required
# sesscovdf <- expand.grid(Year = 1:3, Season = 1:2, Grid = levels(covariates(pkmdets)$Gd))
# sesscovdf$densid <- as.factor(1:nrow(sesscovdf))
# head(sesscovdf,12)
# nrow(sesscovdf)
# #    Year Season Grid
# # 1    1      1    A
# # 2    2      1    A
# # 3    3      1    A
# # 4    1      2    A
# # 5    2      2    A
# # 6    3      2    A
# # etc.
# 
# # #***Try adding habitat covariates...   KTS: I don't think this works, since habitat covars are at the 
# #  detector level and not the session level. I think we can access these covariates already anyway
# # sesscovdf2 <- expand.grid(Year = 1:3, Season = 1:2, Grid = levels(covariates(pkmdets)$Gd), 
# # 	Snd = covariates(pkmdets)$S_Snd)
# # head(sesscovdf2,12)
# # #  Year Season Grid      Snd
# # #1    1      1    A 92.87566
# # #2    2      1    A 92.87566
# # #3    3      1    A 92.87566
# # #4    1      2    A 92.87566
# # #5    2      2    A 92.87566
# # #6    3      2    A 92.87566	#***seemed to work at this point from what I can tell	 
# 
# 
# 
# 
# ## make a set of 108 masks!
# mask <- make.mask(traps(CH), type="trapbuffer", spacing=5,buffer = 50)  ## nx = 32 to save time  # nx = 32, 
# plot(mask[[2]])
# 
# nrow(mask[[1]])   # not sure what the spacing should be...
# spacing(mask[[1]])
# spacing(mask[[2]])
# spacing(mask[[3]])
# 
# 
# plot(traps(CH$B))
# plot(mask$B,add=T)
# 
# class(traps(CH[[1]]))
# 
# str(traps(CH[[1]]))
# 
# 
# ?make.mask
# 
# ## KTS: note- could just specify the buffer here, no need to make the mask...
# ##      note- this takes a very long time to run. Haven't gotten any errors yet, but wow I wasn't expecting it to be this slow
# ##      note- could probably specify D~session (I think this is the same thing...)
# ##      note- might be faster if we have fewer points in the habitat mask... we should look into this...
# ##      note- use simulate.secr to get better confidence intervals around density?
# ##      note- need to make sure initial parameter estimates are good
# ##      note- need to make sure that the buffer width is appropriate
# 
# 
# 
# 
# ## try it out
# test1 <- secr.fit(CH, mask = mask, model = list(D~densid, g0~S_Snd, sigma~1), start = list(D=1,g0=0.1,sigma=8), 
#                    sessioncov = sesscovdf, verify = FALSE)
# 
# getwd()
# 
# save(test1,file="Test1.RData")    
# 
# summary(test1)
# 
# ## KTS: test1 model seemed to run, although with a couple warnings...  took 26 hours on my computer
# 
# 
# ### Print out densities etc. ... 
# coef(test1)     # it ran! Looks like density estimates (here on log scale) were obtained for all grid/session combinations
# 
# test1     # this gives density estimates and sigma estimates on the real scale. 
# 
# ### get densities for export
# 
# a<-coef(test1)
# exp(a$beta[1:N_DENSITIES])
# 
# 
# esa.plot(test1)   # check the buffer width
# abline(v=40,lty=2,col="black")
# 
# plot(test1, limits = TRUE)    # takes longer with the "limits" argument
#  
# 
# ###
# 
# # ##*** now test with sand on g0
# # test2 <- secr.fit (CH, mask = mask, model = list(D~1, g0~Snd, sigma~1), start = list(D=1,g0=0.1,sigma=8), 
# #                    sessioncov = sesscovdf2, verify = FALSE)
# #***it doesn't like the added row for Snd covariate:
# 	#>Error in secr.design.MS(capthist, model, timecov, sessioncov, groups,  : 
#   	#>number of rows in 'sessioncov' should equal number of sessions
# 
# #####***Was not able to go further than this...
# 
# 
# 
# 
# ####################
# 
# ## extract estimates for the levels we want, not all 108 sessions:
# predict(test1, newdata = data.frame(Season = 1:2))
# 
# # $`Season = 1`
# # link   estimate SE.estimate        lcl        ucl
# # D       log  5.8114196 0.274412488  5.2979917  6.3746037
# # g0    logit  0.1434496 0.009753807  0.1253718  0.1636464
# # sigma   log 14.8551105 0.430047663 14.0358658 15.7221727
# # 
# # $`Season = 2`
# # link   estimate SE.estimate        lcl       ucl
# # D       log  5.8114196 0.274412488  5.2979917  6.374604
# # g0    logit  0.1033473 0.008321152  0.0881391  0.120832
# # sigma   log 14.8551105 0.430047663 14.0358658 15.722173
# 
# #############
# #Can we also extract estimates for each grid, not each season?
# predict(test1, newdata = data.frame(Grid = A))
# #gives this error...not sure how to recall each grid based on how he renamed them earlier..??
# #$> predict(test1, newdata = data.frame(Grid = A))
# #$Error in data.frame(Grid = A) : object 'A' not found









