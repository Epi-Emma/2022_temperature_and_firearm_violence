
################################################################
## DLNM - Maximum Daily Heat and Firearm Injuries (2-Day Lag) ##
################################################################

#Written by: Emma Gause
#Date: 02/18/22
#This code was adapted from example code provided by Antonio Gasparinni:
  # https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata


#Load in libraries
library("tidyverse")
library("dplyr")
library("dlnm")
library("tsModel")
library("mvmeta")
library("splines")
library("lubridate")
library("ggplot2")

#Create path to directory
datadir <- "[ADD DIRECTORY PATH, ENDING IN "/"]"

#Read in time series data - temp and injuries
dat <- readRDS(paste0(datadir, "GVA_Top100_Cities_2015_2020.rds"))

##------------------------------------------------------------------------##

## prepare the data ##

#create variable for year
dat$year <- year(dat$temp_date)

#create variable for month
dat$month <- month(dat$temp_date)

#create city label for named identification other than FIPS
dat <- dat %>% mutate(city_label = paste0(NAME, ", ", ST))

#create variable for day of week
dat$dow <- weekdays(dat$temp_date)

#Create NCA Region data for graphs
#create NCA crosswalk by state
st <- c("WA", "OR", "ID", 
        "CA", "NV", "AZ", "CO", "NM", "UT",
        "MT", "NE", "ND", "SD", "WY",
        "KS", "OK", "TX",
        "MN", "IA", "MO", "WI", "IL", "IN", "MI", "OH",
        "AR", "LA" ,"MS", "AL", "GA", "FL", "SC", "NC", "TN", "VA", "KY",
        "DC", "WV", "MD" ,"DE", "NJ", "NY", "PA", "CT", "RI", "MA", "NH", "VT", "ME")
nca_name <- c(rep("Northwest", 3), 
              rep("Southwest", 6), 
              rep("Great Plains", 8),
              rep("Midwest", 8), 
              rep("Southeast", 11), 
              rep("Northeast", 13)) 
nca <- cbind(st, nca_name)

#Now merge data by NCA categories
dat <- merge(dat, nca, by.x = "ST", by.y = "st", all.x = TRUE)
table(dat$nca_name, useNA = "ifany")


#get city metadata file
citynames <- dat %>% dplyr::select(STPLFIPS, ST, NAME, nca_name) %>% unique() %>% filter(!is.na(NAME))
citynames <- citynames %>% mutate(city_label = paste0(NAME, ", ", ST))
saveRDS(citynames, paste0(datadir, "city_metadata.rds"))

#Make lists of time series data for each city....
dat_list = split(dat, dat$city_label)

# Designate order so analysis and results match up
ord <- order(citynames$city_label)
dlist <- dat_list[ord]
citynames <- citynames[ord,]


##------------------------------------------------------------------------##

# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun = "bs"
vardegree = 2
varper <- c(10,75,90)

# SPECIFICATION OF THE LAG FUNCTION
#unnecessary since we've removed these arguments from the model, but a good reminder
lag <- 2
lagnk <- 1

# DEGREE OF FREEDOM FOR SEASONALITY
dfseas <- 7  ## Choice of 7 based on 2013 Bhaskaran (Int. J. Epidemiol.)

# COMPUTE TEMPERATURE PERCENTILES
per <- t(sapply(dat_list,function(x) 
  quantile(x$Tmax_C,c(2.5,10,25,50,75,90,97.5)/100,na.rm=T)))

# MODEL FORMULA
formula <- incidents~cb+dow+ns(temp_date,df=dfseas*length(unique(year)))

##------------------------------------------------------------------------##

################################################################################
# CREATE THE OBJECTS TO STORE THE RESULTS

# COEFFICIENTS AND VCOV FOR OVERALL CUMULATIVE SUMMARY
coef <- matrix(NA,nrow(citynames),length(varper)+vardegree,
               dimnames=list(citynames$city_label))
vcov <- vector("list",nrow(citynames))
names(vcov) <- citynames$city_label

################################################################################
# RUN THE LOOP

# LOOP
time <- proc.time()[3]
for(i in seq(length(dat_list))) {
  
  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data <- dat_list[[i]]
  
  # DEFINE THE CROSSBASIS
  argvar <- list(fun=varfun,knots=quantile(data$Tmax_C,varper/100,na.rm=T),
                 degree=vardegree)
  cb <- crossbasis(data$Tmax_C,lag=lag,argvar=argvar,
                   arglag=list(knots=logknots(lag,lagnk)))
  
  summary(cb)
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  model <- glm(formula,data,family=quasipoisson,na.action="na.exclude")
  cen <- mean(data$Tmax_C,na.rm=T)
  pred <- crosspred(cb,model,cen=cen)
  
  # REDUCTION TO OVERALL CUMULATIVE
  red <- crossreduce(cb,model,cen=cen)
  coef[i,] <- coef(red)
  vcov[[i]] <- vcov(red)
  
}
proc.time()[3]-time

##------------------------------------------------------------------------##

################################################################################
# MULTIVARIATE META-ANALYSIS OF THE REDUCED COEF AND THEN COMPUTATION OF BLUP
################################################################################

# CREATE AVERAGE TEMPERATURE AND RANGE AS META-PREDICTORS
avgTmax_C <- sapply(dat_list,function(x) mean(x$Tmax_C,na.rm=T))
rangeTmax_C <- sapply(dat_list,function(x) diff(range(x$Tmax_C,na.rm=T)))

################################################################################

# META-ANALYSIS
mv <- mvmeta(coef~avgTmax_C+rangeTmax_C,vcov,data=citynames,control=list(showiter=T))
summary(mv)
## I-square stat == 11.7%; Cochran Q-test for residual heterogeneity: p=0.02

################################################################################

# FUNCTION FOR COMPUTING THE P-VALUE OF A WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}

# TEST THE EFFECTS
fwald(mv,"avgTmax_C") #p=0.311
fwald(mv,"rangeTmax_C") #p=0.601

################################################################################

# OBTAIN BLUPS
blup <- blup(mv,vcov=T)

################################################################################

# PREDICT THE POOLED COEFFICIENTS
datanew <- data.frame(avgTmax_C=mean(tapply(avgTmax_C,citynames$city_label,mean)),
                      rangeTmax_C=mean(tapply(rangeTmax_C,citynames$city_label,mean)))

mvpred <- predict(mv,datanew,vcov=T,format="list")

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
predper <- c(0:98,seq(99,100,0.1))
tmeanall <- rowMeans(sapply(dat_list,function(x) quantile(x$Tmax_C, predper/100,na.rm=T)))

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS

# OBTAIN THE CENTERED PREDICTIONS
bvar <- onebasis(tmeanall,fun="bs",degree=2,
                 knots=if(is.null(varper)) NULL else tmeanall[paste(varper,".0%",sep="")])
cen <- tmeanall[50] #median

cp <- crosspred(bvar,coef=mvpred$fit,vcov=mvpred$vcov,model.link="log",
                at=tmeanall,cen=cen)


################################################################################
################################################################################
# RE-CENTERING

#USE THE MEDIAN TEMP

# GENERATE THE MATRIX FOR STORING THE RESULTS
medperccity <- medtempcity <- rep(NA,length(dat_list))
names(medtempcity) <- names(medperccity) <- citynames$city_label

# DEFINE MEDIAN TEMP:
for(i in seq(length(dat_list))) {
  data <- dat_list[[i]]
  medtempcity[i] <- median(data$Tmax_C)

}

# also get 90th perc temp for extreme heat AR calculation
perc90tempcity <- rep(NA,length(dat_list))
names(perc90tempcity) <- citynames$city_label
for(i in seq(length(dat_list))) {
  data <- dat_list[[i]]

  # getting 90th perc
  perc90tempcity[i] <- quantile(data$Tmax_C,90/100,na.rm=T)
}

##------------------------------------------------------------------------##

################################################################################
# COMPUTE THE ATTRIBUTABLE INCIDENTS FOR EACH CITY, WITH EMPIRICAL CI
################################################################################

# LOAD THE FUNCTION FOR COMPUTING THE ATTRIBUTABLE RISK MEASURES
source("[DIRECTORY PATH]/attrdl.R")

# CREATE THE VECTORS TO STORE THE TOTAL FIREARM INCIDENTS
totinc <- rep(NA,nrow(citynames))
names(totinc) <- citynames$city_label

# CREATE THE MATRIX TO STORE THE ATTRIBUTABLE INCIDENTS
#include components of total heat (moderate and extreme heat)
matsim <- matrix(NA,nrow(citynames),3,dimnames=list(citynames$city_label,
                                                 c("heat","mod_heat","ext_heat")))

# NUMBER OF SIMULATION RUNS FOR COMPUTING EMPIRICAL CI
nsim <- 1000

# CREATE THE ARRAY TO STORE THE CI OF ATTRIBUTABLE DEATHS
arraysim <- array(NA,dim=c(nrow(citynames),3,nsim),dimnames=list(citynames$city_label,
                                                              c("heat","mod_heat","ext_heat")))

################################################################################

# RUN THE LOOP
for(i in seq(dat_list)){
  
  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data <- dat_list[[i]]
  
  # DERIVE THE CROSS-BASIS
  argvar <- list(x=data$Tmax_C,fun=varfun,knots=quantile(data$Tmax_C,
                                                         varper/100,na.rm=T),degree=vardegree)
  cb <- crossbasis(data$Tmax_C,lag=lag,argvar=argvar,
                   arglag=list(knots=logknots(lag,lagnk)))
  
  # COMPUTE THE ATTRIBUTABLE INCIDENTS
  # THE REDUCED COEFFICIENTS ARE USED HERE (NECESSARY FOR "FORW")
  matsim[i,"heat"] <- attrdl(data$Tmax_C,cb,data$incidents,coef=blup[[i]]$blup,
                             vcov=blup[[i]]$vcov,type="an",dir="forw",cen=medtempcity[i],
                             range=c(medtempcity[i],100))
  matsim[i,"mod_heat"] <- attrdl(data$Tmax_C,cb,data$incidents,coef=blup[[i]]$blup,
                             vcov=blup[[i]]$vcov,type="an",dir="forw",cen=medtempcity[i],
                             range=c(medtempcity[i],perc90tempcity[i]))
  matsim[i,"ext_heat"] <- attrdl(data$Tmax_C,cb,data$incidents,coef=blup[[i]]$blup,
                             vcov=blup[[i]]$vcov,type="an",dir="forw",cen=medtempcity[i],
                             range=c(perc90tempcity[i],100))
  
  # COMPUTE EMPIRICAL OCCURRENCES OF THE ATTRIBUTABLE INCIDENTS
  # USED TO DERIVE CONFIDENCE INTERVALS
  arraysim[i,"heat",] <- attrdl(data$Tmax_C,cb,data$incidents,coef=blup[[i]]$blup,
                                vcov=blup[[i]]$vcov,type="an",dir="forw",cen=medtempcity[i],
                                range=c(medtempcity[i],100),sim=T,nsim=nsim)
  arraysim[i,"mod_heat",] <- attrdl(data$Tmax_C,cb,data$incidents,coef=blup[[i]]$blup,
                                vcov=blup[[i]]$vcov,type="an",dir="forw",cen=medtempcity[i],
                                range=c(medtempcity[i], perc90tempcity[i]),sim=T,nsim=nsim)
  arraysim[i,"ext_heat",] <- attrdl(data$Tmax_C,cb,data$incidents,coef=blup[[i]]$blup,
                                vcov=blup[[i]]$vcov,type="an",dir="forw",cen=medtempcity[i],
                                range=c(perc90tempcity[i],100),sim=T,nsim=nsim)
  
  # STORE THE DENOMINATOR OF ATTRIBUTABLE INCIDENTS, I.E. TOTAL OBSERVED INCIDENTS
  totinc[i] <- sum(data$incidents,na.rm=T)
}

################################################################################
# ATTRIBUTABLE NUMBERS

# CITY-SPECIFIC
ancity <- matsim
ancitylow <- apply(arraysim,c(1,2),quantile,0.025)
ancityhigh <- apply(arraysim,c(1,2),quantile,0.975)
rownames(ancity) <- rownames(ancitylow) <- rownames(ancityhigh) <- citynames$city_label

# TOTAL
antot <- colSums(matsim)
antotlow <- apply(apply(arraysim,c(2,3),sum),1,quantile,0.025)
antothigh <- apply(apply(arraysim,c(2,3),sum),1,quantile,0.975)

################################################################################
# TOTAL INCIDENTS IN DATA
totinctot <- sum(totinc)

################################################################################
# ATTRIBUTABLE FRACTIONS

# CITY-SPECIFIC
afcity <- ancity/totinc*100
afcitylow <- ancitylow/totinc*100
afcityhigh <- ancityhigh/totinc*100

# TOTAL ACROSS 100 CITIES
aftot <- antot/totinctot*100
aftotlow <- antotlow/totinctot*100
aftothigh <- antothigh/totinctot*100

##------------------------------------------------------------------------##

################################################################################
# CALCULATE MAXIMUM INCIDENT TEMPERATURE (BASED ON BLUPs)
################################################################################

################################################################################
# GENERATE THE MATRIX FOR STORING THE RESULTS
maxperccity <- maxtempcity <- rep(NA,length(dat_list))
names(maxtempcity) <- names(maxperccity) <- citynames$city_label

# DEFINE MAXIMUM INCIDENT VALUES: EXCLUDE 5% VERY LOW & VERY HOT TEMPERATURES AT TAILS
for(i in seq(length(dat_list))) {
  data <- dat_list[[i]]
  predvar <- quantile(data$Tmax_C,2.5:97.5/100,na.rm=T)
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  argvar <- list(x=predvar,fun=varfun,
                 knots=quantile(data$Tmax_C,varper/100,na.rm=T),degree=vardegree,
                 Bound=range(data$Tmax_C,na.rm=T))
  bvar <- do.call(onebasis,argvar)
  maxperccity[i] <- (1:99)[which.max((bvar%*%blup[[i]]$blup))]
  maxtempcity[i] <- quantile(data$Tmax_C,maxperccity[i]/100,na.rm=T)
}

##------------------------------------------------------------------------##

################################################################################
# DATA EXPORTS
################################################################################

################################################################################
# ATTRIBUTABLE FRACTION OVERALL
t(cbind(aftot,aftotlow,aftothigh))

# ATTRIBUTABLE NUMBER OVERALL
t(cbind(antot,antotlow,antothigh))

# ATTRIBUTABLE FRACTIONS BY CITY, W/ TOTAL INCIDENTS, MEDIAN AND MAX INCIDENT TEMPS
  #1 is all heat, 2 is moderate heat and 3 is extreme heat
(export <- cbind(totinc, medtempcity, afcity[,1],afcitylow[,1],afcityhigh[,1],
                 afcity[,2],afcitylow[,2],afcityhigh[,2],
                 afcity[,3],afcitylow[,3],afcityhigh[,3],maxperccity, maxtempcity))

#make it a dataframe for easier manipulation
export <- data.frame(export)
colnames(export) <- c("total_inc", "median_temp_c", "heat_ARf", "heat_ARfLB", "heat_ARfUB",
                      "heatmod_ARf", "heatmod_ARfLB", "heatmod_ARfUB",
                      "heatext_ARf", "heatext_ARfLB", "heatext_ARfUB", "max_inc_perc", "max_inc_temp_c")

#CONVERT C to F FOR TEMPS (WHY NOT?)
export$median_temp_f <- (export$median_temp_c*1.8)+32
export$max_inc_temp_f <- (export$max_inc_temp_c*1.8)+32

#write.csv(export, file = "[DIRECTORY PATH]/City_Results_2DayLag.csv",
#          na = "", row.names = TRUE)

##------------------------------------------------------------------------##

################################################################################
# PLOT
################################################################################

################################################################################
# Overall pooled association

indlab <- predper%in%c(0,1,10,50,90,99,100)

tiff("[DIRECTORY PATH]/Pooled_Assoc_2DayLag.tiff", 
     width=5.5,height=5, units = "in", res = 300)

plot(cp,ylab="RR",xlab="Temperature Percentile",axes=F,
     ylim=c(0.4,1.5),lwd=2,col=2,main="Pooled Association")
axis(1,at=tmeanall[indlab],labels=predper[indlab],cex.axis=0.9)
axis(2,cex.axis=0.9, at = c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6))
abline(v=c(tmeanall[c("10.0%","50.0%", "90.0%")]),
       lty=c(2,3,2))

dev.off()

################################################################################
# CITY-SPECIFIC BLUP GRAPHS W/ TEMPERATURE HISTOGRAM
  #ablines at 2.5, 97.5, 50 and maximum incident temperature (option to comment out)

xlab <- expression(paste("Temperature (",degree,"C)"))

pdf("[DIRECTORY PATH]/DLNM_2DayLag_City_Graphs.pdf",width=8,height=9)
layout(matrix(c(0,1,1,2,2,0,rep(3:8,each=2),0,9,9,10,10,0),ncol=6,byrow=T))
par(mar=c(4,3.8,3,2.4),mgp=c(2.5,1,0),las=1)

for(i in seq(length(dat_list))) {
  data <- dat_list[[i]]
  argvar <- list(x=data$Tmax_C,fun=varfun,degree=vardegree,
                 knots=quantile(data$Tmax_C,varper/100,na.rm=T))
  bvar <- do.call(onebasis,argvar)
  # CENTERING POINT IS AD MEDIAN TEMP OF EACH CITY
  pred <- crosspred(bvar,coef=blup[[i]]$blup,vcov=blup[[i]]$vcov,
                    model.link="log",by=0.1,cen=medtempcity[i])
  plot(pred,type="n",ylim=c(0,2),yaxt="n",lab=c(6,5,7),xlab=xlab,ylab="RR",
       main=citynames$city_label[i])
  ind1 <- pred$predvar<=medtempcity[i]
  ind2 <- pred$predvar>=medtempcity[i]
  lines(pred$predvar[ind1],pred$allRRfit[ind1],col=4,lwd=1.5)
  lines(pred$predvar[ind2],pred$allRRfit[ind2],col=2,lwd=1.5)
  mtext(citynames$nca_name[i],cex=0.7,line=0)
  axis(2,at=1:5*0.5)
  breaks <- c(min(data$Tmax_C,na.rm=T)-1,seq(pred$predvar[1],
                                            pred$predvar[length(pred$predvar)],length=30),
              max(data$Tmax_C,na.rm=T)+1)
  hist <- hist(data$Tmax_C,breaks=breaks,plot=F)
  hist$density <- hist$density/max(hist$density)*0.7
  prop <- max(hist$density)/max(hist$counts)
  counts <- pretty(hist$count,3)
  plot(hist,ylim=c(0,max(hist$density)*3.5),axes=F,ann=F,col=grey(0.95),
       breaks=breaks,freq=F,add=T)
  axis(4,at=counts*prop,labels=counts,cex.axis=0.7)
  abline(v=medtempcity[i],lty=3)
  abline(v=maxtempcity[i],lty=5, col = "grey")
  abline(v=c(per[i,c("2.5%","97.5%")]),lty=2)
}

dev.off()

##------------------------------------------------------------------------##