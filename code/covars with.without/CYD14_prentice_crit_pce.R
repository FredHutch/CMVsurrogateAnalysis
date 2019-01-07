### Code for Prentice criteria proportion of treatment effect captured
### Coded by: Zoe Moodie
### 25AUG2014

library(osDesign)
library(survival)

pce <- function (dataf,endpoint) {
  
  if(endpoint=="Overall"){
    #2phase data:    
    dat <- subset(dataf,!is.na(IMPSTLOG.AUCMB))
    
    flrstatus <- dat$ofstatus_m13 
    IMPSTLOG <- dat$IMPSTLOG.AUCMB
    nn0 <- nrow(subset(dataf,ofstatus_m13==0))
    nn1 <- nrow(subset(dataf,ofstatus_m13==1))
    
    group <- rep(1,length(flrstatus)) ####### group should indicate 2phase data: cases + controls with immuno data...
  }
  
  if(endpoint!="Overall"){
    
    if(endpoint=="Sero1"){i<-1}
    if(endpoint=="Sero2"){i<-2}
    if(endpoint=="Sero3"){i<-3}
    if(endpoint=="Sero4"){i<-4}
    
    #2phase data:       
    #exclude subjects with  missing serospecific titers (need to exclude for tps)
    dat <- dataf[!is.na(dataf$IMPSTLOG.AUCMB),]  
    dat <- na.omit(dat[,c("SUBJID","VACC","MALE","SEX","AGE","COUNTRY", paste("IMPSTLOG.Sero",i,sep=""),paste("s",i,"fstatus_m13",sep=""),"wts")])
    
    flrstatus <- dat[,paste("s",i,"fstatus_m13",sep="")]
    IMPSTLOG <- dat[,which( colnames(dat)==paste("IMPSTLOG.Sero",i,sep="") )]
    
    #control + non-sero-i cases, excluding subjects with missing casesero 
    nn0 <- nrow(subset(dataf,ofstatus_m0==0))+sum(subset(dataf,ofstatus_m13==1)[,paste("s",i,"fstatus_m13",sep="")]==0)  
    
    # sero-i-cases:
    nn1 <- sum(subset(dataf,ofstatus_m13==1)[,paste("s",i,"fstatus_m13",sep="")]==1)
    
    group <- rep(1,length(flrstatus)) ####### group should indicate 2phase data: cases + controls with immuno data...
  
    
  }
  
  dat$IMPSTLOG <- IMPSTLOG
  dat$flrstatus <- flrstatus
  
  #### Now calculate PCE=proportion of captured treatment effect:
  
  fit <- tps(flrstatus ~ IMPSTLOG+VACC+MALE+factor(AGE)+factor(COUNTRY),data=dat,nn0=nn0, nn1=nn1,group=group,method="PL",cohort=TRUE)
  
  # sz indicates whether IMPSTLOG from V or P should be used, 'male' is 0 or 1, 'age', 'country' are character strings
  pY <- function(dat, fit, sz, vacc, male, age, country){
    datXV <- subset(dat, VACC==sz & MALE==male & AGE==age & COUNTRY==country)
    newdataTps <- cbind(1, datXV$IMPSTLOG,vacc, male, as.numeric(age==">11"), as.numeric(age==">5,<=11"),
                        as.numeric(country=="MYS"), as.numeric(country=="PHL"), as.numeric(country=="THA"), as.numeric(country=="VNM"))
    return(tpsPredict(fit, newdataTps)) # vector    
  }
  
  dS <- function(dat, vacc, male, age, country){
    datVX <- subset(dat, VACC==vacc & MALE==male & AGE==age & COUNTRY==country)
    nvx <- sum(datVX$wts)
    return(datVX$wts/nvx) # vector
  }
  
  pX <- function(dataf, male, age, country){
    datafX <- subset(dataf, MALE==male & AGE==age & COUNTRY==country)
    propX <- nrow(datafX)/nrow(dataf)
    return(propX) # vector
  }
  
  cpi <- 0
  for (male in 0:1){
    for (age in levels(dat$AGE)){
      for (country in levels(dat$COUNTRY)){
        cpi <- cpi + sum(pY(dat, fit, sz=1, vacc=0, male, age, country)*dS(dat, vacc=1, male, age, country)) * pX(dataf, male,age,country) # here 'pX' is a scalar
      }
    }
  }
  
  cpii <- 0
  for (male in 0:1){
    for (age in levels(dat$AGE)){
      for (country in levels(dat$COUNTRY)){
        cpii <- cpii + sum(pY(dat, fit, sz=0, vacc=0, male, age, country)*dS(dat, vacc=0, male, age, country)) * pX(dataf, male,age,country) # here 'pX' is a scalar
      }
    }
  }
  
  cp <- cpi - cpii
  
  ncp <- 0
  for (male in 0:1){
    for (age in levels(dat$AGE)){
      for (country in levels(dat$COUNTRY)){
        ncp <- ncp + sum((pY(dat, fit, sz=1, vacc=1, male, age, country) - pY(dat,fit,sz=1,vacc=0,male,age,country))*dS(dat, vacc=1, male, age, country)) * pX(dataf, male,age,country) # here 'pX' is a scalar
      }
    }
  }
  
  pce <- cp^2/(cp^2 + ncp^2)
  
  list(pce=pce,cp=cp,ncp=ncp,nvc=sum(dat$flrstatus[dat$VACC==1]),npc=sum(dat$flrstatus[dat$VACC==0]))
  
}


## Code to calculate bootstrap CI for pce 

# 'bPCE' returns the pce estimate based on the observed data, the 95% bootstrapped CI, and the size of the bootstrapped immuno sets
# 'data' is assumed to be the ITT set (will sample cases and controls from those at risk at Mon13)
# 'endpoint' is one of "Overall", "Sero1", "Sero2", "Sero3", "Sero4"
# 'iter' is the number of bootstrap iterations
# 'seed' is an integer for set.seed()
bPCE <- function(dataf, endpoint, iter, seed=45378934){
  
  dataF <- dataf 
  
  if (endpoint!="Overall"){
    
    if(endpoint=="Sero1"){ i <- 1 }
    if(endpoint=="Sero2"){ i <- 2 }
    if(endpoint=="Sero3"){ i <- 3 }
    if(endpoint=="Sero4"){ i <- 4 }
    
    dataF$ofstatus_m0 <- dataF[,paste("s",i,"fstatus_m0",sep="")]
    dataF$ofstatus_m13 <- dataF[,paste("s",i,"fstatus_m13",sep="")]
    dataF$oftime_m13 <- dataF[,paste("s",i,"ftime_m13",sep="")]
    dataF$IMPSTLOG.AUCMB <- dataF[,paste("IMPSTLOG.Sero",i,sep="")]
    
  }
  
  # extract subsets of controls ('dataControls') and cases ('dataCases') to be used for resampling
  
  dataControls <- subset(dataF, ofstatus_m0==0)
  dataCases <- subset(dataF, ofstatus_m13==1)  
  
  if(!is.null(seed)){ set.seed(seed) }
  
  # the overall numbers of controls and cases for resampling
  nControls <- NROW(dataControls)
  nCases <- NROW(dataCases)
  
  bSampleControls <- matrix(sample(1:nControls, nControls*iter, replace=TRUE), nrow=nControls, ncol=iter)
  bSampleCases <- matrix(sample(1:nCases, nCases*iter, replace=TRUE), nrow=nCases, ncol=iter)
   
  # 'bPCE' is a list of PCE estimates
  #  and
  # 'bnI' the size of the bootstrapped immunogenicity set
  bPCElist <- lapply(1:iter, function(j){
    
    # create a bootstrap sample
    bdata <- rbind(dataControls[bSampleControls[,j],], dataCases[bSampleCases[,j],])
    # extract the bootstrapped immunogenicity set
  #  bdataI <- subset(bdata, !is.na(IMPSTLOG.AUCMB))
    # size of the bootstrapped immunogenicity set   
  #  bnI <- NROW(bdataI)
    
    pce.out <- pce(bdata,endpoint)$pce
    
  #  return(list(pce=pce.out, bnI=bnI))
  return(list(pce=pce.out))
  
  })
  
  # combine all PCE estimates:
  bPCE <- drop(do.call(rbind, lapply(bPCElist,"[[","pce")))
 # bnI <- do.call(c, lapply(bPCElist,"[[","bnI"))
  ci <- quantile(bPCE,probs=c(0.025,0.975))
  pce.est <- pce(dataf,endpoint)$pce
 
 # bList <- list(pce=pce.est,ci=ci, bnI=bnI)
  bList <- list(pce=pce.est,ci=ci)
  
  return(bList)
  
  rm(dataF,bdata)
}




