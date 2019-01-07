### Code for Prentice criteria proportion of treatment effect captured
### Coded by: Zoe Moodie
### 25AUG2014
### Modified by: Bhavesh Borate
### 30NOV2018

library(dplyr)

dataf <- read.csv("H:/CMV/CMV surrogate analysis/data/110318_weekly_vl_baseline.csv")
colnames(dataf) <- toupper(colnames(dataf))
dataf$wts <- 1

## Code to calculate bootstrap CI for pce 
# 'bPCE' returns the pce estimate based on the observed data and the 95% bootstrapped CI
# 'dataf' is the dataset
# 'iter' is the number of bootstrap iterations. 
# NOTE: Iterations to 11000 instead of 10000. This is because (depending on the seed) around 8 percent of bootstrapped datasets throw errors due to the dataset subsets having no records. (filtering on covariates and GCV results in datasets with no records!!) Such bootstrapped datasets are removed and the error caught without breaking the loop.   
# 'seed' is an integer for set.seed()
# DIS56 is indicator for CMV disease
# VL.0 is surrogate 
# GCV indicates Treatment with Gancyclovir
# AGVH and CMV.D are covariates

# To run the code
bPCE(dataf, iter=11000)


bPCE <- function(dataf, iter, seed=20181211){
  
  dataF <- dataf 
  
  # extract subsets of controls ('dataControls') and cases ('dataCases') to be used for resampling
  dataControls <- subset(dataF, DIS56==0)
  dataCases <- subset(dataF, DIS56==1)  
  
  if(!is.null(seed)){ set.seed(seed) }
  
  # the overall numbers of controls and cases for resampling
  nControls <- NROW(dataControls)
  nCases <- NROW(dataCases)
  
  bSampleControls <- matrix(sample(1:nControls, nControls*iter, replace=TRUE), nrow=nControls, ncol=iter)
  bSampleCases <- matrix(sample(1:nCases, nCases*iter, replace=TRUE), nrow=nCases, ncol=iter)
  
  # 'bPCE' is a list of PCE estimates
  bPCElist <- lapply(1:iter, function(j){
    # print(j)
    # create a bootstrap sample
    bdata <- rbind(dataControls[bSampleControls[,j],], dataCases[bSampleCases[,j],])
    
    tryCatch({
      pce.out <- pce(bdata)$pce
    }, error=function(e){cat("ERROR: Bootstrapped dataset was dropped as it had no records upon subsetting.\n")})
    
    return(list(pce=pce.out))
  })
  
  # combine all PCE estimates:
  bPCE <- drop(do.call(rbind, lapply(bPCElist,"[[","pce")))
  ci <- quantile(bPCE,probs=c(0.025,0.975))
  pce.est <- pce(dataf)$pce
  
  bList <- list(pce=pce.est,ci=ci)
  
  return(bList)
  
  rm(dataF, bdata)
}



pce <- function (dataf) {
  
  dat <- subset(dataf, !is.na('VL.0'))

  # Calculate PCE=proportion of captured treatment effect:
  fit <- glm(DIS56 ~ VL.0 + GCV + AGVH + CMV.D, data=dat, family="binomial")
  
  # sz indicates whether VL.0 from V or P should be used, agvh is 0 or 1, cmv.d is 0 or 1
  pY <- function(dat, fit, sz, gcv, agvh, cmv.d){
    datXV <- subset(dat, GCV==sz & AGVH==agvh & CMV.D==cmv.d)
    newdata <- data.frame(VL.0=datXV$VL.0, GCV=gcv, AGVH=agvh, CMV.D=cmv.d)
    return(predict(fit, newdata, type="response")) # vector 
  }
  
  dS <- function(dat, gcv, agvh, cmv.d){
    datVX <- subset(dat, GCV==gcv & AGVH==agvh & CMV.D==cmv.d)
    nvx <- sum(datVX$wts)
    return(datVX$wts/nvx) # vector
  }
  
  pX <- function(dataf, agvh, cmv.d){
    datafX <- subset(dataf, AGVH==agvh & CMV.D==cmv.d)
    propX <- nrow(datafX)/nrow(dataf)
    return(propX) # vector
  }
  
  cpi <- 0
  for (agvh in 0:1){
    for (cmv.d in 0:1){
        cpi <- cpi + sum(pY(dat, fit, sz=1, gcv=0, agvh, cmv.d)*dS(dat, gcv=1, agvh, cmv.d)) * pX(dataf, agvh, cmv.d) # here 'pX' is a scalar
    }
  }
  
  cpii <- 0
  for (agvh in 0:1){
    for (cmv.d in 0:1){
        cpii <- cpii + sum(pY(dat, fit, sz=0, gcv=0, agvh, cmv.d)*dS(dat, gcv=0, agvh, cmv.d)) * pX(dataf, agvh, cmv.d) # here 'pX' is a scalar
    }
  }
  
  cp <- cpi - cpii
  
  ncp <- 0
  for (agvh in 0:1){
    for (cmv.d in 0:1){
        ncp <- ncp + sum((pY(dat, fit, sz=1, gcv=1, agvh, cmv.d) - pY(dat,fit,sz=1,gcv=0, agvh, cmv.d))*dS(dat, gcv=1, agvh, cmv.d)) * pX(dataf, agvh, cmv.d) # here 'pX' is a scalar
     }
  }
  
  pce <- cp^2/(cp^2 + ncp^2)
  list(pce=pce,cp=cp,ncp=ncp,nvc=sum(dat$DIS56[dat$GCV==1]),npc=sum(dat$DIS56[dat$GCV==0]))
}






