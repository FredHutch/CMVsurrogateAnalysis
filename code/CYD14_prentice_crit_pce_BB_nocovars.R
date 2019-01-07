### Code for Prentice criteria proportion of treatment effect captured
### Coded by: Zoe Moodie
### 25AUG2014
### Modified by: Bhavesh Borate
### 30NOV2018

library(dplyr)

dataf <- read.csv("../data/110318_weekly_vl_baseline.csv")
colnames(dataf) <- toupper(colnames(dataf))
dataf$wts <- 1
pce.out <- list(0)

## Code to calculate bootstrap CI for pce 
# 'bPCE' returns the pce estimate based on the observed data and the 95% bootstrapped CI
# 'dataf' is the dataset
# 'iter' is the number of bootstrap iterations. 
# NOTE: I have iterations to 11000 instead of 10000. This is because depending on the seed around 800ish or so throw errors due to the dataset being small. (filtering on covariates and GCV results in datasets with no records!!) I remove these bootstrapped datasets by catching them without breaking the loop.   
# 'seed' is an integer for set.seed()
# DIS56 is indicator for CMV disease
# VL.1 is surrogate 
# GCV indicates Treatment with Gancyclovir


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
    }, error=function(e){
      write(paste("ERROR: ", e$message), "logfile.csv", append= TRUE)
      })
    
    return(list(pce=pce.out))
  })
  
  # combine all PCE estimates:
  bPCE <- drop(do.call(rbind, lapply(bPCElist,"[[","pce"))) %>% 
    unlist(recursive = TRUE, use.names = TRUE)

  ci <- quantile(bPCE,probs=c(0.025,0.975))
  pce.est <- pce(dataf)$pce
  
  bList <- list(pce=pce.est,ci=ci)
  
  return(bList)
  
  rm(dataF, bdata)
}



pce <- function (dataf) {
  
  dat <- dataf %>% filter(!is.na(VL.1))
  
  # Calculate PCE=proportion of captured treatment effect:
  fit <- glm(DIS56 ~ VL.1 + GCV, data=dat, family="binomial")
  
  # sz indicates whether VL.1 from V or P should be used, agvh is 0 or 1, cmv.d is 0 or 1
  pY <- function(dat, fit, sz, gcv){
    datXV <- subset(dat, GCV==sz)
    newdata <- data.frame(VL.1=datXV$VL.1, GCV=gcv)
    return(predict(fit, newdata, type="response")) # vector 
  }
  
  dS <- function(dat, gcv){
    datVX <- subset(dat, GCV==gcv)
    nvx <- sum(datVX$wts)
    return(datVX$wts/nvx) # vector
  }
  
  # pX <- function(dataf){
  #   datafX <- subset(dataf, AGVH==agvh & CMV.D==cmv.d)
  #   propX <- nrow(datafX)/nrow(dataf)
  #   return(propX) # vector
  # }
  
  cpi <- sum(pY(dat, fit, sz=1, gcv=0)*dS(dat, gcv=1)) 
  cpii <- sum(pY(dat, fit, sz=0, gcv=0)*dS(dat, gcv=0)) 
  cp <- cpi - cpii
  ncp <- ncp + sum((pY(dat, fit, sz=1, gcv=1) - pY(dat,fit,sz=1,gcv=0))*dS(dat, gcv=1)) 
  
  pce <- cp^2/(cp^2 + ncp^2)
  list(pce=pce,cp=cp,ncp=ncp,nvc=sum(dat$DIS56[dat$GCV==1]),npc=sum(dat$DIS56[dat$GCV==0]))
}




# To run the code
bPCE(dataf, iter=10000, seed=20181211)






