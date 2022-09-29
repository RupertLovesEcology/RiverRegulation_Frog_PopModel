##############################################################################################################
##############################################################################################################
##############################################################################################################
## Short run attenuated code designed for sensitivity analysis using a Latin Hypercube
## Uses Large wetlands under regulated flow conditions at a sill height of 7.5 m
##############################################################################################################
# Nesting is 1000 runs then wet between the dry then duration of the dry then population size then T years  

## Variables for sensitivity
# clutch.size.lr <- 1000 to 6000
# hatch.dur <- 0.5 to 14
# hatch.pr <- 0.5 - .98
# tadpole.dur <- 50 to 90
# tadpole.mn <- 0.03 to 0.6
# ad.s.yr.mn <- 0.03 <- 0.35
# maxAge <- 5 to 10
# K.rel.tad.dem.vec[3] <- 500 to 5000
# K.rel.juv.dem.vec[3] <- 200 to 5000


###########################################################
## Aspatial southern bell frog sensitivity analysis code
##########################################################

## CJA Bradshaw & Rupert Mathwin
## June 2021

## Remove everything
rm(list = ls())

##Changed source file location below
## set source NOTE from Rversion 3.6.2 must be set to the U:drive on the flinders network and not to the local documents folder
## Subsequent versions choke and fail to connect to U, changed back to C to run smoothly (with Raph 12/7/2020)
source("C:/workspace/math0286/R/win-library/3.6/matrixOperators.r")
setwd("C:/workspace/math0286/R/")
.libPaths("C:/workspace/math0286/R/win-library/4")

## libraries
library(DescTools)
library(pracma)
library(DataCombine)
library(grid) #, lib.loc = "C:/Program Files/R/R-4.0.2/library")
library(doSNOW)
library(iterators)
library(snow)
library(foreach)
library(lhs)
library(data.table)
library(dismo)
library(gbm)
library(VGAM)
library(raster)



## Set up parallel processing (nproc is the number of processing cores to use)
nproc <- 6
cl.tmp = makeCluster(rep('localhost', nproc), type='SOCK')
registerDoSNOW(cl.tmp)
getDoParWorkers()


## Encapsulate the model within a simulation function to run in parallel once the LHC is complete
sbf_sim <- function(input, dir.nm, rowNum) {
  
  ## assign all parameter values
  for (d in 1:ncol(input)) {assign(names(input)[d], input[,d])}
  
    ####################################################
    ## necessary input calculations for stochastic model
    ####################################################
  
    #######################################################################################
    ## Functions
    # beta distribution shape parameter estimator function
    ##Generates an alpha and a beta value to inform the beta distribution 
    estBetaParams <- function(mu, var) {
      alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
      beta <- alpha * (1 / mu - 1)
      return(params = list(alpha = alpha, beta = beta)) }
    
    # complementary log-log
    cloglog <- function(x) log(-log(1-x))
    
    ##  Not used in this Southern Bell Frog Model
    AICc.glm <- function(...) {
      models <- list(...)
      num.mod <- length(models)
      AICcs <- numeric(num.mod)
      ns <- numeric(num.mod)
      ks <- numeric(num.mod)
      AICc.vec <- rep(0,num.mod)
      for (i in 1:num.mod) {
        if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
        if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
        AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
        ns[i] <- n
        ks[i] <- k
        AICc.vec[i] <- AICcs[i]}
      return(AICc.vec)}
    
    
    ####################################################################################
    ## Variables for sensitivity
    ####################################################################################
    
    ##                    To vary the model run
    # Number of iterations
    iter <- 1000
    itdiv <- iter/10
    
    ## The output for use in the Latin HyperCube will be population survival to gen 80 stored in the outPut array
    #rows are 1) total runs, 2) surviving runs, 3) can add at the end for % Sv or similar
    outPut <- array(data = 0, dim = 3 * iter)
    dim(outPut) <-c(3,iter)
    
    ## Clutch size could reasonably range from 1000 (a low estimate) - 6000 (a large clutch for a closely-related L. aurea)
    #clutch.size.lr <- c(1885,3893) # L. raniformis
    clutch.size.lr.sd.prop <- 502/mean(c(1885,3893))
    
    #duration of egg stage. Could range from 0.5 to 14 days? Some spp can go up to weeks and months but not this sp
    #hatch.dur <- c(2,4)
    hatch.dur.sd.prop <- 0.5/3
    
    ## survival data: Probability of hatch
    ## Could range from 
    #hatch.pr <- c(0.933, 1) # 2-4 days, L. raniformis
    
    ## tadpole duration could range from 50 to 90 days
    #tadpole.dur <- c(70,80) # duration (days) L. raniformis 23 deg
    tadpole.dur.sd.prop <- 2.5/75
    
    ## tadpole survival to metamorphosis could be as low as 5% or as high as 60%
    ##the below uses survival TO METAMORPHOSIS figures from Bull (C.signifera)
    tadpole.mn.1 <- mean(c(.15,.26))
    tadpole.sd.1 <- ((.26-tadpole.mn.1)+(tadpole.mn.1-.15))/2/1.96
    tadpole.mn.2 <- mean(c(.07,.56))
    tadpole.sd.2 <- ((.56-tadpole.mn.1)+(tadpole.mn.1-.07))/2/1.96
    #tadpole.mn <- mean(c(tadpole.mn.1, tadpole.mn.2))
    tadpole.sd <- sqrt(tadpole.sd.1^2 + tadpole.sd.2^2)
    tp.s.alpha <- estBetaParams(tadpole.mn, tadpole.sd/10)$alpha
    tp.s.beta <- estBetaParams(tadpole.mn, tadpole.sd/10)$beta
    
    #Adult annual survival from Pickett's PhD on L. AUREA
    ## These could reasonable range between means of 5% to 35%  ?
    #ad.s.yr.mn <- 0.2172 # Litoria aurea (Pickett et al.)
    ad.s.yr.sd <- 0.087 # Litoria aurea (Pickett et al.)
    
    ## SBF in this system breed and die during their <maxAge> year (keeping in mind we have a 0 yr)
    maxAge <- 5 
    #This could range from 5 (oldest recorded from skeletochronology) to 10 (oldest in captivity)
    
    
    ########################################################################
    
    WDTracker <- ""
    
    ## set time limit for projection in 1-yr increments
    yr.now <- 2020
    yr.end <- 2020 + 86
    
    #************************
    ##Note I have constrained t to 86 years
    t <- (yr.end - yr.now)
    yrs <- seq(yr.now,yr.end,1)
    longev <- 5
    age.vec <- seq(0,longev,1)
    lage <- length(age.vec)
    sex.ratio <- 0.5
    stages <- 4
    ## set population storage matrices n.mat
    n.mat <- array(data = 0, dim = ((maxAge + 1) * (t + 1) * 4))
    dim(n.mat) <- c((maxAge + 1),(t + 1),4)
    popmat <- matrix(0,nrow=maxAge+2,ncol=maxAge + 1)
    #colnames(popmat) <- age.vec[1:stages]
    #rownames(popmat) <- age.vec[1:stages]
    #rownames(popmat) <- c("Fecundity","Survival to 1 year","Survival to 2 year","Survival to 3 year","Survival to 4 year","Survival to 5 year")
    #,"Survival to 4 year","Survival to 5 year")
    
    
    ## fertility data 
    #clutch.size.la <- c(4124,6178) # L. aurea
    
    prop.breeding <- c(0,rep(1,maxAge))
    fert.mn <- mean(clutch.size.lr)*prop.breeding
    fert.mn
    
    
    tadpole.s.cs <- 0.10 # C. signifera
    tadpole.s.ra <- 0.05 # R. aurora
    ad.s.daily <- c(0.975, 0.996) # daily survival L. raniformis
    
    
    #tomet.dur.iter <- round(sum(c(runif(1, hatch.dur[1], hatch.dur[2]), runif(1, tadpole.dur[1], tadpole.dur[2]))), 0)
    tomet.dur.iter <- round((rnorm(1, hatch.dur, hatch.dur.sd.prop*hatch.dur) + rnorm(1, tadpole.dur, tadpole.dur.sd.prop*tadpole.dur)), 0)
    
    toad.dur.iter <- 365 - tomet.dur.iter
    #tomet.s.iter <- (runif(1, min=hatch.pr[1], max=hatch.pr[2])) * (runif(1, min=tadpole.s.ra, max=tadpole.s.cs)) 
    
    tomet.dur.mn <- round(sum(c(mean(hatch.dur), mean(tadpole.dur))), 0)
    toad.dur.mn <- 365 - tomet.dur.mn
    
    #tomet.s.mn <- mean(hatch.pr) * mean(c(tadpole.s.ra, tadpole.s.cs))
    #juv.s.daily.vec <- runif(toad.dur.iter, min=ad.s.daily[1], max=ad.s.daily[2])
    #juv.s.daily.mn <- rep(mean(ad.s.daily), toad.dur.mn)
    #ad.s.iter <- prod(runif(365, min=ad.s.daily[1], max=ad.s.daily[2]))
    #ad.s.mn <- (mean(ad.s.daily))^365
    
    #NEW ADULT ANNUAL SURVIVAL VALUES FROM PICKETT L. AUREA
    #ad.s.yr.mn <- 0.2172 # Litoria aurea (Pickett et al.)
    ad.s.yr.sd <- 0.087 # Litoria aurea (Pickett et al.)
    
    ##FROM PICKETT
    ad.s.yr.alpha <- estBetaParams(ad.s.yr.mn, ad.s.yr.sd/10)$alpha
    ad.s.yr.beta <- estBetaParams(ad.s.yr.mn, ad.s.yr.sd/10)$beta
    
    ##Calculate our survivals
    
    #Likelihood of egg hatching and tadpole surviving to metamorphosis
    tomet.s.iter <- (rbeta(1, tp.s.alpha, tp.s.beta)) * (runif(1, min=0.933, max=1))
    
    #sample a daily probability of survival
    toad.daily.s.iter <- nthroot(rbeta(1,ad.s.yr.alpha, ad.s.yr.beta) , 365)
    
    toad.s.season.iter <- toad.daily.s.iter ^ toad.dur.mn
    ## change the above line to the below line just need the syntax to resample
    ## toad.s.season.iter <- prod(runif(toad.dur.mn, toad.daily.s.iter))
    toad.s.iter <- tomet.s.iter * toad.s.season.iter                     
    
    # TO RESAMPLE annual survival of an adult
    #  ad.s.iter <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)
    
    #Create the survival vector, to adult survival then 4 adult survivals
    ad.s.vec.iter <- rep(NA,maxAge-2)
    for (s in 1:maxAge-1) {
      ad.s.vec.iter[s] <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)  }
    surv.iter <- c(toad.s.iter, ad.s.vec.iter,0)

    # We hereafter consider 4 population sizes, each represented with 4 density feedback values (one for each of egg, tad, juve, adult)
    # Small, medium, large, very large.
    # Assign the number of spawning masses that each site holds and convert it to the number of female egges (uses mean clutch size)
    K.rel.egg.dem.vec <- c(50, 150, 500, 1000)
    K.rel.egg.dem.vec <- K.rel.egg.dem.vec * 1444
    
    # tadpoles surviving to 1 year old will inhibit themselves NOTE THIS EQUATION ONLY CONSIDERS FEMALES SO THE ACTUAL NUMBER WILL BE DOUBLE  
    K.rel.tad.dem.vec <- c(100, 300, 800, 1000)
    #K.rel.tad.dem <- K.rel.tad.dem.vec[3]
    
    # number of first years present will influence survival from 1 to 2 years old but more strongly than for older age brackets         
    K.rel.juv.dem.vec <- c(200, 600, 1800, 3000)
    #K.rel.juv.dem <- K.rel.juv.dem.vec[3]
    
    # number of first years present will influence adult survival NOTE THIS EQUATION ONLY CONSIDERS FEMALES SO ACTUAL NUMBER WILL BE DOUBLE         
    K.rel.adult.dem.vec <- c(1000000,1000000,1000000,1000000)
    
    initSpawn <- c(10,30,85,150)
    
    #set founding fem pop sizes 
    init.fem.pop <- c(20,65,150,500) 
    
    init.vec <- rep(NA,((maxAge+1) * 4))
    dim(init.vec) <- c((maxAge+1),4)
    
    ##Populate the Matix (popmat) and create a failure matrix(popmat.fail) for years with no breeding
    ## removed  diag(popmat[2:(stages), ]) <- surv.mn[-stages]
    stages <- maxAge + 2
    diag(popmat[2:(maxAge+2), ]) <- surv.iter[-stages]
    #   diag(popmat[2:(stages), ]) <- surv.iter[-stages]
    #popmat[5,4] <- 0 # surv.mn[stages] 
    popmat[1,] <- fert.mn * sex.ratio
    popmat.orig <- popmat ## save original matrix as popmat.orig
    popmat <- popmat.orig
    
    #Create init.vec note 2nd dimension is 1 = S, 2 = M, 3 = L, 4 = XL
    totalSurv <- sum(popmat[2:(maxAge+2),1:(maxAge+1)])
    for (size in 1:4) {
      for (x in 1:(maxAge+1)) { init.vec[x,size] <- 0 }
      #init.vec[1,size] <- round(initSpawn[size] * (runif(1, min=clutch.size.lr[1], max=clutch.size.lr[2])))
      init.vec[1,size] <- round(initSpawn[size] * (rnorm(n=1, mean=clutch.size.lr, sd=clutch.size.lr.sd.prop*clutch.size.lr)))
      for (dd in 2:(maxAge+1)) { init.vec[dd,size] <- round((sum(popmat[2:(maxAge+1),(dd - 1)])/totalSurv) * (init.fem.pop[size]))    } }
    
    for (size in 1:4) { n.mat[,1,size] <- init.vec[,size] }
    
    # Create popmat.all [row,column,population size] to hold (density dependence * popmat) for each population size (uses popmat.orig as a placeholder)
    #popmat.all <- matrix(0,nrow=(4 * stages),ncol=stages)
    #dim(popmat.all) <- c(stages,stages,4)
    #for (size in 1:4) {
    #  colnames(popmat.all[,,size]) <- age.vec[1:stages]
    #  rownames(popmat.all[,,size]) <- c("Fecundity","Survival to 1 year","Survival to 2 year","Survival to 3 year","Survival to 4 year","Survival to 5 year")
    #  popmat.all[,,size] <- popmat.orig }
    popmat.current <- popmat
    popmat.fail <- popmat
    popmat.fail[1,] <- 0 
    
    #create the egg density correcting function and variables
    #   line equation y = 1.01 - ((x ^ (9/x)) * (x ^ 2)))
    surv.mult.egg.up <- 1.0
    surv.mult.egg.upmid <- 1.0
    surv.mult.egg.mid <- 0.94
    surv.mult.egg.lo <- 0.8
    surv.mult.egg.lo.lo <- 0.2
    surv.mult.egg.lo.lo.lo <- 0.05
    
    K.egg.up <- 1
    K.egg.upmid <- 0.98
    K.egg.mid <- 0.90
    K.egg.lo <- 0.7
    K.egg.lo.lo <- 0.3
    K.egg.lo.lo.lo <- 0.01
    
    K.egg.vec <- c(K.egg.up,K.egg.upmid, K.egg.mid,K.egg.lo, K.egg.lo.lo, K.egg.lo.lo.lo)
    surv.mult.egg.vec <- rev(c(surv.mult.egg.up, surv.mult.egg.upmid, surv.mult.egg.mid, surv.mult.egg.lo, surv.mult.egg.lo.lo, surv.mult.egg.lo.lo.lo))
    #(K.egg.vec, surv.mult.egg.vec, pch=19)
    ##plot(K.egg.vec, surv.mult.egg.vec, type="n")
    DD.dat <- data.frame(K.egg.vec, surv.mult.egg.vec)
    param.init <- c(1.01, 9, 2)
    fit.expd.egg <- nls(surv.mult.egg.vec ~ (a - ((K.egg.vec ^ (b/K.egg.vec)) * (K.egg.vec ^ d))),
                        algorithm = "port",
                        start = c(a = param.init[1], b = param.init[2], d = param.init[3]),
                        trace = TRUE,      
                        nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
    K.pred.egg.vec <- seq(0.01, 1, 0.01)
    pred.surv.egg.mult <- (as.numeric(coef(fit.expd.egg)[1]) - (K.pred.egg.vec ^ (as.numeric(coef(fit.expd.egg)[2])/K.pred.egg.vec) * (K.pred.egg.vec ^ (coef(fit.expd.egg)[3]))))
    #lines(K.pred.egg.vec, pred.surv.egg.mult, lty=2, lwd=1, col="red")
    
    eggfunc <- function(x) (1.01 - ((x ^ (9/x)) * (x ^ 2)))
    intArea <- integrate(eggfunc,lower=0,upper=0.99)
    area1 <- (as.numeric(intArea[1])/0.99)
    s.mult.egg.iter <- (as.numeric(coef(fit.expd.egg)[1]) - (0.999 ^ (as.numeric(coef(fit.expd.egg)[2])/0.999) * (0.999 ^ (coef(fit.expd.egg)[3]))))
    
    #plot(K.egg.vec, surv.mult.egg.vec, type="n", xlab="Carrying capacity of eggs(K)", ylab="Reduction in egg laying success")
    #lines(K.pred.egg.vec, pred.surv.egg.mult, lty=1, lwd=1, col="black")
    
    #plot(K.juv.vec, surv.mult.juv.vec, type="n", xlab="Proportion of the site carrying capacity (K)", ylab="Reduction in survival rate")
    #lines(K.pred.juv.vec, pred.surv.juv.mult, lty=1, lwd=1, col="black")
    
    ## invoke a density-feedback function on tadpole survival to year 1
    # density feedback survival multiplier for tadpoles hinges on the density of other tadpoles in the pond
    surv.mult.up <- 1.0
    surv.mult.upmid <- 0.58
    surv.mult.mid <- 0.19
    surv.mult.lo <- 0.10
    
    K.up <- 1
    K.upmid <- 0.83
    K.mid <- 0.45
    K.lo <- 0.3
    
    K.tad.vec <- c(K.up,K.upmid, K.mid,K.lo)
    surv.mult.tad.vec <- rev(c(surv.mult.up, surv.mult.upmid, surv.mult.mid, surv.mult.lo))
    #plot(K.tad.vec, surv.mult.tad.vec, pch=19)
    
    # Bleasdale
    # y = (a + bx)^(-1/c)
    DD.dat <- data.frame(K.tad.vec, surv.mult.tad.vec)
    param.init <- c(-2.41e-01, 1.54, 1.17)
    fit.expd.tad <- nls(surv.mult.tad.vec ~ (a + (b*K.tad.vec))^(-1/c), 
                        data = DD.dat,
                        algorithm = "port",
                        start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                        trace = TRUE,      
                        nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
    #plot(K.tad.vec, surv.mult.tad.vec, pch=19, xlab="K", ylab="reduction Tadpole survival to 1 yr")
    K.pred.tad.vec <- seq(K.lo,1,0.01)
    pred.surv.tad.mult <- (as.numeric(coef(fit.expd.tad)[1]) + (K.pred.tad.vec * as.numeric(coef(fit.expd.tad)[2])))^(-1/as.numeric(coef(fit.expd.tad)[3]))
    #lines(K.pred.tad.vec, pred.surv.tad.mult, lty=2, lwd=1, col="red")
    
    ## invoke a density-feedback function on Juvenile survival from year 1 to year 2
    # density feedback survival multiplier for juveniles hinges on the density of themselves (but more strongly than adults)
    surv.mult.up <- 1.0
    surv.mult.upmid <- 0.58
    surv.mult.mid <- 0.19
    surv.mult.lo <- 0.10
    
    K.up <- 1
    K.upmid <- 0.83
    K.mid <- 0.45
    K.lo <- 0.3
    
    K.juv.vec <- c(K.up,K.upmid, K.mid,K.lo)
    surv.mult.juv.vec <- rev(c(surv.mult.up, surv.mult.upmid, surv.mult.mid, surv.mult.lo))
    #plot(K.juv.vec, surv.mult.juv.vec, pch=19)
    
    # Bleasdale
    # y = (a + bx)^(-1/c)
    DD.dat <- data.frame(K.juv.vec, surv.mult.juv.vec)
    param.init <- c(-2.41e-01, 1.54, 1.17)
    fit.expd.juv <- nls(surv.mult.juv.vec ~ (a + (b*K.juv.vec))^(-1/c), 
                        data = DD.dat,
                        algorithm = "port",
                        start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                        trace = TRUE,      
                        nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
    
    K.pred.juv.vec <- seq(K.lo,1,0.01)
    pred.surv.juv.mult <- (as.numeric(coef(fit.expd.juv)[1]) + (K.pred.juv.vec * as.numeric(coef(fit.expd.juv)[2])))^(-1/as.numeric(coef(fit.expd.juv)[3]))
    
    
    ## The Markov chain PrChange values
    
    transitionLikelihood <- read.csv("C:/Workspace/AspatialtransitionLikelihood.csv", header = TRUE, sep = ",", dec = ".")
    transitionLikelihood <- transitionLikelihood[,2:ncol(transitionLikelihood)]
    rownames(transitionLikelihood) <- c("W2D","D2W")
    
    
    
    ###############################################################################################################
    ## And NOW for the actual model!
    ###############################################################################################################
    
    ## sillFlow tracks the Likelihood of transition under sill height 7 m - 11 m and regulated vs natural as 1:18  as found in transitionLikelihood 
    ## For the Latin Hypercube we will use only the 7.5 m sill height under regulated flow conditions (4)
    for (sillFlow in 4:4) {
     
      
      ## The Outermost Loop: normally we would cycle through population sizes (1,2,3,4/S,M,L,XL)
      ## For the Latin Hypercube we will use only the Large size wetland (3)
      for (size in 3:3) {
        
        ## The Second Loop: iterate the process iter times  *defined at the top of the code
        for (e in 1:iter) {
          #reset the population matrix except for the starting populations (n.mat[,1,])
          n.mat[1:maxAge+1,2:ncol(n.mat),1:4] <- 0
          
          # track the string of 0101010 for wet and dry  
          WDTracker <- 0
          # track the generation we reached outside of the i loop  
          Gen <- 1   
          
          ## The Third Loop: run the current projection set up for 80 years then update to outPut and popChanges
          for (i in 1:t) {
            #iterate the outPut array 
            if (i==1) { 
            outPut[1,e] <- outPut[1,e] + 1
            Wetness <- 1
            }
            
            
            # if there are no frogs alive break from this loop
            if (sum(n.mat[,i,size]) == 0) {  
             
              break  }
            
            # break on 80
            if (i==80) { 
              if (sum(n.mat[1:6,i,3]) > 0) {
                outPut[2,e] <- outPut[2,e] + 1
              }
              break 
              }
            
            # resample the duration of egg to tadpole
            #tomet.dur.iter <- round(sum(c(runif(1, hatch.dur[1], hatch.dur[2]), runif(1, tadpole.dur[1], tadpole.dur[2]))), 0)
            tomet.dur.iter <- round((rnorm(1, hatch.dur, hatch.dur.sd.prop*hatch.dur) + rnorm(1, tadpole.dur, tadpole.dur.sd.prop*tadpole.dur)), 0)
            
            toad.dur.iter <- 365 - tomet.dur.iter
            
            ##Calculate our survivals
            #Likelihood of egg hatching and tadpole surviving to metamorphosis
            tomet.s.iter <- (rbeta(1, tp.s.alpha, tp.s.beta)) * (runif(1, min=0.933, max=1))
            toad.daily.s.iter <- nthroot(rbeta(1,ad.s.yr.alpha, ad.s.yr.beta) , 365)
            toad.s.season.iter <- toad.daily.s.iter ^ toad.dur.iter
            toad.s.iter <- tomet.s.iter * toad.s.season.iter                     
            ad.s.iter <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)
            
            #Create the survival vector (popmat) for the year (density dependence not considered)
           
          
             ad.s.vec.iter <- rep(NA,(maxAge-1))
             for (s in 1:(maxAge-1)) {
               ad.s.vec.iter[s] <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)  }
             surv.iter <- c(toad.s.iter, ad.s.vec.iter,0)
             
            
            min.clutch.size.lr <- as.numeric(quantile(rnorm(n=1000, mean=clutch.size.lr, sd=clutch.size.lr.sd.prop*clutch.size.lr), probs=0.025, na.rm=T))
            max.clutch.size.lr <- as.numeric(quantile(rnorm(n=1000, mean=clutch.size.lr, sd=clutch.size.lr.sd.prop*clutch.size.lr), probs=0.975, na.rm=T))
             
            fert.iter <- round((runif((maxAge+1), min=min.clutch.size.lr, max=max.clutch.size.lr)) * prop.breeding, 0)
            popmat[1,] <- fert.iter * sex.ratio 
            diag(popmat[2:(maxAge+2), ]) <- surv.iter[-stages]
            
            #note these are placeholder values for popmat.current, they are recalculated below not necesary but retained for safety
            popmat.current <- popmat
            
            # Implement density dependence effects for each of four densities and feed into the popmat.all[,,]
            #note density effect on egg laying is applied after matrix multiplication 
            #  density feedback for tadpoles to first year is strong and driven by the number of tadpoles in the cohort
            s.mult.iter.tad <- 1
            s.mult.iter.juv <- 1
            s.mult.iter.ad <- 1
            
            # instil density dependence for tadpoles growing into year 1 adults
            K.rel.tad <- (n.mat[1,i,size]/K.rel.tad.dem)
            if (!is.nan(K.rel.tad)) {
              if (K.rel.tad > 2.1)  { K.rel.tad <- 2.1 }
              if (K.rel.tad <= 2.1) {
                s.mult.iter.tad <- (as.numeric(coef(fit.expd.tad)[1]) + (K.rel.tad * as.numeric(coef(fit.expd.tad)[2])))^(-1/as.numeric(coef(fit.expd.tad)[3]))
                popmat.current[2,1] <- popmat[2,1] * s.mult.iter.tad   } 
            } 
            
            # instil density dependence for juveniles (1 - 2 years) is driven by the  number of yr 1 present/competing per 
            K.rel.juv <- (n.mat[2,i,size]/K.rel.juv.dem) 
            if (!is.nan(K.rel.juv)) {
              
              if (K.rel.juv > 2.1)  { K.rel.juv <- 2.1 }
              
              if (K.rel.juv <= 2.1) {
                s.mult.iter.juv <- (as.numeric(coef(fit.expd.juv)[1]) + (K.rel.juv * as.numeric(coef(fit.expd.juv)[2])))^(-1/as.numeric(coef(fit.expd.juv)[3]))
                popmat.current[3,2] <- popmat[3,2] * s.mult.iter.juv    } 
            } 
            
            # instill density dependence for adults  driven by the  number of yr 1s emerging from Berven 2009
            # NOTE K.rel.adult.dem.vec[size] values are set so high as to be irrelevant
            # This is done to offset the very low survival values used
            ## NOTE CHNGED ABOVE FROM >0.5 TO <0.999 CONSIDER THIS IF REJIGGNG ADULTS 
            K.rel.adult <- (n.mat[2,i,size]/K.rel.adult.dem.vec[size]) 
            if (!is.nan(K.rel.adult)) {
              if (K.rel.adult > 0.999)  { K.rel.adult <- 0.999 }
              if (K.rel.adult > 0.51) {
                for (adgens in 4:5) {
                  s.mult.iter.ad <- (as.numeric(coef(fit.expd.adult)[1]) + (K.rel.adult * as.numeric(coef(fit.expd.adult)[2])))^(-1/as.numeric(coef(fit.expd.adult)[3]))
                  popmat.current[adgens,(adgens-1),size] <- (popmat[adgens,(adgens - 1)] * s.mult.iter.ad)
                } } }   
            
            # set popmat.fail for use during dry years
            popmat.fail <- popmat.current
            popmat.fail[1,] <- 0  
            
            
            ## The matrix multiplication step
            if (Wetness == 1) {
              swap <- rbinom(1, 1, transitionLikelihood[1,sillFlow]) 
              if (swap == 1) {   Wetness <- 0 }
            } else if (Wetness == 0) {
              swap <- rbinom(1, 1, transitionLikelihood[2,sillFlow]) 
              if (swap == 1) {   Wetness <- 1 }
            } else {
              print("Wetness is neither 0 nor 1")
              break
            }
            
            if (Wetness == 1) {
            hiya  <- popmat.current %*% n.mat[,i,size]
               n.mat[,i+1,size] <- hiya[1:(maxAge+1)]
              WDTracker <- paste(c(WDTracker,"1"),collapse = "")
            } else if (Wetness == 0) {
              hiya  <- popmat.fail %*% n.mat[,i,size]
               n.mat[,i+1,size] <- hiya[1:(maxAge+1)]
              WDTracker <- paste(c(WDTracker,"0"),collapse = "")
            } else {
              stop("hmm winFail is neither 1 or 0") 
            }
            
            Gen <- Gen + 1
            
            for (ii in 2:maxAge+1) {       n.mat[ii,i+1,size] <- round(n.mat[ii,i+1,size], digits = 0) }   
            
            #  density feedback for eggs
            #note evidence says there is no density dependence on egg laying in female L. aurea  But there must be some inhibition re maximum volume if nothing else
            # I now use integration. i.e. early in the curve females will lay with 100% success. As it tends towards the pond limit 
            # successive females lay  with diminishing success 
            # above the egg limit for the pond, laying is possible but with a huge inhibition
            K.rel.egg <- (n.mat[ 1, i+1, size]/K.rel.egg.dem.vec[size]) 
            if (!is.nan(K.rel.egg) && (K.rel.egg > 0)) {
              if (K.rel.egg <= 0.99) {  
                intArea <-  integrate(eggfunc,lower=0,upper=K.rel.egg)
                area <- as.numeric(intArea[1])/K.rel.egg
                if (area >= 1) { area <- 1 }
                n.mat[ 1, i+1, size] <- (round(n.mat[ 1, i+1, size] * area))
              } else if (K.rel.egg > 0.99) {
                # if  > than the limit for the pond then calculate what happens to the first 99% of the pond limit (egg1) then apply the 99th%ile inhibition on the remaining eggs(remain)
                egg1 <- (K.rel.egg.dem.vec[size] * area1)
                remain <- (n.mat[ 1, i+1, size] - K.rel.egg.dem.vec[size]) 
                n.mat[ 1, i+1, size] <- (round(egg1 + (s.mult.egg.iter * remain)))  
              } else { 
                stop("Crashed at line 1894: the egg conversion value K.rel.egg is misbehaving")
              }  
            }
            
            # theoretically unnecesary but included for stability
            for (clean in 1:maxAge+1) {
              if (is.nan(n.mat[clean,i+1,size])) {
                n.mat[clean,i+1,size] <- 0 
              }  } 
            
            
            if (i > 80){
              stop("failed to break at the 80th year")
            }
    
            ##  Last line of the generation loop (86 years)
          }
          
          # save
          
          input$PrExt <- 1 - (sum(outPut[2,], na.rm=T) / iter)
          save.nm <- paste0('res',sprintf("%09.0f", rowNum))
          assign(save.nm, input)
          save(list=save.nm,file=paste(dir.nm,save.nm,sep='/'))
          
          #print(e)
          #if (e %% itdiv==0) print(e) 
          
          ##  Last line of the 'iter' Loop  
        }
        ### Last line of the Size Loop, 
      }
      ## Last line of the sillFlow loop 
}

    print("*******************")
    print(d) 
    print("*******************")
} # end lh loop




## parameter ranges
ranges <- list()

ranges$clutch.size.lr <- c(1000,6000)
ranges$hatch.dur <- c(0.5, 14)
ranges$hatch.pr <- c(0.5, 0.98)
ranges$tadpole.dur <- c(50, 90)
ranges$tadpole.mn <- c(0.03, 0.6)
ranges$ad.s.yr.mn <- c(0.03, 0.35)
ranges$maxAge <- c(5, 10)
ranges$K.rel.tad.dem <- c(500, 5000)
ranges$K.rel.juv.dem <- c(200, 5000)

## create hypercube
nSamples <- 1000
lh <- data.frame(randomLHS(n=nSamples, k=length(ranges)))
names(lh) <- names(ranges)

## convert parameters to required scale
for (j in 1:ncol(lh)) {
  par <- names(lh)[j]
  lh[,par] <- qunif(lh[,j], min=ranges[[par]][1], max=ranges[[par]][2]) ## continuous
}

## number of iterations for each parameter set
lh$iter <- 1

## folder for saving the results of each row
## we could just store in memory, but then if something breaks we lose these data
dir.nm <- 'GSAsbf'
dir.create(dir.nm)

## rerun the simulation function in parallel using the LH values
res <- foreach(rowNum=1:nrow(lh),.verbose=T) %do% {sbf_sim(input=lh[rowNum,],dir.nm=dir.nm,rowNum=rowNum)}

# IF RERUNNING FROM PREGENERATED FILES * CHANGE LAST FOLDER NAME ----
#dir.nm <- "c:/Workspace/GSAsbf_repasted 1k runs"

## retrieve results
res.nms <- list.files(dir.nm)
res.list <- lapply(res.nms, function(x) {load(paste(dir.nm,x,sep='/')) ; print(x) ; return(eval(as.name(x)))})
#res.list1 <- lapply(res.nms1, function(x) {load(paste('GSA50kSparatest5Dmax50pc1',x,sep='/')) ; print(x) ; return(eval(as.name(x)))})
dat <- data.frame(rbindlist(res.list))
#dat <- rbind(dat1,dat2,dat3,dat4)
head(dat)
dim(dat)[1]
tail(dat)
sum(is.na(dat$PrExt))







#########
## BRT ##
#########
# Take the outputs generated in the simulation loop above and fit a boosted regression tree to the data 
dat.nona <- data.frame(na.omit(dat[!is.infinite(rowSums(dat)),]))
dat.nona <- dat.nona[,-9]
dat.nona$cllPrExt <- cloglog(dat.nona$PrExt)
dim(dat.nona)[1]
head(dat.nona)

brt.fit <- gbm.step(dat.nona, gbm.x = attr(dat.nona, "names")[1:8], gbm.y = attr(dat.nona, "names")[9], family="laplace", max.trees=100000, 
                    tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.75, tree.complexity = 2)

summary(brt.fit)
dim(dat.nona)[1]
D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit)
gbm.plot.fits(brt.fit)

CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
CV.cor
CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
CV.cor.se
print(c(CV.cor, CV.cor.se))

eq.sp.points <- 100
RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=12)

## output average predictions
for (p in 1:8) {
  RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
  RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
}
RESP.val.dat <- as.data.frame(RESP.val)
colnames(RESP.val.dat) <- brt.fit$var.names
RESP.pred.dat <- as.data.frame(RESP.pred)
colnames(RESP.pred.dat) <- brt.fit$var.names
RESP.val.dat
RESP.pred.dat

## plot the highest-influence variables only
plot(RESP.val.dat[,6], RESP.pred.dat[,6], type="l", xlab="mean adult annual survival", ylab="Pr(ext)")
plot(RESP.val.dat[,5], RESP.pred.dat[,5], type="l", xlab="tadpole survival to metamorphosis", ylab="Pr(ext)")

# save the plots and the image
setwd("C:/workspace/math0286/R/")
write.table(RESP.val.dat,file="BRT.val.GSAsbf.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.pred.dat,file="BRT.pred.GSAsbf.csv",sep=",", row.names = T, col.names = T)
save.image(paste("GSAsbf",".RData",sep=""))
