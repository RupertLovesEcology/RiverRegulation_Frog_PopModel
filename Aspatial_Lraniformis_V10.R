
##############################################################################################################
## V10.0 
## 4 populations of varying sizes all density feedback
## probability of wet or dry for natural and regulated flow scenarios over a range of sill heights
## at the end of each run, revise the results to find the desired sequence (e.g. WWDDD) and record the extinction/survival 
## store the % reduction over each of the sequences
##############################################################################################################
#Nesting is 1000 runs then wet between the dry then duration of the dry then  population size then T years  
## This version includes no pre-calculation of arrays and uses very large arrays for data storage
## If memory is limiting consider splitting the large storage arrays
## if runtime is unacceptable consider parrallelisation and/or precalculate some arrays and/or split the larger storage arrays 

##################################################
## southern bell frog population projection model
#################################################

## CJA Bradshaw & Rupert Mathwin
## September 2022

## Remove everything
rm(list = ls())

##Changed source file location below
source("C:/workspace/math0286/R/win-library/3.6/matrixOperators.r")
setwd("C:/workspace/math0286/R")
.libPaths("C:/workspace/math0286/R/win-library/4")

## libraries
library(DataCombine)
library(DescTools)
library(ggplot2)
library(plotly)
library(pracma)
library(grid, lib.loc = "C:/Program Files/R/R-4.0.2/library")
library(gridExtra)
library(cowplot)
library(magick)

# beta distribution shape parameter estimator function
##Generates an alpha and a beta value to inform the beta distribution 
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta)) }


WDTracker <- ""

## create the final outPut array
## Third Dimension is SillFlow
## 4th Dimension: 1 - count of occurrence, 2- mean change % adults, 3 - StDev % change adults, 4 - number of times extinct (adults),
# 5 - mean change % all (includes eggs), 6 - StDev change % all (includes eggs), 7 - number of times extinct (all), 
# 8 - % times ended in extinction (adult), 9 - % times ended in extinction (adult)
# The remaining slots (10 - 20) are placeholders for additional datapoints
# Output is preserved and data are passed to Output2 for manipulation
outPut2 <- outPut <- array(data = 0, dim = 64 * 80 * 18)
dim(outPut) <- c(8,8,18,4,20)
rownames(outPut) <- c("DDD","DDW","DWD","DWW","WDD","WDW","WWD","WWW")
colnames(outPut) <- c("DDD","DDW","DWD","DWW","WDD","WDW","WWD","WWW")

# consecDry tracks the extinctions and population decrease by consecutive dry years
# consecDry[size,DD - DDDDDD, ->] 
## Third Dimension is SillFlow
## 4rd dimension: 1- count of all occurrences, 2- count of extinctions in adults, 3 - count of extinctions in all, 4 -% extinction adult, 
# 5 - % extinction all, 6 - % decrease in adults, 7 - % decrease in all, 8 - 9 placeholders, 10 we keep track of the runs which we couldn't use for adult extinction
consecDry <- array(data = 0, dim = 200 * 18)
dim(consecDry) <- c(4,5,18,10)
rownames(consecDry) <- c("Small", "Medium", "Large", "Extra Large")
colnames(consecDry) <- c("DD", "DDD", "DDDD", "DDDDD", "DDDDDD")

# this hold the changes for each consecutive dry duration for averaging at the end
#durationHolder[size,DryDuration DD - DDDDDD, 1=adults only/2=all ages,sillFlow,50000]
durationHolder <- array(data = NA, dim = 2000000 * 18)
dim(durationHolder) <- c(4,5,2,18,50000)

###  STORAGE FOR ERRORCHECKING (population size - S,M,L,XL)
runCounter <- c(0,0,0,0)

## Store the individual population changes for each wet/Dry pattern in popChanges
## Note third dimension is 1 - Starting Adult, 2 - finishing adult, 3 - change in adult, 
# 4 - if adults was > 0 at year 1 and is ==0 at completion then store a 1(extinction this run), 5 - Starting population all life stages
# 6 - Final population all life stages, 7 - Change in all ages (in %), 8 - if pop goes extinct 1 otherwise 0
popChanges <- array(data = NA, dim = 512000000)
dim(popChanges) <- c(8,8,10,800000)

## iterate projection 
iter <- 10000
itdiv <- iter/10

## set time limit for projection in 1-yr increments
t <- 86
longev <- 5
age.vec <- seq(0,longev,1)
lage <- length(age.vec)
sex.ratio <- 0.5
stages <- lage

## set population storage matrices n.mat
n.mat <- array(data = 0, dim = (stages * (t + 1) * 4))
dim(n.mat) <- c(stages,(t + 1),4)
popmat <- matrix(0,nrow=stages,ncol=stages)
colnames(popmat) <- age.vec[1:stages]
rownames(popmat) <- c("Fecundity","Survival to 1 year","Survival to 2 year","Survival to 3 year","Survival to 4 year","Survival to 5 year")


## fertility data 
clutch.size.lr <- c(1885,3893) # L. raniformis
prop.breeding <- c(0,rep(1,5))
fert.mn <- mean(clutch.size.lr)*prop.breeding
fert.mn

#duration data (eggs and tadpoles)
hatch.dur <- c(2,4) # 2-4 days, L. raniformis
tadpole.dur <- c(70,80) # duration (days) L. raniformis 23 deg

## survival data
hatch.pr <- c(0.933, 1) 

##the below uses the more complete survival TO METAMORPHOSIS figures from Bull (C.signifera)
tadpole.mn.1 <- mean(c(.15,.26))
tadpole.sd.1 <- ((.26-tadpole.mn.1)+(tadpole.mn.1-.15))/2/1.96
tadpole.mn.2 <- mean(c(.07,.56))
tadpole.sd.2 <- ((.56-tadpole.mn.1)+(tadpole.mn.1-.07))/2/1.96
tadpole.mn <- mean(c(tadpole.mn.1, tadpole.mn.2))
tadpole.sd <- sqrt(tadpole.sd.1^2 + tadpole.sd.2^2)
tp.s.alpha <- estBetaParams(tadpole.mn, tadpole.sd/10)$alpha
tp.s.beta <- estBetaParams(tadpole.mn, tadpole.sd/10)$beta

# ADULT ANNUAL SURVIVAL VALUES FROM PICKETT L. AUREA
ad.s.yr.mn <- 0.2172 # Litoria aurea (Pickett et al.)
ad.s.yr.sd <- 0.087 # Litoria aurea (Pickett et al.)

##FROM PICKETT
ad.s.yr.alpha <- estBetaParams(ad.s.yr.mn, ad.s.yr.sd/10)$alpha
ad.s.yr.beta <- estBetaParams(ad.s.yr.mn, ad.s.yr.sd/10)$beta

tomet.dur.iter <- round(sum(c(runif(1, hatch.dur[1], hatch.dur[2]), runif(1, tadpole.dur[1], tadpole.dur[2]))), 0)
toad.dur.iter <- 365 - tomet.dur.iter

##Calculate some valuies for the population matrix
# probability of egg hatching and tadpole surviving to metamorphosis
tomet.s.iter <- (rbeta(1, tp.s.alpha, tp.s.beta)) * (runif(1, min=0.933, max=1))

#sample a daily probability of survival
toad.daily.s.iter <- nthroot(rbeta(1,ad.s.yr.alpha, ad.s.yr.beta) , 365)

toad.s.season.iter <- toad.daily.s.iter ^ toad.dur.iter
toad.s.iter <- tomet.s.iter * toad.s.season.iter                     

# TO RESAMPLE annual survival of an adult
#ad.s.iter <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)

#Create the survival vector, to adult survival then 4 adult survivals
ad.s.vec.iter <- rep(NA,5)
for (s in 1:5) {
  ad.s.vec.iter[s] <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)  }
surv.iter <- c(toad.s.iter, ad.s.vec.iter)



# We hereafter consider 4 population sizes, each represented with 4 density feedback values (one for each of egg, tad, juve, adult)
# Small, medium, large, very large.
# Assign the number of spawning masses that each site holds and convert it to the number of female egges (uses mean clutch size)
K.rel.egg.dem.vec <- c(50, 150, 500, 1000)
K.rel.egg.dem.vec <- K.rel.egg.dem.vec * 1444

# tadpoles surviving to 1 year old will inhibit themselves NOTE THIS EQUATION ONLY CONSIDERS FEMALES SO THE ACTUAL NUMBER WILL BE DOUBLE  
K.rel.tad.dem.vec <- c(100, 300, 800, 1000)

# number of first years present will influence survival from 1 to 2 years old but more strongly than for older age brackets         
K.rel.juv.dem.vec <- c(200, 600, 1800, 3000)

# number of first years present will influence adult survival NOTE THIS EQUATION ONLY CONSIDERS FEMALES SO ACTUAL NUMBER WILL BE DOUBLE         
K.rel.adult.dem.vec <- c(1000000,1000000,1000000,1000000)

initSpawn <- c(10,30,85,150)

#set founding fem pop sizes 
init.fem.pop <- c(20,65,150,500) 

init.vec <- rep(NA,24)
dim(init.vec) <- c(6,4)

##Populate the Matix (popmat) and create a failure matrix(popmat.fail) for years with no breeding
diag(popmat[2:(stages), ]) <- surv.iter[-stages]
popmat[stages,stages] <- 0 # surv.mn[stages] 
popmat[1,] <- fert.mn * sex.ratio
popmat.orig <- popmat ## save original matrix as popmat.orig

n.sums.mat <- array(data = 0, dim = (iter * (t + 1) * 4)) 
dim(n.sums.mat) <-c(iter, (t+1), 4)

#Create init.vec note 2nd dimension is 1 = S, 2 = M, 3 = L, 4 = XL
totalSurv <- sum(popmat[2:6,1:6])
for (size in 1:4) {
  for (x in 1:6) { init.vec[x,size] <- 0 }
  init.vec[1,size] <- round(initSpawn[size] * (runif(1, min=clutch.size.lr[1], max=clutch.size.lr[2])))
  for (dd in 2:6) { init.vec[dd,size] <- round((sum(popmat[2:6,(dd - 1)])/totalSurv) * (init.fem.pop[size]))    } }



for (size in 1:4) { n.mat[,1,size] <- init.vec[,size] }

# Create population matrices for Leslie modelling
popmat.fail <- popmat.current <- popmat
popmat.fail[1,] <- 0 

# Create three-dimensional storage matrices & vectors for graphing (note third dimension S = 1, M = 2, L = 3 and XL = 4)
n.sums.mat <- array(data = 0, dim = (iter * (t + 1) * 4))
dim(n.sums.mat) <-c(iter, (t+1), 4)
r.mat <- array(data = 0, dim = (iter * (t + 1) * 4))
dim(r.mat) <-c(iter, (t + 1), 4)

# use temporary variables K.egg.vec, surv.mult.egg.vec, pred.surv.egg.mult, and K.pred.egg.vec to calculate fit.exp.egg
# fit.exp.egg is used to invoke a density-feedback function on egg laying 
# I  use integration. i.e. early in the curve females will lay with 100% success. As it tends towards the pond limit 
# successive females lay  with diminishing success 
# above the egg limit for the pond, laying is possible but with a strong inhibition
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
plot(K.egg.vec, surv.mult.egg.vec, pch=19)
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
lines(K.pred.egg.vec, pred.surv.egg.mult, lty=2, lwd=1, col="red")

eggfunc <- function(x) (1.01 - ((x ^ (9/x)) * (x ^ 2)))
intArea <- integrate(eggfunc,lower=0,upper=0.99)
area1 <- (as.numeric(intArea[1])/0.99)
s.mult.egg.iter <- (as.numeric(coef(fit.expd.egg)[1]) - (0.999 ^ (as.numeric(coef(fit.expd.egg)[2])/0.999) * (0.999 ^ (coef(fit.expd.egg)[3]))))

plot(K.egg.vec, surv.mult.egg.vec, type="n", xlab="Carrying capacity of eggs(K)", ylab="Reduction in egg laying success")
lines(K.pred.egg.vec, pred.surv.egg.mult, lty=1, lwd=1, col="black")

plot(K.juv.vec, surv.mult.juv.vec, type="n", xlab="Proportion of the site carrying capacity (K)", ylab="Reduction in survival rate")
lines(K.pred.juv.vec, pred.surv.juv.mult, lty=1, lwd=1, col="black")

# use temporary variables K.tad.vec, surv.mult.tad.vec, pred.surv.tad.mult, and K.pred.tad.vec to calculate fit.exp.tad
## fit.exp.tad is used to invoke a density-feedback function on tadpole survival in the first year after hatch
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
plot(K.tad.vec, surv.mult.tad.vec, pch=19)

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
plot(K.tad.vec, surv.mult.tad.vec, pch=19, xlab="K", ylab="reduction Tadpole survival to 1 yr")
K.pred.tad.vec <- seq(K.lo,1,0.01)
pred.surv.tad.mult <- (as.numeric(coef(fit.expd.tad)[1]) + (K.pred.tad.vec * as.numeric(coef(fit.expd.tad)[2])))^(-1/as.numeric(coef(fit.expd.tad)[3]))
lines(K.pred.tad.vec, pred.surv.tad.mult, lty=2, lwd=1, col="red")

# use temporary variables K.juv.vec, surv.mult.juv.vec, pred.surv.juv.mult, and K.pred.juv.vec to calculate fit.exp.juv
## fit.exp.juv is used to invoke a density-feedback function on Juvenile survival from year 1 to year 2
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
plot(K.juv.vec, surv.mult.juv.vec, pch=19)

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
plot(K.juv.vec, surv.mult.juv.vec, pch=19, xlab="Proportion of the site carrying capacity (K)", ylab="Reduction in survival rate")
plot(K.juv.vec, surv.mult.juv.vec, type="n", xlab="Proportion of the site carrying capacity (K)", ylab="Reduction in survival rate")
K.pred.juv.vec <- seq(K.lo,1,0.01)
pred.surv.juv.mult <- (as.numeric(coef(fit.expd.juv)[1]) + (K.pred.juv.vec * as.numeric(coef(fit.expd.juv)[2])))^(-1/as.numeric(coef(fit.expd.juv)[3]))
lines(K.pred.juv.vec, pred.surv.juv.mult, lty=2, lwd=1, col="red")
plot(K.juv.vec, surv.mult.juv.vec, type="n", xlab="Proportion of the site carrying capacity (K)", ylab="Reduction in survival rate")
lines(K.pred.juv.vec, pred.surv.juv.mult, lty=1, lwd=1, col="black")



######################################################
## Process the Wet/Dry transitions probabilities for each of 18 sill height/regulations combinations to inform the Markov-Chain 
## 7N, 7R, 7.5N, 7.5R - 11N, 11R
## Store the probability of transition 
######################################################

# upload the observed rating curve data to convert flow to river height in this reach (for the natural flow scenario)
RCurve <- read.csv("C:/Workspace/ratingCurve.csv", header = TRUE, sep = ",", dec = ".")
x <- RCurve$height
y <- RCurve$flow

# Find the polynome that fits this function (flow2Height.model)
flow2Height.model <- lm(y ~ x + I(x^2) + I(x^3))
# store the coefficients in f2H
f2H <- flow2Height.model$coefficients

# download natural flow scenario
NatFlowRate <- read.csv("C:/Workspace/NaturalFlow.csv", header = TRUE, sep = ",", dec = ".")

# convert natural flow  to river height
rF <- NatFlowRate$Natural
rH <- f2H[1] + (f2H[2]*rF) + (f2H[3]*rF^2) + (f2H[4]*rF^3)
NatFlowRate$Natural <- rH

# paste in regulated river scenario heights
regHeight <- read.csv("C:/Workspace/RegHeight.csv", header = TRUE, sep = ",", dec = ".")
regHeight2 <- matrix(data = regHeight, nrow = 32264, ncol = 1)

Heights <- data.frame(NatFlowRate$Day,NatFlowRate$Month,NatFlowRate$Year,regHeight,NatFlowRate$Natural)


# OR just import the csv of pre-calculated river heights under natural and regulated flow conditions by uncommenting below     d (^.^) b
# Heights <- read.csv("C:/Workspace/CompHeightL3.csv", header = TRUE, sep = ",", dec = ".")


CurrentYr <- 1925
NatTracker <- 0
InvTracker <- 0
OutPut <- c(0,0,0)
dim(OutPut) <- c(1,3)
colnames(OutPut) <- c("Year","Natural","Regulated")

# collect every natural and gauge height data for winter and spring into an array
for (i in 1:nrow(Heights)) {
  if (Heights[i,3] == CurrentYr && Heights[i,2] >= 6 && Heights[i,2] <= 11) {
    NatTracker <- append(NatTracker, Heights[i,5], after = length(NatTracker))
    InvTracker <- append(InvTracker, Heights[i,4], after = length(InvTracker))
  } 
  # order them numerically then store the 10th highest value each year for each of Nat and gauge
  if (as.numeric(Heights[i,3]) > as.numeric(CurrentYr)) {
    
    ## Do the 'find largest' step
    NatTracker <- sort(NatTracker)
    InvTracker <- sort(InvTracker)
    Place <- length(NatTracker) - 9
    Place2 <- length(InvTracker) - 9
    NewRow <- c(Heights[i,3],NatTracker[Place],InvTracker[Place2])
    OutPut <-  InsertRow(OutPut, NewRow, RowNum = NULL)
    CurrentYr <- Heights[i,3]
    NatTracker <- 0
    InvTracker <- 0
  }
  
} 
colnames(OutPut) <- c("Year","Natural","Regulated")
OutPut <- OutPut[-c(1),]
OutPut <- OutPut[,2:3] - 0.1

## Calculate Wet (1) and Dry (0) iteration for different sill heights

HeightsColNm <- c(7,7.5,8,8.5,9,9.5,10,10.5,11)
Scenarios <- c("Natural", "Regulated")
NameList <- array(data = NA, dim = 2 * length(HeightsColNm))
WetDry <- array(data = NA, dim = 2 * length(HeightsColNm) * 83)
dim(WetDry) <- c(83,2 * length(HeightsColNm))
Upper <- 2*length(HeightsColNm)
for (i in 1:Upper) {
  
  if (i %% 2 == 1) { 
    name <- paste(HeightsColNm[ceiling(i/2)], Scenarios[1], sep = " ")
    NameList[i] <- name
    for (j in 1:83) {
      if (OutPut[j,1] >= HeightsColNm[ceiling(i/2)]) {
        WetDry[j,i] <- 1
      } else {
        WetDry[j,i] <- 0
      }
    }
  }  else if (i %% 2 == 0) {
    name <- paste(HeightsColNm[ceiling(i/2)], Scenarios[2], sep = " ")
    NameList[i] <- name
    for (j in 1:83) {
      if (OutPut[j,2] >= HeightsColNm[ceiling(i/2)]) {
        WetDry[j,i] <- 1
      } else {
        WetDry[j,i] <- 0
      }
    }
  } 
} 

colnames(WetDry) <- NameList

## This produces an array WetDry. Data are arranged in columns with each column being a different treatment and each row is a year (from 1926 to 2009, n = 84) 
## re column names "8 Natural" = did a wetland 8 m above the sea fill up this year under a natural flow scenario. 1 = Yes and 0 = No
# write.csv(WetDry,"c:/Workspace/WetDry.csv") # uncomment to save a local copy

## Create an array to store the probability of transition from wet to dry/dry to wet
transitionLikelyhood <- array(data = 0, dim = ncol(WetDry) * 2)
dim(transitionLikelyhood) <- c(2,18)
rownames(transitionLikelyhood) <- c("W2D","D2W")
colnames(transitionLikelyhood) <- colnames(WetDry)

#Probability of transition for each of 2 states (18 columns)
for (r in 1:ncol(WetDry)) {
  #  reset state before starting the next column
  state <- WetDry[1,r]
  Wtrans <- Wsame <- Dtrans <- Dsame <- 0
  for (s in 2:nrow(WetDry)) {
    if (state == 0) {
      if (WetDry[s,r] == 0) { 
        Dsame <- Dsame + 1
      } else {
        Dtrans <- Dtrans + 1
        state <- 1
      }
      
    } else if (state == 1) {
      if (WetDry[s,r] == 1) { 
        Wsame <- Wsame + 1
      } else {
        Wtrans <- Wtrans + 1
        state <- 0
      }
    } else { 
      stop("While calculating likelyhood of transition I have achieved a state that is neither 1 or 0!  Ohmmmmmmmm")
    }
  } 
  #calculate, assign the probability of change, and reset
  transitionLikelyhood[1,r] <- Wtrans/(Wtrans + Wsame)
  transitionLikelyhood[2,r] <- Dtrans/(Dtrans + Dsame) 
}



###############################################################################################################
## And NOW for the actual model!
###############################################################################################################

## The Outermost Loop: sillFlow tracks the likelyhood of transition under sill height 7 m - 11 m and regulated vs natural as 1:18  as found in transitionLikelyhood 
for (sillFlow in 1:18) {
  cat("starting sillFlow ", sillFlow, " which is ", NameList[sillFlow],"\n")
  
  ## Next we cycle through population sizes (1,2,3,4/S,M,L,XL)
  for (size in 1:4) {
    
    ## The Second Loop: iterate the process times
    for (e in 1:iter) {
      cat("starting iteration ", e, "\n") 
      
      ## set value storage vectors Note t is the number of generations we are currently running (note third dimension S = 1, M = 2, L = 3 and XL = 4)
      # This array tracks the rate of population change. I did not use in the final outputs but I have retained it as a useful metric to record
      # However, it must be treated with caution because of the boom-and-bust nature of fertile frog species
      r.stoch <- array(data = 0, dim = ((t + 1) * 4))
      dim(r.stoch) <-c((t + 1), 4)

      #reset the population matrix except for the starting populations (n.mat[,1,])
      n.mat[1:6,2:ncol(n.mat),1:4] <- 0
      
      # track the string of 0101010 for wet and dry  
      WDTracker <- 0
      # track the generation we reached outside of the i loop  
      Gen <- 1   
      
      ## The Third Loop: run the current projection set up for 80 years then update to outPut and popChanges
      for (i in 1:t) {
        cat("starting generation ", i, "\n")
        #iterate the runCounter used for errorchecking
        if (i==1) { runCounter[size] <- (runCounter[size] + 1) 
        Wetness <- 1
        }
   
        # if there are no frogs alive break from this loop
        if (sum(n.mat[,i,size]) == 0) {  break  }
        
        # break on 80
        if (i==80) { break }
        
        # resample the duration of egg to tadpole
        tomet.dur.iter <- round(sum(c(runif(1, hatch.dur[1], hatch.dur[2]), runif(1, tadpole.dur[1], tadpole.dur[2]))), 0)
        toad.dur.iter <- 365 - tomet.dur.iter
        
        ##Calculate our survivals
        #likelyhood of egg hatching and tadpole surviving to metamorphosis
        tomet.s.iter <- (rbeta(1, tp.s.alpha, tp.s.beta)) * (runif(1, min=0.933, max=1))
        toad.daily.s.iter <- nthroot(rbeta(1,ad.s.yr.alpha, ad.s.yr.beta) , 365)
        toad.s.season.iter <- toad.daily.s.iter ^ toad.dur.iter
        toad.s.iter <- tomet.s.iter * toad.s.season.iter                     
        ad.s.iter <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)
        
        #Create the survival vector (popmat) for the year (density dependence not considered)
        ad.s.vec.iter <- rep(NA,5)
        for (s in 1:5) {
          ad.s.vec.iter[s] <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)  }
        surv.iter <- c(toad.s.iter, ad.s.vec.iter)
        
        fert.iter <- round((runif(stages, min=clutch.size.lr[1], max=clutch.size.lr[2])) * prop.breeding, 0)
        popmat[1,] <- fert.iter * sex.ratio 
        diag(popmat[2:(stages), ]) <- surv.iter[-stages]
        
        #note these are placeholder values for popmat.current, they are recalculated below so not strictly necessary but retained for safety
        popmat.current <- popmat
        
        # Implement density dependence effects 
        #note density effect on egg laying is applied after matrix multiplication 
        #  density feedback for tadpoles to first year is strong and driven by the number of tadpoles in the cohort
        s.mult.iter.tad <- 1
        s.mult.iter.juv <- 1
        s.mult.iter.ad <- 1
        
        # instil density dependence for tadpoles growing into year 1 adults
        K.rel.tad <- (n.mat[1,i,size]/K.rel.tad.dem.vec[size])
        if (!is.nan(K.rel.tad)) {
          if (K.rel.tad > 2.1)  { K.rel.tad <- 2.1 }
          if (K.rel.tad <= 2.1) {
            s.mult.iter.tad <- (as.numeric(coef(fit.expd.tad)[1]) + (K.rel.tad * as.numeric(coef(fit.expd.tad)[2])))^(-1/as.numeric(coef(fit.expd.tad)[3]))
            popmat.current[2,1] <- popmat[2,1] * s.mult.iter.tad   } 
        } 
        
        # instil density dependence for juveniles (1 - 2 years) is driven by the  number of yr 1 present/competing per Berven
        K.rel.juv <- (n.mat[2,i,size]/K.rel.juv.dem.vec[size]) 
        if (!is.nan(K.rel.juv)) {
          
          if (K.rel.juv > 2.1)  { K.rel.juv <- 2.1 }
          
          if (K.rel.juv <= 2.1) {
            s.mult.iter.juv <- (as.numeric(coef(fit.expd.juv)[1]) + (K.rel.juv * as.numeric(coef(fit.expd.juv)[2])))^(-1/as.numeric(coef(fit.expd.juv)[3]))
            popmat.current[3,2] <- popmat[3,2] * s.mult.iter.juv    } 
        } 
        
        # instill density dependence for adults  driven by the  number of yr 1s emerging from Berven 2009
        # This is done to offset the very low survival values used
        K.rel.adult <- (n.mat[2,i,size]/K.rel.adult.dem.vec[size]) 
        if (!is.nan(K.rel.adult)) {
          if (K.rel.adult > 0.999)  { K.rel.adult <- 0.999 }
          if (K.rel.adult > 0.51) {
            for (adgens in 4:6) {
              s.mult.iter.ad <- (as.numeric(coef(fit.expd.adult)[1]) + (K.rel.adult * as.numeric(coef(fit.expd.adult)[2])))^(-1/as.numeric(coef(fit.expd.adult)[3]))
              popmat.current[adgens,(adgens-1),size] <- (popmat[adgens,(adgens - 1)] * s.mult.iter.ad)
            } } }   
        
        # set popmat.fail for use during dry years
        popmat.fail <- popmat.current
        popmat.fail[1,] <- 0  
        
        
        ## The matrix multiplication step
        if (Wetness == 1) {
          swap <- rbinom(1, 1, transitionLikelyhood[1,sillFlow]) 
          if (swap == 1) {   Wetness <- 0 }
        } else if (Wetness == 0) {
          swap <- rbinom(1, 1, transitionLikelyhood[2,sillFlow]) 
          if (swap == 1) {   Wetness <- 1 }
        } else {
          print("Wetness is neither 0 nor 1")
          break
        }
        
        if (Wetness == 1) {
          n.mat[,i+1,size] <- popmat.current %*% n.mat[,i,size]
          WDTracker <- paste(c(WDTracker,"1"),collapse = "")
        } else if (Wetness == 0) {
          n.mat[,i+1,size] <- popmat.fail %*% n.mat[,i,size]
          WDTracker <- paste(c(WDTracker,"0"),collapse = "")
        } else {
          stop("hmm winFail is neither 1 or 0") 
        }
        
        Gen <- Gen + 1
        
        for (ii in 2:6) {       n.mat[ii,i+1,size] <- round(n.mat[ii,i+1,size], digits = 0) }   
        
        #  density feedback for eggs
        #note evidence says there is no density dependence on egg laying in female L. aurea  But there must be some inhibition re maximum volume if nothing else
        # I  use integration. i.e. early in the curve females will lay with 100% success. As it tends towards the pond limit 
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
        
        # theoretically unnecessary but included for stability
        for (clean in 1:6) {
          if (is.nan(n.mat[clean,i+1,size])) {
            n.mat[clean,i+1,size] <- 0 
          }  } 
        
        
        if (i > 80){
          stop("failed to break at the 80th year")
        }
        
        ## save r for this iteration' stochastic matrix
        r.running <- log(sum(n.mat[,i+1,size], na.rm=T) / sum(n.mat[,i,size], na.rm=T))
        r.stoch[i,size] <- ifelse(r.running == -Inf, NA, r.running)
        
        #update r.stoch (track instantaneous rate of population change)
        # note this must be treated cautiously due to the explosive boom-and-bust in this species/system 
                for (kkk in 1:(ncol(r.stoch))) {
          for (kk in 1:ncol(r.stoch)) {
            r.mat[e,kkk,kk] <- r.stoch[kkk,kk] }}
        
        # update n.sums.mat using only adult frogs   
        for (kkk in 1:(ncol(n.mat))) {
          n.sums.mat[e,kkk,size] <- as.vector(sum(n.mat[2:6,kkk,size], na.rm=T)) } 
        
        ##  Last line of the generation loop (86 years)
      }
      
      #scan through the binaries for this generational run (allowing a 10 gen run-in time)
      if (Gen >= 16) {
        
        for (j in 10:(nchar(WDTracker) - 5)) {
          #Input must be a string already otherwise 000000 just becomes 0 and NA
          Rrow <- BinToDec(as.numeric(substr(WDTracker, j, j + 2))) + 1
          Ccol <- BinToDec(as.numeric(substr(WDTracker, j + 3, j + 5))) + 1
          
          #Great Lets start with the popChanges Array
          ## z <- length(popChanges[sillFlow,size,Rrow,Ccol,1,][!is.na(popChanges[size,Rrow,Ccol,1,])]) + 1
          z <- length(popChanges[Rrow,Ccol,1,][!is.na(popChanges[Rrow,Ccol,1,])]) + 1
          
          
          popChanges[Rrow,Ccol,1,z] <- sum(n.mat[2:6,j - 1,size])
          print(z)
          
          popChanges[Rrow,Ccol,2,z] <- sum(n.mat[2:6,j + 5,size])
          popChanges[Rrow,Ccol,5,z] <- sum(n.mat[,j - 1,size])
          popChanges[Rrow,Ccol,6,z] <- sum(n.mat[,j + 5,size])
          if (sum(n.mat[,j + 5,size]) == 0) { 
            popChanges[Rrow,Ccol,8,z] <- 1 
          } else {
            popChanges[Rrow,Ccol,8,z] <- 0 
          }
          popChanges[Rrow,Ccol,7,z] <- (100 * popChanges[Rrow,Ccol,6,z]/popChanges[Rrow,Ccol,5,z]) - 100 
          # easy ones done. Next we will do the conditionals % changeAdult and extinctionAdult
          if (popChanges[Rrow,Ccol,2,z] == 0 & popChanges[Rrow,Ccol,1,z] > 0) {
            popChanges[Rrow,Ccol,4,z] <- 1 
          } else {
            popChanges[Rrow,Ccol,4,z] <- 0 
          } 
          if (popChanges[Rrow,Ccol,1,z] != 0) {
            popChanges[Rrow,Ccol,3,z] <- (100 * popChanges[Rrow,Ccol,2,z]/popChanges[Rrow,Ccol,1,z]) - 100 
          }
          
          #Next we will assign what we can to the outPut    
          # add one to the output count for this combination
          outPut[Rrow,Ccol,sillFlow,size,1] <- outPut[Rrow,Ccol,sillFlow,size,1] + 1
        }   
        ##
        # Next we will assign to consecDry which tracks the extinctions and population decrease by consecutive dry years
        # consecDry[size,DD - DDDDDD, ->] 3rd dimension: 1- count of all occurrences, 2- count of extinctions in adults, 3 - count of extinctions in all, 
        # THESE ARE CALCULATED LATER 4 -% extinction adult,  5 - % extinction all, 6 - % decrease in adults, 7 - % decrease in all, 8 - 9 placeholders,
        # 10 we keep track of the runs which we couldn't use for adult extinction  
        # Also uses the holding array durationHolder to keep the individual changes for averaging later 
        J <- 10
        # finds the beginning (J) and end (end) of each DDD strings  
        while  (J <= nchar(WDTracker)) {
          end <- J 
          while (BinToDec(as.numeric(substr(WDTracker, J,end + 1))) == 0) {
            if (end == nchar(WDTracker)) { break }
            end <- end + 1
          }
          if (end + 1 - J > 1) {
            dur <- end - J
            
            consecDry[size,dur,sillFlow,1] <- consecDry[size,dur,sillFlow,1] + 1
            
            if  (sum(n.mat[2:6,J - 1,size]) == 0) {
              ## Sometimes we can't use this run for adult extinction because the starting pop is 0 we will keep count of this in 10
              consecDry[size,dur,sillFlow,10] <- consecDry[size,dur,sillFlow,10] + 1  
            } else if (sum(n.mat[2:6,end,size]) == 0  & (sum(n.mat[2:6,J,size]) > 0)) {
              consecDry[size,dur,sillFlow,2] <- consecDry[size,dur,sillFlow,2] + 1 
            }
            
            if (sum(n.mat[,end,size]) == 0) { 
              consecDry[size,dur,sillFlow,3] <- consecDry[size,dur,sillFlow,3] + 1 }
            ## NOTE FOR ERROR CHECKING THESE SHOULD BE THE SAME
            # THATS ALL WE CAN DO FOR NOW WITH consecDry  Lets track this generations dry duration changes into durationHolder
            #durationHolder[size,DryDuration DD - DDDDDD, 1=adults only 2=all ages,50000]
            Z <- length(durationHolder[size,dur,1,sillFlow,][!is.na(durationHolder[size,dur,1,sillFlow,])]) + 1
            if (sum(n.mat[2:6,J - 1,size]) > 0) {
              durationHolder[size,dur,1,sillFlow,Z] <- (100 * sum(n.mat[2:6,end,size])/sum(n.mat[2:6,J - 1,size])) - 100
            }
            ZZ <- length(durationHolder[size,dur,2,sillFlow,][!is.na(durationHolder[size,dur,2,sillFlow,])]) + 1
            durationHolder[size,dur,2,sillFlow,ZZ] <- (100 * sum(n.mat[,end,size])/sum(n.mat[,J - 1,size])) - 100 
          }
          J <- end + 1
        }
        ##   last line of the assignment to outputs loop    
      }  
      
     #comment this out to run faster (note this .R completes almost all calculations within these loops.)
      # it takes a long time. Pre-processing arrays and matrices will make this much faster. 
      # Can also rationalise the multidimensional arrays (outPut) if memory issues arise
      print(e)
      
      ##  Last line of the 'iter' Loop  
    }
    
    ## all done, with that combination of size and sillFlow
    ##  let's assign the remaining fields to outPut and then reset the popChanges array
    
    for (Rrow in 1:8) {
      for (Ccol in 1:8) {
        
        z <- length(popChanges[Rrow,Ccol,1,][!is.na(popChanges[Rrow,Ccol,1,])])  
        outPut[Rrow,Ccol,sillFlow,size,2] <- round(mean(popChanges[Rrow,Ccol,3,1:z], na.rm = TRUE), digits = 2)
        outPut[Rrow,Ccol,sillFlow,size,3] <- round(sd(popChanges[Rrow,Ccol,3,1:z], na.rm = TRUE), digits = 2)
        outPut[Rrow,Ccol,sillFlow,size,4] <- round(sum(popChanges[Rrow,Ccol,4,1:z], na.rm = TRUE), digits = 2)
        outPut[Rrow,Ccol,sillFlow,size,5] <- round(mean(popChanges[Rrow,Ccol,7,1:z], na.rm = TRUE), digits = 2)
        outPut[Rrow,Ccol,sillFlow,size,6] <- round(sd(popChanges[Rrow,Ccol,7,1:z], na.rm = TRUE), digits = 2)
        outPut[Rrow,Ccol,sillFlow,size,7] <- round(sum(popChanges[Rrow,Ccol,8,1:z], na.rm = TRUE), digits = 2)
        if (outPut[Rrow,Ccol,sillFlow,size,1] > 0) {
          outPut[Rrow,Ccol,sillFlow,size,8] <- round((outPut[Rrow,Ccol,sillFlow,size,4] * 100)/outPut[Rrow,Ccol,sillFlow,size,1], digits = 2)
          outPut[Rrow,Ccol,sillFlow,size,9] <- round((outPut[Rrow,Ccol,sillFlow,size,7] * 100)/outPut[Rrow,Ccol,sillFlow,size,1], digits = 2)
        }
      }
    }
    popChanges[] <- NA
    
    
    ## Last line of the Size Loop, iterates to the next sized population
  }
  ## Last line of the sillFlow loop iterates to the next sill / flow combination population
}

# all done lets assign the remaining results to consecDry[4,5,4 - 7]
# 4 % ext adult, 5- %ext all, 6- % change adult, 7- % change all
for (slflw in 1:18) {
  for (size in 1:4) {
    for (DUR in 1:5) {
      consecDry[size,DUR,slflw,4] <- 100 * consecDry[size,DUR,slflw,2] /(consecDry[size,DUR,slflw,1] - consecDry[size,DUR,slflw,10])
      consecDry[size,DUR,slflw,5] <- 100 * consecDry[size,DUR,slflw,3] /consecDry[size,DUR,slflw,1]
      Z <- length(durationHolder[size,DUR,1,slflw,][!is.na(durationHolder[size,DUR,1,slflw,])]) + 1
      consecDry[size,DUR,slflw,6] <- mean(durationHolder[size,DUR,1,slflw,1:Z], na.rm = TRUE)
      ZZ <- length(durationHolder[size,DUR,2,slflw,][!is.na(durationHolder[size,DUR,2,slflw,])]) + 1
      consecDry[size,DUR,slflw,7] <- mean(durationHolder[size,DUR,2,slflw,1:ZZ], na.rm = TRUE)
    }
  }
}



# Create a temporary array for transfer 8 x 8 
temp <- array(data = outPut[,,1,3,2], dim = 64)
dim(temp) <- c(8,8)
rownames(temp) <- c("DDD","DDW","DWD","DWW","WDD","WDW","WWD","WWW")
colnames(temp) <- c("DDD","DDW","DWD","DWW","WDD","WDW","WWD","WWW")

## reorder the outputs into outPut2
titles <- c("DDD","DDW","DWD","DWW","WDD","WDW","WWD","WWW")
outPut2 <- outPut
for (size in 1:4) {
  for (all in 1:20) {
    for (slflw in 1:18) {
      temp <- outPut[,,slflw,size,all]
      outPut2[,,slflw,size,all] <- temp[c(1,2,5,3,4,7,6,8),c(1,5,3,7,4,2,6,8)]
      rownames(outPut2) <- titles[c(1,2,5,3,4,7,6,8)]
      colnames(outPut2) <- titles[c(1,5,3,7,4,2,6,8)]
    }
  }
}

## Concatenate outputs as "mean(S.Dev)" in output2[,,,10] (adults) and output2[,,,11] (all) to 2 decimal places
for (slflw in 1:18) {
  for (size in 1:4) {
    for (r in 1:8) {
      for (c in 1:8) {
        if (!is.nan(as.numeric(outPut2[r,c,slflw,size,2]))) {
          holder1 <- paste(round(as.numeric(outPut2[r,c,slflw,size,2]), digits = 2)," (",round(as.numeric(outPut2[r,c,slflw,size,3]), digits = 2),")")
          outPut2[r,c,slflw,size,10] <- holder1 }  
        if (!is.nan(as.numeric(outPut2[r,c,slflw,size,5]))) {
          holder1 <- paste(round(as.numeric(outPut2[r,c,slflw,size,5]), digits = 2)," (",round(as.numeric(outPut2[r,c,slflw,size,6]), digits = 2),")")
          outPut2[r,c,slflw,size,11] <- holder1  
        }
      }
    }
  }
}










###############################################################################################
##     Graphs for publication
############################################################################################
# this creates the basic plots which were finished in image manipulation software

sHei <- c(7,	7.5,	8,	8.5,	9,	9.5,	10,	10.5,	11)
NSil <- c(0.987951807,	0.963855422,	0.915662651,	0.86746988,	0.78313253,	0.626506024,	0.56626506,	0.469879518,	0.361445783)
RSil <- c(0.795180723,	0.746987952,	0.65060241,	0.542168675,	0.43373494,	0.385542169,	0.277108434,	0.180722892,	0.14457831)

durations <- data.frame(sHei,NSil,RSil)

qplot



## Make a plot 1
pd <- position_dodge(0.2) # move them .05 to the left and right
ggplot(durations, aes(x=Sill, y=Mean, colour=Flow)) +
  theme(axis.line=element_line(color="black",size=0.5)) +
  theme(aspect.ratio=1) +
  geom_errorbar(aes(ymin=Mean-StDev, ymax=Mean+StDev), width=.2, position = pd) +
  # geom_smooth(aes(ymin=Mean-StDev, ymax=Mean+StDev), position=pd) +
  geom_point(position=pd) +
  labs(y= "Successive dry years", x = "Sill height (m)") +
  geom_line(linetype="dashed",aes(x=Sill,y=Max), show.legend = TRUE) 


ggsave("Display1.png",width=5,height=5)

#  geom_density(data=durations,aes(x=Sill, y=Max))


# panel.grid.major.y = element_line(colour = "grey")


ggplot(durations, aes(x=Sill, y=PercWet, fill=Flow,colour=Flow)) +
  geom_density(stat="identity",aes(x=Sill, y=PercWet)) +
  labs(y= "% Wet years", x = "Sill height (m)") +
  theme(aspect.ratio=1) 

ggsave("Display2.png",width=5,height=5)




######START HERE





ggplot(sills, aes(x=Sill, y=Mean, colour=Flow))

##############################################################################################
##  GRAPH 3
## graph 4 size on sill height x PrExt x SillReg(actually PrDry)  
##  2 plots regulated and unreg

#NatExt <- read.csv("C:/Workspace/PrExtNat.csv", sep = ",", dec = ".",col.names = TRUE,rownames = TRUE)
#RegExt <- read.csv("C:/Workspace/PrExtReg.csv", header = TRUE, sep = ",", dec = ".")

NatExt <- read.csv("C:/Workspace/PrExtNat.csv", check.names = FALSE, row.names = 1,header = TRUE)
RegExt <- read.csv("C:/Workspace/PrExtReg.csv", check.names = FALSE, row.names = 1,header = TRUE)
NatExt <- t(NatExt)
RegExt <- t(RegExt)
## Make Individual graphs

#####Graph3a
############
S <- NatExt[,1]
M <- NatExt[,2]
L <- NatExt[,3]
V <- NatExt[,4]
Heights <- NatExt[,5]
Nframe <- data.frame(S,M,L,V,Heights)


natPlot <- ggplot(data=Nframe, aes(x=Heights,y=), colour=Flow) + 
 # ggplot(data=NatExt) +
  theme(axis.line=element_line(color="black",size=1)) +
# scale_x_discrete(breaks = 1:9,labels = Heights) +
  scale_x_continuous(limits = c(7, 11),breaks=seq(7.5,11,.5),  expand = c(0, 0)) +
   scale_y_continuous(limits = c(0, 1),breaks=seq(.25,1,.25), expand = c(0, 0)) +
 # draw_image(Piccy, x = 8, y = -19.23, width = 40, height = 40, scale = 1,clip = "inherit", interpolate = TRUE, hjust = 0, vjust = 0) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank()) + 
  ggtitle("Natural") +
  theme(plot.title = element_text(size=25, hjust=0.5)) +
  labs(y= "",x="") + 
  theme(aspect.ratio=.7) +
  theme(axis.text.y = element_text(size=12, hjust= -0.5)) +
  
  geom_line(linetype="dotted",color="black",size=1,aes(x=Heights,y=S), show.legend = TRUE) +
  geom_line(linetype="dotdash",color="black",size=1,aes(x=Heights,y=M), show.legend = TRUE) +
  geom_line(linetype="dashed",color="black",size=1,aes(x=Heights,y=L), show.legend = TRUE) +
  geom_line(linetype="solid",color="black",size=1,aes(x=Heights,y=V), show.legend = TRUE) 


#####Graph3b
############
S <- RegExt[,1]
M <- RegExt[,2]
L <- RegExt[,3]
V <- RegExt[,4]
Heights <- RegExt[,5]
Rframe <- data.frame(S,M,L,V,Heights)


regPlot <- ggplot(data=Rframe, aes(x=Heights,y=), colour=Flow) + 
  # ggplot(data=NatExt) +
  theme(axis.line=element_line(color="black",size=1)) +
  #scale_x_discrete(breaks = 1:9,labels = Heights) +
   scale_x_continuous(limits = c(7, 11),breaks=seq(7.5,11,.5),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1),breaks=seq(.25,1,.25), expand = c(0, 0)) +
  # draw_image(Piccy, x = 8, y = -19.23, width = 40, height = 40, scale = 1,clip = "inherit", interpolate = TRUE, hjust = 0, vjust = 0) +
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank()) + 
  #theme(axis.title.x=element_blank()) + 
  ggtitle("Regulated") +
  theme(plot.title = element_text(size=25, hjust=0.5)) +
  labs(y= "",x="") + 
  theme(aspect.ratio=0.7) +
  theme(axis.text.y = element_text(size=12, hjust= -0.5)) +
  theme(axis.text.x = element_text(size=12, vjust= -0.5)) +
  geom_line(linetype="dotted",color="black",size=1,aes(x=Heights,y=S), show.legend = TRUE) +
  geom_line(linetype="dotdash",color="black",size=1,aes(x=Heights,y=M), show.legend = TRUE) +
  geom_line(linetype="dashed",color="black",size=1,aes(x=Heights,y=L), show.legend = TRUE) +
  geom_line(linetype="solid",color="black",size=1,aes(x=Heights,y=V), show.legend = TRUE) 



graph3 <- plot_grid(natPlot,regPlot,nrow = 2,align = "h")

y.grob <- textGrob("Probability of Extinction(PrExt)", vjust = 7, hjust = 0.4, gp=gpar(fontface="bold", col="black", fontsize=17), rot=90)

x.grob <- textGrob("Sill height of wetland (m)", vjust = 0, hjust = 0.32, gp=gpar(fontface="bold", col="black", fontsize=17))

#add to plot

grid.arrange(arrangeGrob(graph3, left = y.grob, bottom = x.grob))










##############################################################################################
##  GRAPH 4
## graph 4 size x drought duration x PrExt x SillReg(actually PrWet)
##Make the Matrix
## Order of sillFlows is
sillFlowOrder <- c(1,3,5,7,9,2,4,11,6,13,8,15,10,12,17,14,16,18)
wetProp <- c(98.7,96.3,91.6,86.4,79.5,78.3,74.7,65,62.6,56.6,54.2,46.9,43.3,38.6,36.1,27.7,18,14.5)
graphingTable2 <- array(data = 0, dim = (18*6))
dim(graphingTable2) <- c(18,6)
graphingTable2[,5] <- wetProp
graphingTable2[,6] <- 100 - wetProp

rownames(graphingTable2) <- wetProp
colnames(graphingTable2) <- c("Small", "Medium","Large","Very Large","Prop(W)","Prop(D)")

Piccy <- image_read("C:/Workspace/math0286/R/SBFInkscape.png")


## MAKE GRAPH 1/4
for (size in 1:4) {
  for (cell in 1:18) {
    sillFlowOrder[cell] 
    ## Note the 1 below is DD in consecDry
  graphingTable2[cell,size] <-  consecDry[size,1,sillFlowOrder[cell],5]
 graphingTable2[,1:4]/100
}
}

Small <- graphingTable2[,1]/100
Medium <- graphingTable2[,2]/100
Large <- graphingTable2[,3]/100
Very.Large <- graphingTable2[,4]/100
Wet <- graphingTable2[,5]
Dry <- graphingTable2[,6]


dfGraph2 <- data.frame(Small,Medium,Large,Very.Large,Wet,Dry)
dfGraph2[is.na(dfGraph2)] = 0


P1 <- ggplot(dfGraph2, aes(x=Dry,y=), colour=Flow) + 
theme(axis.line=element_line(color="black",size=1)) +
 scale_x_continuous(limits = c(0, 90),breaks=seq(20,80,20), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1),breaks=seq(.25,1,.25), expand = c(0, 0)) +
 draw_image(Piccy, x = 8, y = -19.3, width = 40, height = 40, scale = 1,clip = "inherit", interpolate = TRUE, hjust = 0, vjust = 0) +
  # theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank()) + 
  geom_text(x=25, y=0.95, label="2y drought",size = 7, fontface = "plain") +
  labs(y= "",x="") + 
  theme(aspect.ratio=1) +
  geom_line(linetype="dotted",color="black",size=1,aes(x=Dry,y=Small), show.legend = TRUE) +
  geom_line(linetype="dotdash",color="black",size=1,aes(x=Dry,y=Medium), show.legend = TRUE) +
  geom_line(linetype="dashed",color="black",size=1,aes(x=Dry,y=Large), show.legend = TRUE) +
  geom_line(linetype="solid",color="black",size=1,aes(x=Dry,y=Very.Large), show.legend = TRUE) 
  

## MAKE GRAPH 2/4
for (size in 1:4) {
  for (cell in 1:18) {
    sillFlowOrder[cell] 
    graphingTable2[cell,size] <-  consecDry[size,2,sillFlowOrder[cell],5]
    graphingTable2[,1:4]/100
  }
}

Small <- graphingTable2[,1]/100
Medium <- graphingTable2[,2]/100
Large <- graphingTable2[,3]/100
Very.Large <- graphingTable2[,4]/100
Wet <- graphingTable2[,5]
Dry <- graphingTable2[,6]


dfGraph2 <- data.frame(Small,Medium,Large,Very.Large,Wet,Dry)
dfGraph2[is.na(dfGraph2)] = 0

P2 <- ggplot(dfGraph2, aes(x=Dry,y=), colour=Flow) + 
  theme(axis.line=element_line(color="black",size=1)) +
  scale_x_continuous(limits = c(0, 90),breaks=seq(20,80,20),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1),breaks=seq(.25,1,.25), expand = c(0, 0)) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) + 
  theme(aspect.ratio=1) +
  labs(y= "", x = "") + 
  geom_text(x=25, y=0.95, label="3y drought",size = 7, fontface = "plain") +
  theme(plot.title = element_text(size=25)) +
  geom_line(linetype="dotted",color="black",size=1,aes(x=Dry,y=Small), show.legend = TRUE) +
  geom_line(linetype="dotdash",color="black",size=1,aes(x=Dry,y=Medium), show.legend = TRUE) +
  geom_line(linetype="dashed",color="black",size=1,aes(x=Dry,y=Large), show.legend = TRUE) +
  geom_line(linetype="solid",color="black",size=1,aes(x=Dry,y=Very.Large), show.legend = TRUE) 




## MAKE GRAPH 3/4
for (size in 1:4) {
  for (cell in 1:18) {
    sillFlowOrder[cell] 
    graphingTable2[cell,size] <-  consecDry[size,3,sillFlowOrder[cell],5]
    graphingTable2[,1:4]/100
  }
}

Small <- graphingTable2[,1]/100
Medium <- graphingTable2[,2]/100
Large <- graphingTable2[,3]/100
Very.Large <- graphingTable2[,4]/100
Wet <- graphingTable2[,5]
Dry <- graphingTable2[,6]


dfGraph2 <- data.frame(Small,Medium,Large,Very.Large,Wet,Dry)
dfGraph2[is.na(dfGraph2)] = 0

P3 <- ggplot(dfGraph2, aes(x=Dry,y=), colour=Flow) + 
  theme(axis.line=element_line(color="black",size=1)) +
  scale_x_continuous(limits = c(0, 90),breaks=seq(20,80,20),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1),breaks=seq(.25,1,.25), expand = c(0, 0)) +
  theme(aspect.ratio=1) +
  labs(y= "",x="",size=10) + 
  geom_text(x=65, y=0.1, label="4y drought",size = 7, fontface = "plain") +
  theme(plot.title = element_text(size=25)) +
  geom_line(linetype="dotted",color="black",size=1,aes(x=Dry,y=Small), show.legend = TRUE) +
  geom_line(linetype="dotdash",color="black",size=1,aes(x=Dry,y=Medium), show.legend = TRUE) +
  geom_line(linetype="dashed",color="black",size=1,aes(x=Dry,y=Large), show.legend = TRUE) +
  geom_line(linetype="solid",color="black",size=1,aes(x=Dry,y=Very.Large), show.legend = TRUE) 


## MAKE GRAPH 4/4
for (size in 1:4) {
  for (cell in 1:18) {
    sillFlowOrder[cell] 
    graphingTable2[cell,size] <-  consecDry[size,4,sillFlowOrder[cell],5]
    graphingTable2[,1:4]/100
  }
}

Small <- graphingTable2[,1]/100
Medium <- graphingTable2[,2]/100
Large <- graphingTable2[,3]/100
Very.Large <- graphingTable2[,4]/100
Wet <- graphingTable2[,5]
Dry <- graphingTable2[,6]

dfGraph2 <- data.frame(Small,Medium,Large,Very.Large,Wet,Dry)
dfGraph2[is.na(dfGraph2)] = 0

Piccy4b <- image_read("C:/Workspace/math0286/R/SBFLegend4.png")


P4 <- ggplot(dfGraph2, aes(x=Dry,y=), colour=Flow) + 
  theme(axis.line=element_line(color="black",size=1)) +
  scale_x_continuous(limits = c(0, 90),breaks=seq(20,80,20),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1),breaks=seq(.25,1,.25), expand = c(0, 0)) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) +
  draw_image(Piccy4b, x = 50, y = -19.6, width = 40, height = 40, scale = 1.5,clip = "inherit", interpolate = TRUE, hjust = 0.25, vjust = 0) +
 # theme(legend.position = c(0.95,0.95)) +
  theme(aspect.ratio=1) +
  labs(y= "", x = "") + 
  geom_text(x=65, y=0.1, label="5y drought",size = 7, fontface = "plain") +
  geom_line(linetype="dotted",color="black",size=1,aes(x=Dry,y=Small), show.legend = TRUE) +
  geom_line(linetype="dotdash",color="black",size=1,aes(x=Dry,y=Medium), show.legend = TRUE) +
  geom_line(linetype="dashed",color="black",size=1,aes(x=Dry,y=Large), show.legend = TRUE) +
  geom_line(linetype="solid",color="black",size=1,aes(x=Dry,y=Very.Large), show.legend = TRUE) 



#grid.arrange(P1,P2,P3,P4,nrow = 2,bottom="Frequency of Dry years at the site (%)", left = "Probability of extinction(PrExt)")





graph4 <- plot_grid(P1,P2,P3,P4,nrow = 2,align = "hv")

y.grob <- textGrob("extinction probability", vjust = 3, gp=gpar(col="black", fontsize=17), rot=90)

x.grob <- textGrob("background frequency of dry years", vjust = -0.5, gp=gpar(col="black", fontsize=17))

#add to plot

grid.arrange(arrangeGrob(graph4, left = y.grob, bottom = x.grob))







fg <- c(25,138,629,3048,7495)
demog <- array (data = fg, dim = 5)
dim(demog) <- c(5,1)
rownames(demog) <- c('5-6','4-5','3-4','2-3','1-2')
  











































df <- data.frame(x=1:10, a=rnorm(10), b=rnorm(10), c=rnorm(10))
mdf <- reshape2::melt(df, id.var = "x")







ggplot(dfGraph2, aes(x=Dry,y=), show.legend=TRUE) + 
  theme(axis.line=element_line(color="black",size=1)) +
  scale_x_continuous(limits = c(0, 90),breaks=seq(20,80,20),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1),breaks=seq(.25,1,.25), expand = c(0, 0)) +
  #theme(axis.title.y=element_blank(),axis.text.y=element_blank()) +
  #theme(legend.position = c(0.5,0.5)) +
  theme(aspect.ratio=1) +
  labs(y= "", x = "") + 
  geom_text(x=65, y=0.1, label="5y drought",size = 7, fontface = "plain") +
  geom_line(linetype="dotted",color="black",size=1,aes(x=Dry,y=Small), show.legend = TRUE) +
  geom_line(linetype="dotdash",color="black",size=1,aes(x=Dry,y=Medium), show.legend = TRUE) +
  geom_line(linetype="dashed",color="black",size=1,aes(x=Dry,y=Large), show.legend = TRUE) +
  geom_line(linetype="solid",color="black",size=1,aes(x=Dry,y=Very.Large), show.legend = TRUE)
  
    
    
    
    scale_linetype_manual("", breaks = c("Small", "Medium", "Large","Very Large"), values = c("dotted", "dotdash", "dashed","solid")) +
theme(legend.position = c(0.5,0.5))



plot(dfGraph2)
dfGraph2


dfGraph2 <- data.frame(Dry,Small,Medium,Large,Very.Large)
dfGraph2[is.na(dfGraph2)] = 0
rownames(dfGraph2) <- Dry

df <- data.frame(Dry=Dry, Small=rnorm(18), Medium=rnorm(18), Large=rnorm(18), Very.Large=rnorm(18))
mdf <- reshape2::melt(df, id.var = "Dry")


plot(mdf)

mdf


ggplot(mdf, aes(x=Dry,y=value), show.legend=TRUE, colour = variable)



ggplot(mdf, aes(x=Dry,y=Val)) + 
  geom_line(linetype="dotted",color="black",size=1,aes(x=Dry,y=Small), show.legend = TRUE) +
  geom_line(linetype="dotdash",color="black",size=1,aes(x=Dry,y=Medium), show.legend = TRUE) +
  geom_line(linetype="dashed",color="black",size=1,aes(x=Dry,y=Large), show.legend = TRUE) +
  geom_line(linetype="solid",color="black",size=1,aes(x=Dry,y=Very.Large), show.legend = TRUE)






ggplot(data = dfGraph2 , aes(x = Dry)) +
  geom_line(aes(x=Dry,y = Small), linetype="dotted") +
  geom_line(aes(x=Dry,y = Medium), linetype="dotdash") +
  geom_line(aes(x=Dry,y = Large), linetype="dashed") +
  geom_line(aes(x=Dry,y = Very.Large), linetype="solid") +
  scale_linetype_manual("", breaks = c("Small", "Medium", "Large","Very.Large"), values = c("dotted", "dotdash", "dashed","solid")) +
  xlab(" ") 















