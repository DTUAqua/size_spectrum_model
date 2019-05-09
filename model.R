library(TMB)
compile("model.cpp")
dyn.load(dynlib("model"))

## Get data
## N = matrix of counts
## times = survey times
## sizes = midpoint of size intervals (assuming equidistant cells)
## N[i,j] = Number caught of size 'size[i]' at time 'time[j]'
dat <- dget("dat.txt")
matplot(dat$sizes, dat$N, type="l")
legend("topright", lty = 1:5, col=1:4, legend=round(dat$times,1))

## Survey period is 1998 to 2012.
## Choose recruitment year-classes needed to model the data
t0 <- (1985:2011)-.5

## List of parameters - see cpp file for explanation
parameters <- list(
    L0 = 0,
    meanLinf = 50,
    sdLinf = 4,
    k = .14, ## Growth = ca. 7 cm per year (k * meanLinf)
    beta = .7,
    l50f = 20,
    M0 = .2,
    l50 = 12,
    alpha = .6,
    minsel = 0,
    logp = 0,
    Finf = 1,
    logR = rep(8.5,length(t0)),
    t0 = t0,
    dt0 = -0.3,
    spawnsd = rep(.03,length(t0))
)

## Map to fix all parameters
map <- lapply(parameters, function(x)factor(x*NA))

## Free selected ones
map$k <- NULL ## Free
map$sdLinf <- NULL ## Free
map$logR <- factor(rep(1, length(parameters$logR))) ## Constant recruitment

obj <- MakeADFun(
    data = dat,
    parameters = parameters,
    map = map                 
)

## Fit model
opt <- nlminb(obj$par, obj$fn, obj$gr)

rep <- obj$report()
sdrep <- sdreport(obj)

## Fitted size distributions
matplot(rep$mat, type="l")

## Compare with observations
layout(matrix(1:24, 4))
par(mar=c(2,2,2,2))
for(i in 1:length(dat$times)){
    matplot(dat$sizes,cbind(rep$N[,i],rep$mat[,i]),type="b",xlab="",ylab="")
    title(round(dat$times[i],2))
}
