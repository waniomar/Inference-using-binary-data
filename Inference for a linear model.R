#######################################################
#######################################################
# Small didactic example for: inference using likelihood function for binary observations

# O.Wani, Feb 2019

# For more details on the method, please refer to:
# https://doi.org/10.1016/j.watres.2017.05.038
#######################################################
#######################################################


#################################
# Set library address
#################################

#assign(".lib.loc", "E:\\02 Office\\10 R\\library", envir = environment(.libPaths))

#################################
# Get libraries
#################################

if ( !require(mvtnorm) )    { install.packages("mvtnorm");    library(mvtnorm)} # for computing  normal distribution for arbitrary limits
if ( !require(IDPmisc) )    { install.packages("IDPmisc");    library(IDPmisc) } # for bivariate plots
if ( !require(adaptMCMC) )  { install.packages("adaptMCMC");    library(adaptMCMC) } # for posterior probability density sampling 

#################################
# Define model and data
#################################

model=function(par,P){a=par[1]*P+par[2] 
return(a)} #par is the parameter vector, P is the input vector


#################################
# Generate data
#################################

P=sin(seq(1,200,1)*0.1)+1 # define input time series



# define a covariance functionbetween two time points t1 and t2
cov.exp <- function(t1, t2, par) {
  r <- abs(t1 - t2)
  a=par["B"]^2*exp(-r/(2*par["c"]))
  return(a)
}

# define an observation generating process

likelihood=function(par){
  Y.det=model(par[1:2],P)
  ## construct covariance of observation points
  t1=seq(1,length(Y.det),1)
  t2=seq(1,length(Y.det),1)
  n1 <- length(t1)
  n2 <- length(t2)
  Sigma <- matrix(NA, nrow=n1, ncol=n2)
  
  for(j in 1:n2) {
    for(i in 1:n1) {
      Sigma[i, j] <- cov.exp(t1[i], t2[j], par)
    }
  }
  
  Sigma= Sigma+diag(rep((par["E"])^2,length(Y.det)),nrow=length(Y.det),ncol=length(Y.det))
  likeli=mvrnorm(n = 1, mu=Y.det, Sigma=Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  
  return(likeli)
}

# generate an observation time series
data=likelihood(c(1,1,E=0.1,B=0.3,c=1))

# define a threshold to generate binary observations
thresh=2.5

out.b=c() # vector for binary observations

for(i in 1:length(data))
{
  
  if (data[i]>=thresh){out.b[i]=T}
  else {out.b[i]=F}
  
}

# define a subset of data for use in inference
Pn=P[1:100]
datan=data[1:100]
out.bn=out.b[1:100]

#######################################
# define prior parameter distribution
#######################################


prior.pbdis<- list( p1   =c("NormalTrunc", 2.5,1,0.25,4),
                    p2   =c("NormalTrunc", 2.5,1,0.25,4),
                    E    =c("NormalTrunc", 0.08,0.1,0,0.2),
                    B    =c("NormalTrunc", 0.1,0.3,0,2),
                    c     =c("NormalTrunc", 1.3,1,0.1,2)) 


log.prior=function(x,distpar){
  # truncated normal distribution; parameters are mean, sd, min and max
  prob=1
  names(prob)="log.probability"
  for (i in 1:length(names(distpar))){
    mean <- as.numeric(distpar[[i]][2])
    sd   <- as.numeric(distpar[[i]][3])
    min  <- as.numeric(distpar[[i]][4])
    max  <- as.numeric(distpar[[i]][5])
    fact <- 1/(pnorm(q=max,mean=mean,sd=sd)-pnorm(q=min,mean=mean,sd=sd))
    prob=prob*ifelse(x[i]<min|x[i]>max,0,fact*dnorm(x[i],mean=mean,sd=sd))
  }
  return(log(prob)[1])
}


# check .......
par.init=c(p1=2.5,p2=2.5,E=0.08,B=0.1,c=1.3)
log.prior(x=par.init,distpar=prior.pbdis)

##########################################
## The likelidood function
##########################################

log.likelihood<-function(par,out.b) {
  Y.det=model(par[1:2],Pn)
  ## construct covariance of observation points
  t1=seq(1,length(Y.det),1)
  t2=seq(1,length(Y.det),1)
  n1 <- length(t1)
  n2 <- length(t2)
  Sigma <- matrix(NA, nrow=n1, ncol=n2)
  
  for(j in 1:n2) {
    for(i in 1:n1) {
      Sigma[i, j] <- cov.exp(t1[i], t2[j], par)
    }
  }
  
  Sigma= Sigma+diag(rep((par["E"])^2,length(Y.det)),nrow=length(Y.det),ncol=length(Y.det))
  
  Sigma <- Sigma + sqrt(.Machine$double.eps)*diag(n1) # for numerical stability
  
  ## integration boudaries
  lower <- upper <- rep(thresh, length(out.b))
  lower[!out.b] <- -Inf
  upper[out.b] <- Inf
  
  prob = try(log(pmvnorm(lower, upper, mean=Y.det, sigma=Sigma))[1])
  
  if(class(prob) == "try-error") {
    print("-----------")
    print(par)
    print(prob)
    prob <- -Inf
  }
  
  return(prob)
}

#test
log.likelihood(par=par.init,out.b = out.bn)

logposterior<- function(par) # for classical error models
{names(par)=names(par.init)
if(log.prior(x=par,distpar=prior.pbdis)==-Inf){out=-Inf}
else {out=try(log.likelihood(par,out.b = out.bn)+log.prior(x=par,distpar=prior.pbdis))}
if(!is.finite(out)){out=-Inf}

if(rnorm(1, mean = 0, sd = 1)>1.9){  # to monitor the progress during iterations
  print(paste("log post: ", format(out, digits=2)));     print(par)}
return(out)}

#################################
# MCMC
#################################

mod.runs = 5000

RAM <- MCMC(p     = logposterior, 
            init  = c(1.5,1.5,0.2,0.4,1.5),
            scale = (c(2,2,0.3,0.5,1))^2,
            n     = mod.runs, 
            adapt = T,
            acc.rate = 0.3,
            gamma =0.7,
            n.start = 100)


samples=RAM$samples[1000:mod.runs,]
par.optim.m=samples[which.max(RAM$log.p[1000:mod.runs]),] #parameter values corresponding to the maximum posterior probability

# comparisom between the parameter values of the prior and posterior maximum density 
par.optim.m
par.init

# plot the model output updated parameter values
par(mfrow = c(1,1))
plot(y=datan,x=seq(1,100,1), ylab="System Response [--]",xlab="Time[--]",ylim=c(0,10))
lines(model(par=par.optim.m[1:2],Pn), col="red")
lines(model(par=par.init[1:2],Pn),col="blue")
legend("topright",legend = c(paste("Max. posterior"),paste("Max. Prior")),horiz=F, col=c("red","blue"),lty=c(1,1),lwd=c(2,2),pch=c(NA,NA),cex=1.2,bty = "n",pt.cex = c(1,1))

# plot posterior 
ipairs(samples)
