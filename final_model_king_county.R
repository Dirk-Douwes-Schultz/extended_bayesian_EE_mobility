# Final model for King County
rm(list=ls())
setwd("~/Github/extended_bayesian_EE_mobility")
library(nimble)
library(dlnm)
library(dplyr)
library(tidyr)
library(coda)

load("king_county_data.RData")


#add weekday
king_cases_home$wkday <- 1
king_cases_home[king_cases_home$weekdays=="Tuesday",]$wkday <- 2 
king_cases_home[king_cases_home$weekdays=="Wednesday",]$wkday <- 3
king_cases_home[king_cases_home$weekdays=="Thursday",]$wkday <- 4
king_cases_home[king_cases_home$weekdays=="Friday",]$wkday <- 5
king_cases_home[king_cases_home$weekdays=="Saturday",]$wkday <- 6
king_cases_home[king_cases_home$weekdays=="Sunday",]$wkday <- 7

#get time of first nonzero case 
minpos <- min(which(king_cases_home$cases_daily>0))

##deal with the zeros
zero_index <- which(king_cases_home[minpos:nrow(king_cases_home),]$cases_daily==0)+minpos-1
for(j in zero_index){
  print(j)
  king_cases_home$cases_daily[j] <- round(mean(king_cases_home$cases_daily[(j+1):(j+5)]))
  for(i in 1:5){
    king_cases_home$cases_daily[j+i] <- round((4/5)*king_cases_home$cases_daily[j+i])
  }
}

#stay at home spline (from the paper, max_lag=L and min_lag=L_0)
max_lag <- 14
min_lag <- 7

cb1.ca <- crossbasis(king_cases_home$StayAtHomeRate[(minpos-max_lag):nrow(king_cases_home)], lag=c(min_lag,max_lag), argvar=list(fun="lin"), arglag=list(fun="bs", degree=3))
summary(cb1.ca)

##set mT (time model is fit up to) and K (steps ahead for forecasting)
mT <- 320
K <- 7

#JAGS code
cases <- king_cases_home$cases_daily[minpos:nrow(king_cases_home)]
wkday <- king_cases_home$wkday[minpos:nrow(king_cases_home)]
stay_home_rate <- king_cases_home$StayAtHomeRate[minpos:nrow(king_cases_home)]

#need a weekday matrix
wkdaym <- matrix(ncol=6,nrow=length(wkday),data=rep(0,6*length(wkday)))
for(i in 1:length(wkday)){
  if(wkday[i]!=1){
    wkdaym[i,wkday[i]-1] <- 1
  }
}

#need monday effect
monday <- rep(0,length(wkday))
monday[wkday==1] <- 1


#Make week variable
week <- rep(NA,length(wkday))
week[1] <- 1
for(j in 2:length(wkday)){
  
  if(wkday[j]==1){
    week[j] <- week[j-1]+1
  }else{
    week[j] <- week[j-1]
  }
}

#just check it
#View(cbind(wkday,week))

#final week for the random walk
maxweek <- week[mT+K]

# center the spline to reduce posterior correlations
stay_spline <- scale(cb1.ca[(max_lag+1):nrow(cb1.ca),],scale = FALSE,center = TRUE)

len <- ncol(stay_spline)

#need to make custom distribution
dmixunif <- nimbleFunction(
  run = function(x = double(0), p = double(0, default = 7), L_max = double(0, default = 10),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    if(0<=x & x<=1){
      logProb <- -log(p)
    }else if(1<x & x<=p){
      logProb <- log(p-2)-log(p)-log(p-1)
    }else if(p<x & x<=L_max){
      logProb <- -log(p)-log(L_max-p)
    }else{
      logProb <- -Inf
    }
    if(log) return(logProb)
    else return(exp(logProb)) 
  })

rmixunif <- nimbleFunction(
  run = function(n = integer(0),  p = double(0, default = 7), L_max = double(0, default = 10)) {
    returnType(double(0))
    if(n != 1) print("rmixunif only allows n = 1; using n = 1.")
    dev <- runif(n=1, min=0, max=1)
    if(dev<(1/p)){
      xsim <- runif(n=1,min=0,max=1)
    }else if((1/p)<dev & dev < ((p-1)/p)){
      xsim <- runif(n=1,min=1,max=p)
    }else{
      xsim <- runif(n=1,min=p,max=L_max)
    }
    return(xsim)
  })

#test
sims <- 1000000
q_prior <- rep(NA,sims)
for(i in 1:sims){
  q_prior[i] <- rmixunif(n=1)
}
hist(q_prior,freq = FALSE)
sum(0<q_prior & q_prior<1)/sims
dmixunif(.5)
(sum(1<q_prior & q_prior<7)/6)/sims
dmixunif(4)
(sum(7<q_prior & q_prior<10)/3)/sims
dmixunif(9)
#it appears accurate

#---- nimble code ----
nimlbe.data <- list(y=cases)

nimble.const <- list(T=mT,wkday=wkday,stay_spline=stay_spline, 
                     len=ncol(stay_spline),K=K,week=week,maxweek=maxweek,
                     stay_home_rate=stay_home_rate,monday=monday)

  
  
  assign('dmixunif', dmixunif, envir = .GlobalEnv)
  assign('rmixunif', rmixunif, envir = .GlobalEnv)
  
  #note not all variable names are the same as in the paper
  nimble.Code <- nimbleCode({
    
    ## Specify data likelihood
    #weights
    for(d in 1:7){
      #shifted negative binomial
      w[d] <- (exp(loggam(d-1+r))/(exp(logfact(d-1))*exp(loggam(r))))*pow(1-kappa,r)*pow(kappa,d-1)
      weight[d] <- w[d]/sum(w[1:7])
      
      #prior distribution of the weights
      w_prior[d] <- (exp(loggam(d-1+r_prior))/(exp(logfact(d-1))*exp(loggam(r_prior))))*pow(1-kappa_prior,r_prior)*pow(kappa_prior,d-1)
      weight_prior[d] <- w_prior[d]/sum(w_prior[1:7])
    }
    #likelihood
    for(t in 8:T){
      y[t] ~ dnegbin(p[t],r2[t])
      # variance_t = u[t](1+u[t]/r2[t])
      p[t] <- r2[t]/(r2[t]+u[t])
      u[t] <- v[t]+phi[t]*(y[t-1]*weight[1]+y[t-2]*weight[2]+y[t-3]*weight[3]+y[t-4]*weight[4]+y[t-5]*weight[5]+y[t-6]*weight[6]+y[t-7]*weight[7])
      fitted[t] ~ dnegbin(p[t],r2[t])
      #endemic component
      log(v[t]) <-  beta[1] 
      #reproduction number
      log(phi[t]) <- moneff*monday[t]+inprod(beta[2:(len+1)],stay_spline[t,1:(len)])+weekeff[week[t]] 
      r2[t] <- alpha[1]
      res[t]<-(y[t]-u[t])/sqrt(u[t]*(1+u[t]/r2[t]))
    }
    
    ##priors 
    moneff ~ dnorm(0,1/100)
    #weekeff
    weekeff[1] ~ dnorm(0,1)
    for(wk in 2:maxweek){
      weekeff[wk] ~ dnorm(weekeff[wk-1],sd=sigma_wk)
    }
    #
    beta[1] ~  dunif(-2,50)
    for(d in 2:(len+1)){
      beta[d] ~ dnorm(0,1/100)
    }
    for(d in 1:1){
      inverse_alpha[d] ~T(dt(0, 1, 1),0,)
      alpha[d] <- 1/inverse_alpha[d]
    }
    #negative binomial weights
    kappa ~ dunif(0,1)
    r ~ dmixunif(p=7,L_max=10)
    r_prior ~ dmixunif(p=7,L_max=10)
    kappa_prior ~ dunif(0,1)
    #random walk
    sigma_wk ~ T(dt(0, 1, 1),0,)
  })
  
  inits <- list(beta= c(runif(n=1,min=-2,max=5),rnorm(n=(len),mean=0,sd=1)),
                moneff = rnorm(n=1,mean=0,sd=1),
                weekeff = c(0,rnorm(n=(maxweek-1),mean=0,sd=1)),
                kappa = runif(n=1,min=0,max=1),
                r = runif(n=1,min=0,max=1),
                sigma_wk = runif(n=1,min=0,max=5),
                inverse_alpha = runif(n=1,min=0,max=50))
  
  nimble.model <- nimbleModel(nimble.Code, nimble.const, nimlbe.data, inits)
  
  Cnimble.model <- compileNimble(nimble.model)
  
  nimble.Conf <- configureMCMC(nimble.model, print = TRUE)
  
  nimble.Conf$removeSampler(c("beta[2:5]"))
  nimble.Conf$addSampler(target=c("beta[2:5]"),type="AF_slice")
  
  nimble.Conf$removeSampler(c("moneff"))
  nimble.Conf$addSampler(target=c("moneff"),type="slice")
  
  nimble.Conf$removeSampler(c("kappa","r"))
  nimble.Conf$addSampler(target=c("kappa","r"),type="AF_slice")
  
  nimble.Conf$removeSampler(c("beta[1]"))
  nimble.Conf$addSampler(target=c("beta[1]"),type="slice")
  
  nimble.Conf$removeSampler(c("sigma_wk"))
  nimble.Conf$addSampler(target=c("sigma_wk"),type="slice")
  
  print(nimble.Conf)
  
  nimble.Conf$addMonitors(c("fitted","weight","p","r2","u","res","phi","r_prior","weight_prior","alpha"))
  
  nimble.MCMC <- buildMCMC(nimble.Conf)
  
  Cnimble.MCMC <- compileNimble(nimble.MCMC, project = nimble.model,resetFunctions = TRUE)
  
  initsFunction <- function()  list(beta= c(runif(n=1,min=-2,max=5),rnorm(n=(len),mean=0,sd=1)),
                                    moneff = rnorm(n=1,mean=0,sd=1),
                                    weekeff = c(0,rnorm(n=(maxweek-1),mean=0,sd=1)),
                                    kappa = runif(n=1,min=0,max=1),
                                    r = runif(n=1,min=0,max=1),
                                    sigma_wk = runif(n=1,min=0,max=5),
                                    inverse_alpha = runif(n=1,min=0,max=50))
  
  
  start_time <- Sys.time()
  samples <- runMCMC(Cnimble.MCMC,niter = 100000,nchains = 3,nburnin=10000
                     ,samplesAsCodaMCMC = TRUE,thin=9,inits = initsFunction)
  
  end_time <- Sys.time()
  end_time - start_time
  
  
#check convergence
gelman.diag(samples[,c("r","kappa","beta[2]",
                         "beta[3]","beta[4]","beta[5]",
                         "weekeff[2]","beta[1]","moneff",
                       "alpha[1]","sigma_wk","weekeff[1]")])

effectiveSize(samples[,c("r","kappa","beta[2]",
                       "beta[3]","beta[4]","beta[5]",
                       "weekeff[2]","beta[1]","moneff",
                       "alpha[1]","sigma_wk","weekeff[1]")])

min(effectiveSize(samples[,c("r","kappa","beta[2]",
                             "beta[3]","beta[4]","beta[5]",
                             "weekeff[2]","beta[1]","moneff",
                             "alpha[1]","sigma_wk","weekeff[1]")]))

#check some random values of the random walk
gelman.diag(samples[,c("weekeff[1]","weekeff[2]","weekeff[9]","weekeff[21]",
                       "weekeff[28]","weekeff[36]","weekeff[42]","weekeff[46]","weekeff[47]")])


#tracplots
plot(samples[,c("r","kappa")])
plot(samples[,c("beta[1]","beta[2]")])
plot(samples[,c("beta[3]","beta[4]","beta[5]")])
plot(samples[,c("sigma_wk","alpha[1]")])
plot(samples[,c("moneff")])
plot(samples[,c("weekeff[1]","weekeff[2]","weekeff[3]")])
plot(samples[,c("weekeff[26]","weekeff[27]","weekeff[28]")])
plot(samples[,c("weekeff[46]","weekeff[47]","weekeff[48]")])

#calc WAIC
samps <- data.frame(rbind(samples[[1]],samples[[2]],samples[[3]]))
lppd <- 0 
pwaic <- 0
for(t in 8:mT){
  lppd <- lppd + log(mean(dnbinom(x = cases[t], size=as.numeric(unlist(samps[paste0(paste0("r2.",t,"."))])), prob=as.numeric(unlist(samps[paste0(paste0("p.",t,"."))])))))
  pwaic <- pwaic + var(log(dnbinom(x = cases[t], size=as.numeric(unlist(samps[paste0(paste0("r2.",t,"."))])), prob=as.numeric(unlist(samps[paste0(paste0("p.",t,"."))])))))
}

waic <- -2*(lppd-pwaic)

#waic = 3151

#plots
#function for error bars
add.error.bars <- function(X,upper,lower,width,col=par( )$fg,lwd=1){
  segments(X,lower,X,upper,col=col,lwd=lwd,lend=1);
  segments(X-width/2,lower,X+width/2,lower,col=col,lwd=lwd,lend=1);
  segments(X-width/2,upper,X+width/2,upper,col=col,lwd=lwd,lend=1);
}

#Pearson residuals
mean_pears <- apply(samps[,grepl( "res" , names( samps ) ) ],MARGIN = 2,mean)
plot(mean_pears)
acf(mean_pears[8:length(mean_pears)])

#need to make predictions
K 
LS_d <- data.frame(mT=mT,k1=NA,k2=NA,k3=NA,k4=NA,k5=NA,k6=NA,k7=NA)
forward_preds <- matrix(nrow=30000,ncol=K)
rn_forward <- matrix(nrow=30000,ncol=K)
for(k in 1:K){
  r2 <- as.numeric(unlist(samps["alpha.1."]))
  v <- exp(as.numeric(unlist(samps["beta.1."])))
  phi <- exp(as.numeric(unlist(samps["moneff"]*monday[mT+k]+
                                 samps["beta.2."]*stay_spline[mT+k,1]+samps["beta.3."]*stay_spline[mT+k,2]+
                                 samps["beta.4."]*stay_spline[mT+k,3]+samps["beta.5."]*stay_spline[mT+k,4]+
                                 samps[paste0("weekeff.",week[mT+k],".")])))
  sum_weight_cases <- rep(0,30000)
  for(d in 1:7){
    if(mT+k-d <= mT){
      sum_weight_cases <- sum_weight_cases + cases[mT+k-d]*as.numeric(unlist(samps[paste0("weight.",d,".")]))
    }else{
      sum_weight_cases <- sum_weight_cases + forward_preds[,k-d]*as.numeric(unlist(samps[paste0("weight.",d,".")]))
    }
  }
  u <- v+phi*sum_weight_cases
  forward_preds[,k] <- rnbinom(n = 30000,size=r2,prob = r2/(r2+u))
  LS <- -log(mean(dnbinom(x=cases[mT+k],prob = r2/(r2+u),
                          size=r2)))
  LS_d[1,k+1] <- LS
  rn_forward[,k] <- phi
}

med = apply(samps[,grepl( "fitted" , names( samps ) ) ],MARGIN = 2,mean)
upper= apply(samps[,grepl( "fitted" , names( samps ) ) ],MARGIN = 2,function(x) quantile(x,probs=c(.975),na.rm=TRUE))
lower= apply(samps[,grepl( "fitted" , names( samps ) ) ],MARGIN = 2,function(x) quantile(x,probs=c(.025),na.rm=TRUE))

medf <- apply(forward_preds,MARGIN = 2,median)
upperf= apply(forward_preds,MARGIN = 2,function(x) quantile(x,probs=c(.975)))
lowerf= apply(forward_preds,MARGIN = 2,function(x) quantile(x,probs=c(.025))) 

par(mfrow=c(1,1))
plot.new()
par(usr=c(-3,length(cases)+3, -5,max(c(upper,upperf))+mean(cases))); #Set the limits for the plotting window axes.
axis(1); #Draw the horizontal axis.
axis(2,las=1); #Draw the vertical axis.
box(bty=letters[12], lwd=2); #Draw a frame around the plot.
title(main="Fitted Daily Covid-19 Cases in King County",cex.main=2,font.main=1,xlab = "Day since start of epidemic",ylab="Cases"); #Add a title
points(1:length(cases),cases,type = "p")
lines(8:mT,med[8:mT],col="red")
lines(8:mT,lower[8:mT],col="red",lty=2)
lines(8:mT,upper[8:mT],col="red",lty=2)
lines((mT+1):(mT+K),medf,col="blue")
lines((mT+1):(mT+K),lowerf,col="blue",lty=2)
lines((mT+1):(mT+K),upperf,col="blue",lty=2)

#plot of the weights
medw = apply(samps[,grepl( "weight" , names( samps ) ) ],MARGIN = 2,median)
upperw= apply(samps[,grepl( "weight" , names( samps ) ) ],MARGIN = 2,function(x) quantile(x,probs=c(.975)))
lowerw= apply(samps[,grepl( "weight" , names( samps ) ) ],MARGIN = 2,function(x) quantile(x,probs=c(.025)))

plot.new()
par(usr=c(0,8, 0,max(upperw)+mean(medw))); #Set the limits for the plotting window axes.
axis(1,at = 1:7); #Draw the horizontal axis.
axis(2,las=1); #Draw the vertical axis.
box(bty=letters[12], lwd=2); #Draw a frame around the plot.
title(main="Weights King County \n (prior=blue, posterior=red)",cex.main=2,font.main=1,xlab = "lag",ylab="weight"); #Add a title
points(1:7,medw[8:14],col="blue")
add.error.bars(X = 1:7,upper = upperw[8:14],lower = lowerw[8:14],width = .5,col="blue")
points(1:7,medw[1:7],col="red")
add.error.bars(X = 1:7,upper = upperw[1:7],lower = lowerw[1:7],width = .5,col="red")
#legend(40,25,legend = c("Posterior Predictive Median","Posterior Predictive 95% CI","Observed"),lty=c(1,2,NA),col=c("red","red","black"),pch=c(NA,NA,1))

#plot of reproduction number
par(mfrow=c(1,1))
medr = apply(samps[,grepl( "phi" , names( samps ) ) ],MARGIN = 2,median)
upperr= apply(samps[,grepl( "phi" , names( samps ) ) ],MARGIN = 2,function(x) quantile(x,probs=c(.975),na.rm=TRUE))
lowerr= apply(samps[,grepl( "phi" , names( samps ) ) ],MARGIN = 2,function(x) quantile(x,probs=c(.025),na.rm=TRUE))

medfrr <- apply(rn_forward,MARGIN = 2,median)
upperfrr= apply(rn_forward,MARGIN = 2,function(x) quantile(x,probs=c(.975)))
lowerfrr= apply(rn_forward,MARGIN = 2,function(x) quantile(x,probs=c(.025))) 

plot.new()
par(usr=c(-3,length(cases)+3, 0,max(upperr)+mean(medr,na.rm = TRUE))); #Set the limits for the plotting window axes.
axis(1); #Draw the horizontal axis.
axis(2,las=1); #Draw the vertical axis.
box(bty=letters[12], lwd=2); #Draw a frame around the plot.
title(main="Reproduction number phi[t]",cex.main=2,font.main=1,xlab = "Day since start of epidemic",ylab="phi"); #Add a title
lines(8:(mT),medr[8:mT],col="red")
lines(8:(mT),lowerr[8:mT],col="red",lty=2)
lines(8:(mT),upperr[8:mT],col="red",lty=2)
lines((mT+1):(mT+K),medfrr,col="blue")
lines((mT+1):(mT+K),lowerfrr,col="blue",lty=2)
lines((mT+1):(mT+K),upperfrr,col="blue",lty=2)
abline(h=1)


##rr plots

# median
beta.median <- apply(samps[,grepl("beta" ,names(samps))], MARGIN = 2, median)
beta.upper <- apply(samps[,grepl("beta" ,names(samps))], MARGIN = 2, function(x) quantile(x,probs=c(.975)))
beta.lower <- apply(samps[,grepl("beta" ,names(samps))], MARGIN = 2, function(x) quantile(x,probs=c(.025)))

# create an identity matrix LxL to get the C matrix for RR plot
L <- max_lag-min_lag+1
sim.var <- c(rep(0, L-1), 1, rep(0, L-1), rep(0,times=min_lag))
W = crossbasis(sim.var, lag=c(min_lag,max_lag), argvar=list(fun="lin"), arglag=list(fun="bs",  degree=3))
C = W[complete.cases(W),]


# lag effects estimates
e.median <- C %*% beta.median[2:length(beta.median)]
e.upper <- C %*% beta.upper[2:length(beta.median)]
e.lower <- C %*% beta.lower[2:length(beta.median)]


# RR plot
par(mfrow=c(1,1))
plot(min_lag:max_lag,exp(e.median*0.1),col="red",pch=19,cex=0.7,ylim = c(0.4, 2), main=paste(loc_name$name[loc_name$id==king_id]," id=",king_id,", lag:",min_lag,"-",max_lag, sep=""))
abline(h=1, lty=2)
add.error.bars(X = min_lag:max_lag,upper = exp(e.upper*0.1),lower = exp(e.lower*0.1),width = .5)
