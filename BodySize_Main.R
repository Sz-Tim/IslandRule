### Theoretical Ecology
### 2014 Spring
### Tim Szewczyk
### Island body size changes


## w.mean not working when < 5 -- causes a fatal error...

###--- Load libraries ---###
  library(truncnorm)
  library(plyr)
  library(ggplot2)
  library(MASS)
  source(paste0(getwd(), "/BodySize_Funcs.R"))
  
###--- Set parameters ---###

  # Simulations
    sims <- 30
    maxt <- 800
  # Island
    isl.r <- 4      # Number of rows
    isl.c <- 4      # Number of columns
    E.mean <- 100    # Mean initial resource value per cell
    E.sd <- 0    # sd for initial resource value across cells
    prod.mean <- 100    # Productivity mean
    prod.sd <- 0       # Productivity sd
  # Population
    N.init <- 40
    w.mean <- 90   # Initial mean
    w.CV <- 0.01     # sd (initial and for each reproductive event): proportion of w.mean
  # Process parameters
    st <- 0.8   # How does pr(starvation) scale with size? prop: small selects for very small w
    f <- 0.6    # How does fecundity scale with size? Larger favors large w more
    p <- 0.5   # How does pr(predation) scale with size? Larger favors small w more
    prP <- 0.2    # Probability of individual predation; w irrelevant
    predOptimum <- 10
    predEnd <- 800    # Time step to remove predation
    
  # Movement
    m <- 0.5   # Probability of movement

  
###--- Initialize objects ---###
  E0 <- matrix(rnorm(isl.r*isl.c,E.mean,E.sd), nrow=isl.r, ncol=isl.c)
  N0.df <- data.frame(loc.r=sample(1:isl.r, N.init, replace=T),
                     loc.c=sample(1:isl.c, N.init, replace=T),
                     w=rtruncnorm(N.init, a=0, b=Inf, mean=w.mean, sd=w.CV*w.mean))
  mns.postStarve <- matrix(nrow=maxt, ncol=sims)
  mns.postRepro <- matrix(nrow=maxt, ncol=sims)
  pops.postStarve <- matrix(nrow=maxt, ncol=sims)
  pops.postRepro <- matrix(nrow=maxt, ncol=sims)
  sumw.postStarve <- matrix(nrow=maxt, ncol=sims)
  sumw.postRepro <- matrix(nrow=maxt, ncol=sims)
  
  start <- NULL
  predStop <- NULL
  final <- NULL
  
  timing.df <- data.frame(move=rep(NA,maxt*sims),
                          pred=rep(NA,maxt*sims),
                          comp=rep(NA,maxt*sims),
                          repro=rep(NA,maxt*sims),
                          prod=rep(NA,maxt*sims))
  step <- 1

plot(NA, NA, xlim=c(0, maxt), ylim=c(0, E.mean), xlab="Generation", ylab="w post-starve")
  abline(v=predEnd, lty=2)
###--- Run Simulations ---###
for(s in 1:sims) {
  N.df <- N0.df
  E <- E0
  start <- c(start, N.df$w)
  for(t in 1:maxt) {
    ###--- Movement ---###
      N.df <- randmove(N.df, m, isl.r, isl.c)
    ###--- Predation ---###
      if(t < predEnd) {
        N.df <- predation(N.df=N.df, p=p, prP=prP, predOpt=predOptimum, fn="none") # proportional, antilogit, even, diff, none
      }
    
    ###--- Intraspecific Competition ---###
      feedout <- feed(N.df=N.df, E=E, st=st, fn="antilogit") # proportional, outcompete, antilogit, even
      N.df <- feedout$N.df
      E <- feedout$E
    
    ###--- Store Post-Starve Summaries ---###
      mns.postStarve[t,s] <- mean(N.df$w[N.df$w>0])
      pops.postStarve[t,s] <- sum(N.df$w > 0)
      sumw.postStarve[t,s] <- sum(N.df$w)
    
      points(rep(t, sum(N.df$w > 0)), N.df$w[N.df$w>0], col=rgb(0,0,0,0.01), cex=0.1)
    
    ###--- Reproduction ---###
      N.df <- reproduce(N.df=N.df, w.CV=w.CV, f=f, fn="log") # proportional, log, even
    
    ###--- Landscape Regeneration ---###
      E <- grow(E=E, MN=prod.mean, SD=prod.sd)

    ###--- Store Post-Reproduction Summaries ---###
      mns.postRepro[t,s] <- mean(N.df$w)
      pops.postRepro[t,s] <- nrow(N.df)
      sumw.postRepro[t,s] <- sum(N.df$w)
      if(t == predEnd) {
        predStop <- c(predStop, N.df$w)
      }
#        hist(N.df$w, xlim=c(0, 35), main=paste0("t=",t,"; s=",s))
#        Sys.sleep(0.05)

    ###--- Progress Update ---###
    if(t%%25==0) {
     print(paste0("Simulation ",s," Time ",t))
     
    }
    step <- step + 1
  }  # End time loop
  final <- c(final, N.df$w)
}  # End sim loop
  
dist.df <- data.frame(w=c(start, predStop, final), 
                      time=rep(c("Start","Predation Ends", "Final"), 
                               times=c(length(start), length(predStop), length(final))))

###--- Mean w ---### 
  matplot(1:maxt, mns.postStarve, type="l", col=rgb(0,0,0,0.5), lty=1, 
          ylim=c(0, max(mns.postStarve, na.rm=T)), xlab="Generation", ylab="Mean w post-starve")
    abline(v=predEnd, lty=2)
    lines(1:maxt, rowMeans(mns.postStarve), col="green", lwd=2)
  matplot(1:maxt, mns.postRepro, type="l", col=rgb(0,0,0,0.5), lty=1, 
            ylim=c(0, max(mns.postRepro)), xlab="Generation", ylab="Mean w post-reproduction")
    abline(v=predEnd, lty=2)
    lines(1:maxt, rowMeans(mns.postRepro), col="green", lwd=2)

###--- Population Size ---###  
  matplot(1:maxt, pops.postStarve, type="l", col=rgb(0,0,0,0.5), lty=1, 
          ylim=c(0,max(pops.postStarve)), xlab="Generation", ylab="N post-starve")
    abline(v=predEnd, lty=2)
    lines(1:maxt, rowMeans(pops.postStarve), col="green", lwd=2)
  matplot(1:maxt, pops.postRepro, type="l", col=rgb(0,0,0,0.5), lty=1, 
          ylim=c(0,max(pops.postRepro)), xlab="Generation", ylab="N post-reproduction")
    abline(v=predEnd, lty=2)
    lines(1:maxt, rowMeans(pops.postRepro), col="green", lwd=2)


###--- Population Sum of w ---###
  matplot(1:maxt, sumw.postStarve, type="l", col=rgb(0,0,0,0.5), lty=1, 
          ylim=c(0,max(sumw.postStarve)), xlab="Generation", ylab="Sum of w post-starve")
    abline(v=predEnd, lty=2)
    lines(1:maxt, rowMeans(sumw.postStarve), col="green", lwd=2)
  matplot(1:maxt, sumw.postRepro, type="l", col=rgb(0,0,0,0.5), lty=1, 
          ylim=c(0,max(sumw.postRepro)), xlab="Generation", ylab="Sum of w post-reproduction")
    abline(v=predEnd, lty=2)
    lines(1:maxt, rowMeans(sumw.postRepro), col="green", lwd=2)

ggplot(dist.df, aes(x=w, colour=time)) + geom_density(size=1) + xlim(0,max(dist.df$w)) + 
    theme_bw()
  