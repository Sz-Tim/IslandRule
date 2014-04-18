### Theoretical Ecology
### 2014 Spring
### Tim Szewczyk
### Island body size changes


####################
###--- Set Up ---###
####################
  setwd("/Users/fossegrimen/Dropbox/Theoretical Ecology/Project/IslandRule")
  library(truncnorm)
  library(plyr)
  library(ggplot2)
  library(MASS)
  theme_set(theme_bw())
  source(paste0(getwd(), "/BodySize_Funcs.R"))
  source(paste0(getwd(), "/BodySize_SimFunc.R"))
  

############################  
###--- Set parameters ---###
############################

    sims <- 50
    maxt <- 500
    plotAllSims <- FALSE

    island.pars <- list(isl.r=5,            # Number of rows in island
                        isl.c=5,            # Number of columns in island
                        E.mean=20,         # Mean initial resource level per cell
                        E.sd=0,             # sd for initial resource level per cell
                        prod.mean=20,      # Mean productivity per cell per time step
                        prod.sd=0)          # sd for productivity per cell per time step
 
    pop.pars <- list(N.init=40,             # Initial population size
                     w.mean=0.5,             # Population mean w
                     w.CV=0.04)             # w variance: sd = w.mean*w.CV; initial population and offspring vs parent

    move.pars <- list(m=0.5,                # Individual probability of movement
                      move.fn="random")     # Which movement function to use

    pred.pars <- list(p=0.5,                # How predation risk scales with w; larger p favors small w
                      prP=0.2,              # Individual probability of predation (pred.fn=="even" only)
                      predOpt=10,           # Optimal w -- highest predation risk (pred.fn=="diff" only)
                      predEnd=300,          # Time step where predation ends
                      pred.fn="none")       # Predation probability function: none, antilogit, proportional, even, diff
  
    feed.pars <- list(st=0.3,                   # How starvation probability scales with w; how it works depends on feed.fn; max 0.6 for outcompete
                      feed.fn="antilogit")   # Starvation probability function: antilogit, proportional, outcompete, even

    repro.pars <- list(f=1,                 # How reproductive rate scales with w; how it works depends on repro.fn
                       babybump=1,          # Increase lambda for all individuals; necessary for w < ~5 to avoid fatal errors
                       repro.fn="log")      # Reproduction function to calculate individual lambda: proportional, log, even
  
  
#############################  
###--- Run Simulations ---###
#############################

  parSet <- makeParSet(par="w.mean", low=0.1, high=18, len=8, logSeq=TRUE, sims=sims, maxt=maxt)
  
  finalmeans.df <- parSet$finalmeans.df;    finalmeans <- NULL
  meansByTime.df <- parSet$meansByTime.df;     meanWs <- NULL
  if(plotAllSims) {
    plot(NA, NA, xlim=c(0,maxt), ylim=c(0,island.pars$E.mean), xlab="Generation", ylab="w")
  }
  for(i in parSet$par.seq) {
    assignPar(parSet$par, i)
    sim.out <- simIsland(sims, maxt, island.pars, pop.pars, move.pars, pred.pars, feed.pars, 
                         repro.pars, plotIndiv=FALSE, plotAll=plotAllSims)
    distr.df <- data.frame(w=c(sim.out$distr.ls[[1]], sim.out$distr.ls[[2]]),
                           time=c( rep("Initial", length(sim.out$distr.ls[[1]])),
                                   rep("Final", length(sim.out$distr.ls[[2]])) ))
    finalmeans <- c(finalmeans, sim.out$finalmeans)
    meanWs <- c(meanWs, ddply(sim.out$summary.df, .(Time), summarize, MeanW=mean(MeanW.pF))[,2])
    cat("Finished par set",i,"\n \n")                                  
  }
  finalmeans.df$finalmeans <- finalmeans
  meansByTime.df$meanW <- meanWs
    
####################
###--- GRAPHS ---###
####################  
  ###--- across a parameter ---###
    ggplot(finalmeans.df, aes(x=w.mean, y=finalmeans)) + geom_point(size=3, alpha=0.5) + ylim(0,20) + labs(x="Initial w", y="Final mean w", title="w")
    ggplot(meansByTime.df, aes(x=time, y=meanWs, group=w.mean, colour=w.mean)) + geom_line() + ylim(0,20)

  ###--- mean w ---### 
    ggplot(sim.out$summary.df, aes(x=Time, group=Sim)) + geom_line(aes(y=MeanW.pF)) + geom_line(aes(y=MeanW.pR), colour="red")

  ###--- population size ---###  
    ggplot(sim.out$summary.df, aes(x=Time, group=Sim))  + geom_line(aes(y=N.pF)) + geom_line(aes(y=N.pR), colour="red")

  ###--- population sum of w ---###
    ggplot(sim.out$summary.df, aes(x=Time, group=Sim))  + geom_line(aes(y=SumW.pF)) + geom_line(aes(y=SumW.pR), colour="red")

  ###--- starting and final distributions ---###
    ggplot(distr.df, aes(x=w, colour=time)) + geom_density(size=1) + xlim(0,max(distr.df$w)) + theme_bw()
  