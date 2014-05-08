### Theoretical Ecology
### 2014 Spring
### Tim Szewczyk
### Island body size changes


####################
###--- Set Up ---###
####################
  setwd("~/Dropbox/Theoretical Ecology/Project/IslandRule")
  baseDir <- getwd()
  library(truncnorm)
  library(plyr)
  library(ggplot2)
  library(MASS)
  theme_set(theme_bw())
  source(paste0(baseDir, "/BodySize_Funcs.R"))
  source(paste0(baseDir, "/BodySize_GenFuncs.R"))
  source(paste0(baseDir, "/BodySize_SimFunc.R"))
  

############################  
###--- Set parameters ---###
############################

    sims <- 50
    maxt <- 500
    parLen <- 8
    plotAllSims <- FALSE
    writeSims <- TRUE

    island.pars <- list(isl.r=5,            # Number of rows in island
                        isl.c=5,            # Number of columns in island
                        E.mean=20,         # Mean initial resource level per cell
                        E.sd=0,             # sd for initial resource level per cell
                        prod.mean=20,      # Mean productivity per cell per time step
                        prod.sd=0,          # sd for productivity per cell per time step
                        incComp=10,        # Mean productivity after competition is increased
                        incComp.time=500)   # Time to increase competition
 
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
                      feed.fn="proportional")   # Starvation probability function: antilogit, proportional, outcompete, even

    repro.pars <- list(f=1,                 # How reproductive rate scales with w; how it works depends on repro.fn
                       babybump=1,          # Increase lambda for all individuals; necessary for w < ~5 to avoid fatal errors
                       repro.fn="log")      # Reproduction function to calculate individual lambda: proportional, log, even
  	mainland <- FALSE
  
  
#############################  
###--- Run Simulations ---###
#############################

  # Select parameter to vary across simulation sets:
  parSet <- makeParSet(param="f", low=0.5, high=3, len=parLen, logSeq=FALSE, sims=sims, maxt=maxt)
  
  # Run simulations
  dirNum <- 1
  if(plotAllSims) {
    plot(NA, NA, xlim=c(0,maxt), ylim=c(0,island.pars$E.mean), xlab="Generation", ylab="w")
  }
  
  for(i in parSet$par.seq) {

    # Set parameters
      assignPar(parSet$param, i)

    # Run simulations  
      sim.out <- simIsland(sims, maxt, island.pars, pop.pars, move.pars, pred.pars, feed.pars, 
                           repro.pars, plotIndiv=FALSE, plotAll=plotAllSims, plotLims=c(0,25))
    
    # Write data to files
      if(writeSims) {
        writeDataAndParsToFile(sim.out$summary.df, parSet$param, dirNum, island.pars, pop.pars, move.pars, pred.pars, feed.pars, repro.pars, mainland)
        writeDistrToFile(sim.out$distr.ls, parSet$param, dirNum, mainland)
        dirNum <- dirNum + 1
      }
    cat("Finished par set",i,"\n \n")                                  
  }

    
# ####################
# ###--- GRAPHS ---###
# ####################
# 
#     summary.df <- getSummary(feedFunc=feed.pars$feed.fn, param=parSet$param, indices=1:parLen, mainland)
#     byTime.df <- getSimsByTime(feedFunc=feed.pars$feed.fn, param=parSet$param, indices=1:parLen, mainland)
#     
#   ###--- across a parameter ---###
#     ggplot(summary.df, aes(x=w.mean, y=MeanW.pF)) + geom_point(size=3, alpha=0.5) + ylim(0,20) + labs(x="Initial w", y="Final mean w", title="w")
#     ggplot(byTime.df, aes(x=time, y=MeanW.pF, group=w.mean, colour=w.mean)) + geom_line() + ylim(0,20)
#     ggplot(summary.df, aes(x=MeanW.init, y=(MeanW.pF-MeanW.init))) + geom_line() + geom_hline(y=0, linetype=2)
# 
#   ###--- mean w ---### 
#     ggplot(sim.out$summary.df, aes(x=Time, group=Sim)) + geom_line(aes(y=MeanW.pF)) + geom_line(aes(y=MeanW.pR), colour="red")
# 
#   ###--- population size ---###  
#     ggplot(sim.out$summary.df, aes(x=Time, group=Sim))  + geom_line(aes(y=N.pF)) + geom_line(aes(y=N.pR), colour="red")
# 
#   ###--- population sum of w ---###
#     ggplot(sim.out$summary.df, aes(x=Time, group=Sim))  + geom_line(aes(y=SumW.pF)) + geom_line(aes(y=SumW.pR), colour="red")
# 
#   ###--- starting and final distributions ---###
#     ggplot(distr.df, aes(x=w, colour=time)) + geom_density(size=1) + xlim(0,max(distr.df$w)) + theme_bw()
#   