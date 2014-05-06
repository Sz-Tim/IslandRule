### Theoretical Ecology
### 2014 Spring
### Tim Szewczyk
### Island body size changes

# Requires:
# - BodySize_Funcs.R
# - truncnorm
# - MASS

simIsland <- function(sims=10, maxt=100, island.pars, pop.pars, move.pars, pred.pars, 
                      feed.pars, repro.pars, plotIndiv=FALSE, plotAll=FALSE, plotLims=c(0,E.mean), alpha=0.02) {
  
  ###--- UNPACK PARAMETERS ---###
    # Island Parameters
    isl.r <- island.pars$isl.r           # Number of rows in island
    isl.c <- island.pars$isl.c           # Number of columns in island
    E.mean <- island.pars$E.mean         # Mean initial resource level per cell
    E.sd <- island.pars$E.sd             # sd for initial resource level per cell
    prod.mean <- island.pars$prod.mean   # Mean productivity per cell per time step
    prod.sd <- island.pars$prod.sd       # sd for productivity per cell per time step
    incComp <- island.pars$incComp       # Mean productivity after competition is increased
    incComp.time <- island.pars$incComp.time  # Time step to increase competition
    
    # Population Parameters
    N.init <- pop.pars$N.init            # Initial population size
    w.mean <- pop.pars$w.mean            # Population mean w
    w.CV <- pop.pars$w.CV                # w variance: sd = w.mean*w.CV; initial population and offspring vs parent
    
    # Movement Parameters
    m <- move.pars$m                     # Individual probability of movement
    move.fn <- move.pars$move.fn         # Which movement function to use
    
    # Predation Parameters
    p <- pred.pars$p                     # How predation risk scales with w; larger p favors small w
    prP <- pred.pars$prP                 # Individual probability of predation (pred.fn=="even" only)
    predOpt <- pred.pars$predOpt         # Optimal w -- highest predation risk (pred.fn=="diff" only)
    predEnd <- pred.pars$predEnd         # Time step where predation ends
    pred.fn <- pred.pars$pred.fn         # Which predation probability function to use
    
    # Feeding Parameters
    st <- feed.pars$st                   # How starvation probability scales with w; how it works depends on feed.fn
    feed.fn <- feed.pars$feed.fn         # Which starvation probability function to use
    
    # Reproduction Parameters
    f <- repro.pars$f                    # How reproductive rate scales with w; how it works depends on repro.fn
    babybump <- repro.pars$babybump      # Increase lambda for all individuals; necessary for w < ~5 to avoid fatal errors
    repro.fn <- repro.pars$repro.fn      # Which reproduction function to use to calculate lambda for each individual
    
    
    

  ###--- INITIALIZE OBJECTS ---###
    # Generate the island
    E0 <- matrix(rnorm(isl.r*isl.c, E.mean, E.sd), nrow=isl.r, ncol=isl.c)
    # Generate the population
    N0.df <- data.frame(loc.r=sample(1:isl.r, N.init, replace=T),
                    loc.c=sample(1:isl.c, N.init, replace=T),
                    w=rtruncnorm(N.init, a=0, b=Inf, mean=w.mean, sd=w.CV*w.mean))
    # Object for storing summary statistics
    summary.df <- data.frame(Sim=rep(1:sims, each=maxt),
                    Time=rep(1:maxt, sims),
                    MeanW.pF=numeric(sims*maxt),
                    VarW.pF=numeric(sims*maxt),
                    MeanW.pR=numeric(sims*maxt),
                    SumW.pF=numeric(sims*maxt),
                    SumW.pR=numeric(sims*maxt),
                    N.pF=numeric(sims*maxt),
                    N.pR=numeric(sims*maxt))
    # Object for storing initial and final w's
    distr.ls <- list(start=NULL, final=NULL)
    finalmeans <- numeric(sims)
    # Increases with each time step across all simulations (i.e., from 1 to sims*maxt)
    step <- 1
  
  
  
  
  ###--- PLOT FOR ALL INDIVIDUALS ---###
    if(plotIndiv & !plotAll) {
      plot(NA, NA, xlim=c(0,maxt), ylim=plotLims, xlab="Generation", ylab="w")
      if(pred.fn != "none") {
        abline(v=predEnd, lty=2, col="darkred")
      }
      abline(v=incComp.time, lty=2, col="darkblue")
    }
  
  
  
    
  ###--- RUN SIMULATIONS ---###
    for(s in 1:sims) {
      
      ## Initialize population and landscape ##
        N.df <- N0.df
        E <- E0
        distr.ls$start <- c(distr.ls$start, N.df$w)
        
      for(t in 1:maxt) {
        
        ## Movement ##
          N.df <- move(N.df, m, isl.r, isl.c, fn="random")
        
        ## Predation ##
          if((pred.fn != "none") & (t < predEnd)) {
            N.df <- predation(N.df=N.df, p=p, prP=prP, predOpt=predOpt, fn=pred.fn)
          }
          
        ## Feeding ##
          feedout <- feed(N.df=N.df, E=E, st=st, fn=feed.fn)
          N.df <- feedout$N.df
          E <- feedout$E
        
        ## Store post-feed summary statistics, plot all individuals ##
          summary.df[step, c(3, 4, 6, 8)] <- c(mean(N.df$w), var(N.df$w), sum(N.df$w), nrow(N.df))
          if(plotIndiv) {
            if(nrow(N.df) > 100) {
              N.plot <- N.df[sample(1:nrow(N.df), 100, replace=FALSE),]
            } else {
              N.plot <- N.df
            }
            points(rep(t, nrow(N.plot)), N.plot$w, col=rgb(0,0,0,alpha), cex=0.2)
          }
        
        ## Reproduction ##
          N.df <- reproduce(N.df=N.df, w.CV=w.CV, f=f, babybump=babybump, fn=repro.fn)
        
        ## Store post-reproduction summary statistics ##
          summary.df[step, c(5, 7, 9)] <- c(mean(N.df$w), sum(N.df$w), nrow(N.df))
        
        ## Landscape regeneration ##
          E <- grow(E=E, MN=prod.mean, SD=prod.sd, incComp=incComp, incComp.time=incComp.time, t=t)
        
        ## Progress Update ##
          if(t%%25 == 0) {
            cat("Simulation", s, "Time", t, "\n")
          }
          step <- step + 1
        
      } # Close time loop
      
      ## Store final w distribution ##
        distr.ls$final <- c(distr.ls$final, N.df$w)
        finalmeans[s] <- mean(N.df$w)
      
    } # Close simulation loop
  
  return(list(N.df=N.df, E=E, summary.df=summary.df, distr.ls=distr.ls, finalmeans=finalmeans))
} # Close function loop

