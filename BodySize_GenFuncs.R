### Theoretical Ecology
### 2014 Spring
### Tim Szewczyk
### Island body size changes



###############################
###--- General functions ---###
###############################

  logit <- function(x) {log(x/(1-x))}

  antilogit <- function(x) {exp(x)/(1+exp(x))}

  # sample(x) uses a vector of elements or, if given an integer, 1:x. To avoid that:
  resample <- function(x,...) {x[sample.int(length(x),...)]}




##################################
###--- Assigning parameters ---###
##################################

  makeParSet <- function(param, low, high, len, logSeq=FALSE, sims, maxt) {
    if(logSeq) { 
      par.seq <- exp( seq(log(low), log(high), length.out=len) )
    } else {
      par.seq <- seq(low, high, length.out=len)
    }
    finalmeans.df <- data.frame(p=rep(par.seq, each=sims))
    names(finalmeans.df) <- param
    meansByTime.df <- data.frame(time=rep(1:maxt, length(par.seq)), p=rep(par.seq, each=maxt))
    names(meansByTime.df)[2] <- param
    return(list(par.seq=par.seq, finalmeans.df=finalmeans.df, meansByTime.df=meansByTime.df, param=param))
  }


  assignPar <- function(param, value) {
    if( sum(ls(island.pars)==param) > 0 ) {
      island.pars[names(island.pars)==param] <<- value
    } else if( sum(ls(pop.pars)==param) > 0 ) {
      pop.pars[names(pop.pars)==param] <<- value
    } else if(sum(ls(move.pars)==param) > 0 ) {
      move.pars[names(move.pars)==param] <<- value
    } else if(sum(ls(pred.pars)==param) > 0 ) {
      pred.pars[names(pred.pars)==param] <<- value
    } else if(sum(ls(feed.pars)==param) > 0 ) {
      feed.pars[names(feed.pars)==param] <<- value
    } else if(sum(ls(repro.pars)==param) > 0 ) {
      repro.pars[names(repro.pars)==param] <<- value
    }
  }




########################################
###--- Saving data and parameters ---###
########################################

  writeDataAndParsToFile <- function(dataframe, param, dirNum, island.pars=island.pars, pop.pars=pop.pars, 
                                     move.pars=move.pars, pred.pars=pred.pars, feed.pars=feed.pars, repro.pars=repro.pars) {
	
	# Set the working directory
		setwd("/Users/fossegrimen/Dropbox/Theoretical Ecology/Project/IslandRule")
		basedir <- paste0(getwd(), "/SimOutput/", feed.pars$feed.fn, "/")
    rundir <- paste0(param, "_ParSet_", as.character(dirNum))
    dir.create(paste0(basedir,rundir))
    
  # Generate filename for dataframe
    fname <- paste0(basedir, rundir, "/SummaryOut.txt")
    
  # Write dataframe
    write.table(dataframe, file=fname, row.names=FALSE)
	
	# Generate filename for parameters
		fname <- paste0(basedir, rundir, "/parameters.R")
	
	# Save parameters into the script
		sink(fname)
  
      # Parameters
			cat("island.pars <- list(isl.r=", island.pars$isl.r, ", ",
									"isl.c=", island.pars$isl.c, ", ",
									"E.mean=", island.pars$E.mean, ", ",
									"E.sd=", island.pars$E.sd, ", ",
									"prod.mean=", island.pars$prod.mean, ", ",
									"prod.sd=", island.pars$prod.sd, ")",
									"\n", sep="")
			cat("pop.pars <- list(N.init=", pop.pars$N.init, ", ",
								 "w.mean=", pop.pars$w.mean, ", ",
								 "w.CV=", pop.pars$w.CV, ")",
								 "\n", sep="")
			cat("move.pars <- list(m=", move.pars$m, ", ",
								  "move.fn='", move.pars$move.fn, "')",
								  "\n", sep="")
			cat("pred.pars <- list(p=", pred.pars$p, ", ",
								  "prP=", pred.pars$prP, ", ",
								  "predOpt=", pred.pars$predOpt, ", ",
								  "predEnd=", pred.pars$predEnd, ", ",
								  "pred.fn='", pred.pars$pred.fn, "')",
								  "\n", sep="")
			cat("feed.pars <- list(st=", feed.pars$st, ", ",
								  "feed.fn='", feed.pars$feed.fn, "')",
								  "\n", sep="")
			cat("repro.pars <- list(f=", repro.pars$f, ", ",
								   "babybump=", repro.pars$babybump, ", ",
								   "repro.fn='", repro.pars$repro.fn, "')",
								   "\n", sep="")
			
      # Summary statistics
			final.df <- subset(sim.out$summary.df, 
							   sim.out$summary.df$Time==max(sim.out$summary.df$Time))
			init.df <- subset(sim.out$summary.df, 
							   sim.out$summary.df$Time==min(sim.out$summary.df$Time))
			cat("nsims <- ", nrow(final.df), "\n", sep="")
			cat("MeanW.init <- ", mean(init.df$MeanW.pF), "\n", sep="")
			cat("MeanW.pF <- ", mean(final.df$MeanW.pF), "\n", sep="")
			cat("MeanW.pR <- ", mean(final.df$MeanW.pR), "\n", sep="")
			cat("SumW.init <- ", mean(init.df$SumW.pF), "\n", sep="")
			cat("SumW.pF <- ", mean(final.df$SumW.pF), "\n", sep="")
			cat("SumW.pR <- ", mean(final.df$SumW.pR), "\n", sep="")
			cat("N.init <- ", mean(init.df$N.pF), "\n", sep="")
			cat("N.pF <- ", mean(final.df$N.pF), "\n", sep="")
			cat("N.pR <- ", mean(final.df$N.pR), "\n", sep="")
		sink()
			
  }



  writeDistrToFile <- function(distr.ls, param, dirNum) {
    # Set directories
      basedir <- paste0(getwd(), "/SimOutput/", feed.pars$feed.fn, "/")
      rundir <- paste0(param, "_ParSet_", as.character(dirNum))    
    # Write starting and final distributions
      write.matrix(distr.ls[[1]], file=paste0(basedir, rundir, "/initDistr.txt"))
      write.matrix(distr.ls[[2]], file=paste0(basedir, rundir, "/finalDistr.txt"))
  }




########################################
###--- Retrieving data from files ---###
########################################

  getSummary <- function(feedFunc, param, indices) {
    # Set up directories
      basedir <- paste0(getwd(), "/SimOutput/", feedFunc, "/")
    # Number of rows
      n <- length(indices)
    # Set up dataframe and vectors
      df <- data.frame(index=numeric(n), 
                       isl.r=numeric(n), 
                       isl.c=numeric(n),
                       E.mean=numeric(n),
                       E.sd=numeric(n),
                       prod.mean=numeric(n),
                       prod.sd=numeric(n),
                       N.init=numeric(n),
                       w.mean=numeric(n),
                       w.CV=numeric(n),
                       m=numeric(n),
                       p=numeric(n),
                       prP=numeric(n),
                       predOpt=numeric(n),
                       predEnd=numeric(n),
                       st=numeric(n),
                       f=numeric(n),
                       babybump=numeric(n),
                       MeanW.init=numeric(n),
                       MeanW.pF=numeric(n),
                       MeanW.pR=numeric(n),
                       SumW.init=numeric(n),
                       SumW.pF=numeric(n),
                       SumW.pR=numeric(n),
                       N.pF=numeric(n),
                       N.pR=numeric(n))
      movefn <- character(n)
      predfn <- character(n)
      feedfn <- character(n)
      reprofn <- character(n)
    
    # Load data into objects
      for(i in 1:n) {
        ri <- indices[i]
        source(paste0(basedir, param, "_ParSet_", as.character(i), "/parameters.R"))

        df[i,1] <- ri
        df[i,2] <- island.pars$isl.r
        df[i,3] <- island.pars$isl.c
        df[i,4] <- island.pars$E.mean
        df[i,5] <- island.pars$E.sd
        df[i,6] <- island.pars$prod.mean
        df[i,7] <- island.pars$prod.sd
        df[i,8] <- pop.pars$N.init
        df[i,9] <- pop.pars$w.mean
        df[i,10] <- pop.pars$w.CV
        df[i,11] <- move.pars$m
        df[i,12] <- pred.pars$p
        df[i,13] <- pred.pars$prP
        df[i,14] <- pred.pars$predOpt
        df[i,15] <- pred.pars$predEnd
        df[i,16] <- feed.pars$st
        df[i,17] <- repro.pars$f
        df[i,18] <- repro.pars$babybump
        df[i,19] <- MeanW.init
        df[i,20] <- MeanW.pF
        df[i,21] <- MeanW.pR
        df[i,22] <- SumW.init
        df[i,23] <- SumW.pF
        df[i,24] <- SumW.pR
        df[i,25] <- N.pF
        df[i,26] <- N.pR
        movefn[i] <- move.pars$move.fn
        predfn[i] <- pred.pars$pred.fn
        feedfn[i] <- feed.pars$feed.fn
        reprofn[i] <- repro.pars$repro.fn
      }
      df$move.fn <- movefn
      df$pred.fn <- predfn
      df$feed.fn <- feedfn
      df$repro.fn <- reprofn
    return(df)
  }



  getSimsByTime <- function(feedFunc, param, indices) {

    # Set up directories
      basedir <- paste0(getwd(), "/SimOutput/", feedFunc, "/")
    # Number of rows
      n <- length(indices)
    # Set up dataframe
     df <- data.frame(index=numeric(n), 
                     isl.r=numeric(n), 
                     isl.c=numeric(n),
                     E.mean=numeric(n),
                     E.sd=numeric(n),
                     prod.mean=numeric(n),
                     prod.sd=numeric(n),
                     N.init=numeric(n),
                     w.mean=numeric(n),
                     w.CV=numeric(n),
                     m=numeric(n),
                     p=numeric(n),
                     prP=numeric(n),
                     predOpt=numeric(n),
                     predEnd=numeric(n),
                     st=numeric(n),
                     f=numeric(n),
                     babybump=numeric(n),
                     time=numeric(n),
                     MeanW.pF=numeric(n),
                     MeanW.pR=numeric(n),
                     SumW.pF=numeric(n),
                     SumW.pR=numeric(n),
                     N.pF=numeric(n),
                     N.pR=numeric(n))
    
    for(i in 1:n) {
      ri <- indices[i]
      # Load data
          simOut <- read.table(paste0(basedir, param, "_ParSet_", as.character(i), "/SummaryOut.txt"), header=T)
          source(paste0(basedir, param, "_ParSet_", as.character(i), "/parameters.R"))
        
        # Calculate means
          i.df <- ddply(simOut, .(Time), summarize, mnW.pF=mean(MeanW.pF), mnW.pR=mean(MeanW.pR),
                        sumW.pF=mean(SumW.pF), sumW.pR=mean(SumW.pR), n.pF=mean(N.pF), n.pR=mean(N.pR))
        # Fill database with parameters
          maxt <- max(i.df$Time)
          dfrows <- (1 + maxt*(i-1)):(maxt*i)
          df[dfrows,1] <- ri
          df[dfrows,2] <- island.pars$isl.r
          df[dfrows,3] <- island.pars$isl.c
          df[dfrows,4] <- island.pars$E.mean
          df[dfrows,5] <- island.pars$E.sd
          df[dfrows,6] <- island.pars$prod.mean
          df[dfrows,7] <- island.pars$prod.sd
          df[dfrows,8] <- pop.pars$N.init
          df[dfrows,9] <- pop.pars$w.mean
          df[dfrows,10] <- pop.pars$w.CV
          df[dfrows,11] <- move.pars$m
          df[dfrows,12] <- pred.pars$p
          df[dfrows,13] <- pred.pars$prP
          df[dfrows,14] <- pred.pars$predOpt
          df[dfrows,15] <- pred.pars$predEnd
          df[dfrows,16] <- feed.pars$st
          df[dfrows,17] <- repro.pars$f
          df[dfrows,18] <- repro.pars$babybump
          df[dfrows,19] <- i.df$Time
          df[dfrows,20] <- i.df$mnW.pF
          df[dfrows,21] <- i.df$mnW.pR
          df[dfrows,22] <- i.df$sumW.pF
          df[dfrows,23] <- i.df$sumW.pR
          df[dfrows,24] <- i.df$n.pF
          df[dfrows,25] <- i.df$n.pR
      }
    return(df)
  }














