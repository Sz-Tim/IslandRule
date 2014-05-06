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
                                     move.pars=move.pars, pred.pars=pred.pars, feed.pars=feed.pars, repro.pars=repro.pars, mainland=FALSE) {
	
	# Working directory should contain all scripts and a folder called "SimOutput" that contains a folder for each feeding function
		basedir <- paste0(getwd(), "/SimOutput/", feed.pars$feed.fn, "/")
	if(mainland) {
		rundir <- paste0(param, "_MainlandParSet_", as.character(dirNum))
	} else {	
    	rundir <- paste0(param, "_ParSet_", as.character(dirNum))
    }
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
									"prod.sd=", island.pars$prod.sd, ", ",
                  "incComp=", island.pars$incComp, ", ",
                  "incComp.time=", island.pars$incComp.time, ")",
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
      cat("VarW.pF <- ", mean(final.df$VarW.pF), "\n", sep="")
			cat("MeanW.pR <- ", mean(final.df$MeanW.pR), "\n", sep="")
			cat("SumW.init <- ", mean(init.df$SumW.pF), "\n", sep="")
			cat("SumW.pF <- ", mean(final.df$SumW.pF), "\n", sep="")
			cat("SumW.pR <- ", mean(final.df$SumW.pR), "\n", sep="")
			cat("N.init <- ", mean(init.df$N.pF), "\n", sep="")
			cat("N.pF <- ", mean(final.df$N.pF), "\n", sep="")
			cat("N.pR <- ", mean(final.df$N.pR), "\n", sep="")
		sink()
			
  }



  writeDistrToFile <- function(distr.ls, param, dirNum, mainland=FALSE) {
    # Set directories
      basedir <- paste0(getwd(), "/SimOutput/", feed.pars$feed.fn, "/")
      	if(mainland) {
		rundir <- paste0(param, "_MainlandParSet_", as.character(dirNum))
	} else {	
    	rundir <- paste0(param, "_ParSet_", as.character(dirNum))
    }  
    # Write starting and final distributions
      write.matrix(distr.ls[[1]], file=paste0(basedir, rundir, "/initDistr.txt"))
      write.matrix(distr.ls[[2]], file=paste0(basedir, rundir, "/finalDistr.txt"))
  }




########################################
###--- Retrieving data from files ---###
########################################

  getSummary <- function(feedFunc, param, indices, mainland=FALSE) {
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
                       incComp=numeric(n),
                       incComp.time=numeric(n),
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
                       VarW.pF=numeric(n),
                       MeanW.pF=numeric(n),
                       VarW.pF=numeric(n),
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
        if(mainland) {
	        source(paste0(basedir, param, "_MainlandParSet_", as.character(i), "/parameters.R"))
		} else {
			source(paste0(basedir, param, "_ParSet_", as.character(i), "/parameters.R"))
		}
        df[i,1] <- ri
        df[i,2] <- island.pars$isl.r
        df[i,3] <- island.pars$isl.c
        df[i,4] <- island.pars$E.mean
        df[i,5] <- island.pars$E.sd
        df[i,6] <- island.pars$prod.mean
        df[i,7] <- island.pars$prod.sd
        df[i,8] <- island.pars$incComp
        df[i,9] <- island.pars$incComp.time
        df[i,10] <- pop.pars$N.init
        df[i,11] <- pop.pars$w.mean
        df[i,12] <- pop.pars$w.CV
        df[i,13] <- move.pars$m
        df[i,14] <- pred.pars$p
        df[i,15] <- pred.pars$prP
        df[i,16] <- pred.pars$predOpt
        df[i,17] <- pred.pars$predEnd
        df[i,18] <- feed.pars$st
        df[i,19] <- repro.pars$f
        df[i,20] <- repro.pars$babybump
        df[i,21] <- MeanW.init
        df[i,22] <- MeanW.pF
        df[i,23] <- Var.pF
        df[i,24] <- MeanW.pR
        df[i,25] <- SumW.init
        df[i,26] <- SumW.pF
        df[i,27] <- SumW.pR
        df[i,28] <- N.pF
        df[i,29] <- N.pR
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



  getSimsByTime <- function(feedFunc, param, indices, mainland=FALSE) {

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
                     incComp=numeric(n),
                     incComp.time=numeric(n),
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
                     VarW.pF=numeric(n),
                     MeanW.pR=numeric(n),
                     SumW.pF=numeric(n),
                     SumW.pR=numeric(n),
                     N.pF=numeric(n),
                     N.pR=numeric(n))
    
    for(i in 1:n) {
      ri <- indices[i]
      # Load data
      	if(mainland) {
      		simOut <- read.table(paste0(basedir, param, "_MainlandParSet_", as.character(i), "/SummaryOut.txt"), header=T)
	        source(paste0(basedir, param, "_MainlandParSet_", as.character(i), "/parameters.R"))
		} else {
			simOut <- read.table(paste0(basedir, param, "_ParSet_", as.character(i), "/SummaryOut.txt"), header=TRUE)
			source(paste0(basedir, param, "_ParSet_", as.character(i), "/parameters.R"))
		}
        
        # Calculate means
          i.df <- ddply(simOut, .(Time), summarize, mnW.pF=mean(MeanW.pF), varW.pF=mean(VarW.pF), mnW.pR=mean(MeanW.pR),
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
          df[dfrows,8] <- island.pars$incComp
          df[dfrows,9] <- island.pars$incComp.time
          df[dfrows,10] <- pop.pars$N.init
          df[dfrows,11] <- pop.pars$w.mean
          df[dfrows,12] <- pop.pars$w.CV
          df[dfrows,13] <- move.pars$m
          df[dfrows,14] <- pred.pars$p
          df[dfrows,15] <- pred.pars$prP
          df[dfrows,16] <- pred.pars$predOpt
          df[dfrows,17] <- pred.pars$predEnd
          df[dfrows,18] <- feed.pars$st
          df[dfrows,19] <- repro.pars$f
          df[dfrows,20] <- repro.pars$babybump
          df[dfrows,21] <- i.df$Time
          df[dfrows,22] <- i.df$mnW.pF
          df[dfrows,23] <- i.df$varW.pF
          df[dfrows,24] <- i.df$mnW.pR
          df[dfrows,25] <- i.df$sumW.pF
          df[dfrows,26] <- i.df$sumW.pR
          df[dfrows,27] <- i.df$n.pF
          df[dfrows,28] <- i.df$n.pR
      }
    return(df)
  }


  summarizeAll <- function(directory) { # Set up for ../IslandRule/SimOutput
    
    # Initialize vectors
      isl.r <- NULL
      isl.c <- NULL
      E.mean <- NULL
      E.sd <- NULL
      prod.mean <- NULL
      prod.sd <- NULL
      incComp <- NULL
      incComp.time <- NULL
      N.init <- NULL
      w.mean <- NULL
      w.CV <- NULL
      m <- NULL
      movefn <- NULL
      p <- NULL
      prP <- NULL
      predOpt <- NULL
      predEnd <- NULL
      predfn <- NULL
      st <- NULL
      feedfn <- NULL
      f <- NULL
      babybump <- NULL
      reprofn <- NULL
      meanW.init <- NULL
      meanW.pF <- NULL
      varW.pF <- NULL
      meanW.pR <- NULL
      sumW.init <- NULL
      sumW.pF <- NULL
      sumW.pR <- NULL
      n.pF <- NULL
      n.pR <- NULL
    
    # Fill vectors
      for(subDir in dir(directory, full.names=TRUE)) {
        for(parDir in dir(subDir, full.names=TRUE)) {
          source(paste0(parDir, "/parameters.R"))
          isl.r <- c(isl.r, island.pars$isl.r)
          isl.c <- c(isl.c, island.pars$isl.c)
          E.mean <- c(E.mean, island.pars$E.mean)
          E.sd <- c(E.sd, island.pars$E.sd)
          prod.mean <- c(prod.mean, island.pars$prod.mean)
          prod.sd <- c(prod.sd, island.pars$prod.sd)
          incComp <- c(incComp, island.pars$incComp)
          incComp.time <- c(incComp.time, island.pars$incComp.time)
          N.init <- c(N.init, pop.pars$N.init)
          w.mean <- c(w.mean, pop.pars$w.mean)
          w.CV <- c(w.CV, pop.pars$w.CV)
          m <- c(m, move.pars$m)
          movefn <- c(movefn, move.pars$move.fn)
          p <- c(p, pred.pars$p)
          prP <- c(prP, pred.pars$prP)
          predOpt <- c(predOpt, pred.pars$predOpt)
          predEnd <- c(predEnd, pred.pars$predEnd)
          predfn <- c(predfn, pred.pars$pred.fn)
          st <- c(st, feed.pars$st)
          feedfn <- c(feedfn, feed.pars$feed.fn)
          f <- c(f, repro.pars$f)
          babybump <- c(babybump, repro.pars$babybump)
          reprofn <- c(reprofn, repro.pars$repro.fn)
          meanW.init <- c(meanW.init, MeanW.init)
          meanW.pF <- c(meanW.pF, MeanW.pF)
          varW.pF <- c(varW.pF, VarW.pF)
          meanW.pR <- c(meanW.pR, MeanW.pR)
          sumW.init <- c(sumW.init, SumW.init)
          sumW.pF <- c(sumW.pF, SumW.pF)
          sumW.pR <- c(sumW.pR, SumW.pR)
          n.pF <- c(n.pF, N.pF)
          n.pR <- c(n.pR, N.pR)
        } # close parDir loop
      } # close subDir loop
    df <- data.frame(isl.r, isl.c, E.mean, E.sd, prod.mean, prod.sd, incComp, incComp.time,
                     N.init, w.mean, w.CV,
                     m, move.fn=movefn,
                     p, prP, predOpt, predEnd, pred.fn=predfn,
                     st, feed.fn=feedfn,
                     f, babybump, repro.fn=reprofn,
                     MeanW.init=meanW.init, MeanW.pF=meanW.pF, VarW.pF=varW.pF, MeanW.pR=meanW.pR,
                     SumW.init=sumW.init, SumW.pF=sumW.pF, SumW.pR=sumW.pR,
                     N.pF=n.pF, N.pR=n.pR)
    return(df)
  }



  getDistr <- function(directory, time="final", sampsize=500) {
    # Initialize vectors
    isl.r <- NULL
    isl.c <- NULL
    E.mean <- NULL
    E.sd <- NULL
    prod.mean <- NULL
    prod.sd <- NULL
    incComp <- NULL
    incComp.time <- NULL
    n.init <- NULL
    w.mean <- NULL
    w.CV <- NULL
    m <- NULL
    movefn <- NULL
    p <- NULL
    prP <- NULL
    predOpt <- NULL
    predEnd <- NULL
    predfn <- NULL
    st <- NULL
    feedfn <- NULL
    f <- NULL
    babybump <- NULL
    reprofn <- NULL
    meanW.init <- NULL
    meanW.pF <- NULL
    varW.pF <- NULL
    meanW.pR <- NULL
    sumW.init <- NULL
    sumW.pF <- NULL
    sumW.pR <- NULL
    n.pF <- NULL
    n.pR <- NULL
    w <- NULL
    
    # Fill vectors
    for(subDir in dir(directory, full.names=TRUE)) {
      for(parDir in dir(subDir, full.names=TRUE)) {
        source(paste0(parDir, "/parameters.R"))
        distr.new <- read.table(paste0(parDir, "/", time, "Distr.txt"))[,1]
        l <- length(distr.new)
        if(l > sampsize) {
          distr.new <- distr.new[sample(1:l, size=sampsize, replace=FALSE)]
          l <- sampsize
        }
        isl.r <- c(isl.r, rep(island.pars$isl.r,l))
        isl.c <- c(isl.c, rep(island.pars$isl.c,l))
        E.mean <- c(E.mean, rep(island.pars$E.mean,l))
        E.sd <- c(E.sd, rep(island.pars$E.sd,l))
        prod.mean <- c(prod.mean, rep(island.pars$prod.mean,l))
        prod.sd <- c(prod.sd, rep(island.pars$prod.sd,l))
        incComp <- c(incComp, rep(island.pars$incComp,l))
        incComp.time <- c(incComp, rep(island.pars$incComp.time,l))
        n.init <- c(n.init, rep(pop.pars$N.init,l))
        w.mean <- c(w.mean, rep(pop.pars$w.mean,l))
        w.CV <- c(w.CV, rep(pop.pars$w.CV,l))
        m <- c(m, rep(move.pars$m,l))
        movefn <- c(movefn, rep(as.character(move.pars$move.fn),l))
        p <- c(p, rep(pred.pars$p,l))
        prP <- c(prP, rep(pred.pars$prP,l))
        predOpt <- c(predOpt, rep(pred.pars$predOpt,l))
        predEnd <- c(predEnd, rep(pred.pars$predEnd,l))
        predfn <- c(predfn, rep(as.character(pred.pars$pred.fn),l))
        st <- c(st, rep(feed.pars$st,l))
        feedfn <- c(feedfn, rep(as.character(feed.pars$feed.fn),l))
        f <- c(f, rep(repro.pars$f,l))
        babybump <- c(babybump, rep(repro.pars$babybump,l))
        reprofn <- c(reprofn, rep(as.character(repro.pars$repro.fn),l))
        meanW.init <- c(meanW.init, rep(MeanW.init,l))
        meanW.pF <- c(meanW.pF, rep(MeanW.pF,l))
        varW.pF <- c(varW.pF, rep(VarW.pF,l))
        meanW.pR <- c(meanW.pR, rep(MeanW.pR,l))
        sumW.init <- c(sumW.init, rep(SumW.init,l))
        sumW.pF <- c(sumW.pF, rep(SumW.pF,l))
        sumW.pR <- c(sumW.pR, rep(SumW.pR,l))
        n.pF <- c(n.pF, rep(N.pF,l))
        n.pR <- c(n.pR, rep(N.pR,l))
        w <- c(w, distr.new)
      } # close parDir
    } # close subDir
    
    df <- data.frame(isl.r, isl.c, E.mean, E.sd, prod.mean, prod.sd, incComp, incComp.time,
                     N.init=n.init, w.mean, w.CV,
                     m, move.fn=movefn,
                     p, prP, predOpt, predEnd, pred.fn=predfn,
                     st, feed.fn=feedfn,
                     f, babybump, repro.fn=reprofn,
                     MeanW.init=meanW.init, MeanW.pF=meanW.pF, VarW.pF=varW.pF, MeanW.pR=meanW.pR,
                     SumW.init=sumW.init, SumW.pF=sumW.pF, SumW.pR=sumW.pR,
                     N.pF=n.pF, N.pR=n.pR,
                     w=w)
    return(df)
  }









### Some of the simulations did not include 'incComp', 'incComp.time', or 'VarW.pF', so the above
### functions gave errors when reading in the data. This inserts NA's into the appropriate parameters.R files

# directory <- "~/Dropbox/Theoretical Ecology/Project/IslandRule/SimOutput"
# for(subDir in dir(directory, full.names=TRUE)) {
#   for(parDir in dir(subDir, full.names=TRUE)) {
#     source(paste0(parDir, "/parameters.R"))
#     
#     if(sum(names(island.pars)=="incComp")==0) {
#       temp.pars <- island.pars
#       island.pars <- vector("list", 8)
#       names(island.pars) <- c(names(temp.pars), "incComp", "incComp.time")
#       island.pars[1:6] <- temp.pars
#       island.pars[7:8] <- c(NA, NA)
#     }
#     if(sum(ls()=="VarW.pF")==0) {
#       VarW.pF <- NA
#     }
#     
#     # Generate filename for parameters
#     fname <- paste0(parDir, "/parameters.R")
#     
#     # Save parameters into the script
#     sink(fname)
#     
#     # Parameters
#     cat("island.pars <- list(isl.r=", island.pars$isl.r, ", ",
#         "isl.c=", island.pars$isl.c, ", ",
#         "E.mean=", island.pars$E.mean, ", ",
#         "E.sd=", island.pars$E.sd, ", ",
#         "prod.mean=", island.pars$prod.mean, ", ",
#         "prod.sd=", island.pars$prod.sd, ", ",
#         "incComp=", island.pars$incComp, ", ",
#         "incComp.time=", island.pars$incComp.time, ")",
#         "\n", sep="")
#     cat("pop.pars <- list(N.init=", pop.pars$N.init, ", ",
#         "w.mean=", pop.pars$w.mean, ", ",
#         "w.CV=", pop.pars$w.CV, ")",
#         "\n", sep="")
#     cat("move.pars <- list(m=", move.pars$m, ", ",
#         "move.fn='", move.pars$move.fn, "')",
#         "\n", sep="")
#     cat("pred.pars <- list(p=", pred.pars$p, ", ",
#         "prP=", pred.pars$prP, ", ",
#         "predOpt=", pred.pars$predOpt, ", ",
#         "predEnd=", pred.pars$predEnd, ", ",
#         "pred.fn='", pred.pars$pred.fn, "')",
#         "\n", sep="")
#     cat("feed.pars <- list(st=", feed.pars$st, ", ",
#         "feed.fn='", feed.pars$feed.fn, "')",
#         "\n", sep="")
#     cat("repro.pars <- list(f=", repro.pars$f, ", ",
#         "babybump=", repro.pars$babybump, ", ",
#         "repro.fn='", repro.pars$repro.fn, "')",
#         "\n", sep="")
#     
#     # Summary statistics
#     cat("nsims <- ", nsims, "\n", sep="")
#     cat("MeanW.init <- ", MeanW.init, "\n", sep="")
#     cat("MeanW.pF <- ", MeanW.pF, "\n", sep="")
#     cat("VarW.pF <- ", VarW.pF, "\n", sep="")
#     cat("MeanW.pR <- ", MeanW.pR, "\n", sep="")
#     cat("SumW.init <- ", SumW.init, "\n", sep="")
#     cat("SumW.pF <- ", SumW.pF, "\n", sep="")
#     cat("SumW.pR <- ", SumW.pR, "\n", sep="")
#     cat("N.init <- ", N.init, "\n", sep="")
#     cat("N.pF <- ", N.pF, "\n", sep="")
#     cat("N.pR <- ", N.pR, "\n", sep="")
#     sink()
#       
#     rm(list=c("temp.pars", "feed.pars", "island.pars", "MeanW.init", "MeanW.pF", "MeanW.pR", "move.pars", "N.init",
#               "N.pF", "N.pR", "nsims", "pop.pars", "pred.pars", "repro.pars", "SumW.init", "SumW.pF", "SumW.pR", "VarW.pF"))
#   }
# }
    
    
    
    
    
    
    
