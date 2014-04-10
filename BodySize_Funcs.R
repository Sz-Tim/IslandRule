### Theoretical Ecology
### 2014 Spring
### Tim Szewczyk
### Island body size changes



###############################################
###--- General functions used throughout ---###
###############################################

  logit <- function(x) {log(x/(1-x))}
  antilogit <- function(x) {exp(x)/(1+exp(x))}
  # sample(x) uses a vector of elements or, if given an integer, 1:x. To avoid that:
  resample <- function(x,...) {x[sample.int(length(x),...)]}
  
  # ALTER THIS FUNCTION
  writeDataAndParams <- function(theData, pars) {
    rundir <- paste0("Run", as.character(dirNum)) # Make a name for a directory using the directory number
    dir.create(rundir)  # Create a directory
    filename <- paste0(rundir, "/simoutput.txt")    # Create a name for a file
    write.matrix(theData, file=filename)    # Write the data to the file
  }


################################
###--- Movement Functions ---###
################################

  move <- function(N.df=N.df, m=m, rmax=isl.r, cmax=isl.c, fn="random") {
    
    if(fn=="random") {
      # Move
      N.df$loc.r <- N.df$loc.r + rbinom(nrow(N.df), 1, m)*sample(c(1,-1), nrow(N.df), replace=T)
      N.df$loc.c <- N.df$loc.c + rbinom(nrow(N.df), 1, m)*sample(c(1,-1), nrow(N.df), replace=T)
  
      # Enforce boundaries
      N.df$loc.r[which(N.df$loc.r < 1)] <- 1
      N.df$loc.r[which(N.df$loc.r > rmax)] <- rmax
      N.df$loc.c[which(N.df$loc.c < 1)] <- 1
      N.df$loc.c[which(N.df$loc.c > cmax)] <- cmax
    }
    return(N.df)
  }

  #moveToResource
  #moveFromOthers


#################################
###--- Predation Functions ---###
#################################

  ###--- Main predation function ---###
  predation <- function(N.df, p, prP, predOpt, fn="antilogit") {
    
    if(fn == "proportional") {
      predprob <- 1 - antilogit(N.df$w^p - mean(N.df$w^p))
      N.df$w <- N.df$w * rbinom(nrow(N.df), 1, prob=predprob)
      
    } else if(fn == "diff") {
      predprob <- (predOpt - N.df$w)^2/max((predOpt - N.df$w)^2)
      N.df$w <- N.df$w * rbinom(nrow(N.df), 1, prob=predprob)
      
    } else if(fn == "antilogit") {
      predprob <- antilogit(N.df$w^p - mean(N.df$w^p))
      N.df$w <- N.df$w * rbinom(nrow(N.df), 1, prob=predprob)
      
    } else if(fn == "even") {
      N.df$w <- N.df$w * rbinom(nrow(N.df), 1, prob=(1-prP))
      
    } else {print("Error: Unrecognized predation function")}
    
    N.df <- subset(N.df, N.df$w>0)
    return(N.df)
  }



###############################
###--- Feeding Functions ---###
###############################

###--- Starvation probability functions ---###
  starveProportional <- function(w.cell, st=st) {
    # Starvation probability scales with resource requirements
    # Larger individuals are more likely to starve
    prob <- (w.cell^st)/sum(w.cell^st)
    return(prob)
  }

  starveOutcompete <- function(w.cell, st=st) {
    # Starvation probability scales inversely with resource requirements
    # Smaller individuals are more likely to starve
    prob <- 1 - antilogit(w.cell^st - mean(w.cell^st))
    return(prob)
  }
  
  starveAntilogit <- function(w.cell, st=st) {
    # Starvation probability scales with resource requirements
    # Larger individuals are more likely to starve
    prob <- antilogit(w.cell^st - mean(w.cell^st))
    return(prob)
  }
  
  starveRandom <- function(w.cell) {
    prob <- rep(1/length(w.cell), length(w.cell))
    return(prob)
  }

###--- Main feeding function ---###
  feed <- function(N.df=N.df, E=E, st=st, fn="proportional") {
    for(i in 1:nrow(E)) {
      for(j in 1:ncol(E)) {
        inCell <- which(N.df$loc.r==i & N.df$loc.c==j)
        
        # Are there enough resources?
        if(fn == "proportional") {
          while( sum(N.df$w[inCell]) > E[i,j] ) {
            w.cell <- N.df$w[inCell]
            starveprob <- starveProportional(w.cell, st)
            starve <- resample(inCell, 1, prob=starveprob)
            N.df$w[starve] <- 0
          }
        } else if(fn == "outcompete") { 
          while( sum(N.df$w[inCell]) > E[i,j] ) {
            w.cell <- N.df$w[inCell]
            starveprob <- starveOutcompete(w.cell, st) 
            starve <- resample(inCell, size=1, prob=starveprob)
            N.df$w[starve] <- 0
          }
        } else if(fn == "antilogit") {
          while( sum(N.df$w[inCell]) > E[i,j] ) {
            w.cell <- N.df$w[inCell]
            starveprob <- starveAntilogit(w.cell, st) 
            starve <- resample(inCell, size=1, prob=starveprob)
            N.df$w[starve] <- 0
          }
        } else if(fn == "even") {
          while( sum(N.df$w[inCell]) > E[i,j] ) {
            w.cell <- N.df$w[inCell]
            starveprob <- starveRandom(w.cell) 
            starve <- resample(inCell, size=1, prob=starveprob)
            N.df$w[starve] <- 0
          }
        } else {print("Error: Unrecognized feeding function")}
          
        # Set post-feeding resource level
        E[i,j] <- E[i,j] - sum(N.df$w[inCell]) 
        N.df <- subset(N.df, N.df$w>0)
      }
    }
    return(list(E=E, N.df=N.df))
  }



####################################
###--- Reproduction Functions ---###
####################################

###--- fecundity = f(w) ---###
  #  proportional: lambda = N.df$w^f
  

###--- Main reproduction function ---###
  reproduce <- function(N.df=N.df, w.CV=w.CV, f=f, babybump=0, fn="proportional") {
    
    if(fn == "proportional") {
      lam <- N.df$w*f
      
    } else if(fn == "log") {
      lam <- log(N.df$w + 1)*f + babybump
      
    } else if(fn == "even") {
      lam <- runif(nrow(N.df), min=0, max=10)
      lam[N.df$w==0] <- 0
    }
    
    offspring <- rpois(nrow(N.df), lam)
    new.df <- data.frame(loc.r=rep(N.df$loc.r, times=offspring),
                         loc.c=rep(N.df$loc.c, times=offspring),
                         w=rtruncnorm(sum(offspring), a=0, b=Inf, mean=rep(N.df$w, times=offspring), sd=w.CV*rep(N.df$w, times=offspring)))
    return(N.df=new.df)
  }
  

####################################
###--- Productivity Functions ---###
####################################

  grow <- function(E=E, MN=prod.mean, SD=prod.sd) {
    E <- E + rnorm(prod(dim(E)), MN, SD)
    return(E)
  }
  
  


