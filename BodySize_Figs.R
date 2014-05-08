### Theoretical Ecology
### 2014 Spring
### Tim Szewczyk
### Island body size changes


####################
###--- Set Up ---###
####################
setwd("~/Dropbox/Theoretical Ecology/Project/IslandRule")
baseDir <- getwd()
library(truncnorm); library(plyr); library(ggplot2); library(MASS); library(gridExtra)
theme_set(theme_bw())
source(paste0(baseDir, "/BodySize_Funcs.R"))
source(paste0(baseDir, "/BodySize_GenFuncs.R"))
source(paste0(baseDir, "/BodySize_SimFunc.R"))

############################  
###--- Set parameters ---###
############################
sims <- 50
maxt <- 1000
parLen <- 8
mainland <- FALSE

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
                  feed.fn="antilogit")   # Starvation probability function: antilogit, proportional, outcompete, even

repro.pars <- list(f=1,                 # How reproductive rate scales with w; how it works depends on repro.fn
                   babybump=1,          # Increase lambda for all individuals; necessary for w < ~5 to avoid fatal errors
                   repro.fn="log")      # Reproduction function to calculate individual lambda: proportional, log, even


parSet <- makeParSet(param="incComp", low=1, high=18, len=parLen, logSeq=TRUE, sims=sims, maxt=maxt)
dirNum <- 1
summary.df <- summarizeAll(paste0(getwd(),"/SimOutput"))
byTime.df <- getSimsByTime(feedFunc=feed.pars$feed.fn, param=parSet$param, indices=1:parLen)#, mainland=mainland)
init.distr <- getDistr(paste0(getwd(),"/SimOutput"), time="init", sampsize=1000)
final.distr <- getDistr(paste0(getwd(),"/SimOutput"), time="final", sampsize=2000) 

ggplot(summary.df, aes(x=incComp, y=MeanW.pF, colour=feed.fn)) + geom_point(size=3, alpha=0.9) + stat_smooth(method="loess", se=FALSE)
ggplot(summary.df, aes(x=SumW.pF, colour=feed.fn)) + geom_density()

time.plot <- ggplot(byTime.df, aes(x=time, group=incComp, colour=incComp))
time.plot + geom_line(aes(y=MeanW.pF)) + ylim(0,40)
time.plot + geom_line(aes(y=VarW.pF))
time.plot + geom_line(aes(y=SumW.pF))
time.plot + geom_line(aes(y=N.pF))
# ggplot(summary.df, aes(x=MeanW.init, y=(MeanW.pF-MeanW.init))) + geom_line() + geom_hline(y=0, linetype=2)

ggplot(final.distr, aes(x=feed.fn, y=w)) + geom_boxplot(fill=NA)


######################################
### plots of probability functions ###
######################################

  N.df <- data.frame(w=seq(0.01, 18, by=0.01))
  p <- pred.pars$p
  prP <- pred.pars$prP
  predOpt <- pred.pars$predOpt
  st <- feed.pars$st

  # Predation
    pred.proportional <- antilogit(N.df$w^p - mean(N.df$w^p))
    pred.antilogit <- 1 - antilogit(N.df$w^p - mean(N.df$w^p))
    pred.diff <- 1 - (predOpt - N.df$w)^2/max((predOpt - N.df$w)^2)
    pred.even <- rep(prP, length(N.df$w))
    plot(N.df$w, pred.proportional, type="l", ylim=c(0,1), xlab="w", ylab="Pr(predation)")
    plot(N.df$w, pred.antilogit, type="l", ylim=c(0,1), xlab="w", ylab="Pr(predation)")
    plot(N.df$w, pred.diff, type="l", ylim=c(0,1), xlab="w", ylab="Pr(predation)")
    plot(N.df$w, pred.even, type="l", ylim=c(0,1), xlab="w", ylab="Pr(predation)")

  # Feeding
    feed.proportional <- starveProportional(N.df$w, st)
    feed.antilogit <- starveAntilogit(N.df$w, st)
    feed.outcompete <- starveOutcompete(N.df$w, st)
    feed.even <- starveRandom(N.df$w)
    plot(N.df$w, feed.proportional, type="l", xlab="w", ylab="Pr(Starvation)")
    plot(N.df$w, feed.antilogit, type="l", xlab="w", ylab="Pr(Starvation)")
    plot(N.df$w, feed.outcompete, type="l", xlab="w", ylab="Pr(Starvation)")
    plot(N.df$w, feed.even, type="l", xlab="w", ylab="Pr(Starvation)")



######################################
### plots of w.mean with st varied ###
######################################
st.df <- subset(summary.df, is.na(summary.df$incComp) & summary.df$w.mean==1 & summary.df$m==0.5 & summary.df$isl.r==5)
ggplot(st.df, aes(x=st, y=MeanW.pF, colour=feed.fn)) + geom_point(size=3) + labs(y="Mean W")

alog.df <- getSimsByTime(feedFunc="antilogit", param="st", indices=1:8)
prop.df <- getSimsByTime(feedFunc="proportional", param="st", indices=1:8)
even.df <- getSimsByTime(feedFunc="even", param="st", indices=1)
outcomp.df <- getSimsByTime(feedFunc="outcompete", param="st", indices=1:7)

alog.g <- ggplot(alog.df, aes(x=time, group=st, colour=st, y=MeanW.pF)) + geom_line() + ylim(0,30) + labs(y="Mean W") 
prop.g <- ggplot(prop.df, aes(x=time, group=st, colour=st, y=MeanW.pF)) + geom_line() + ylim(0,30) + labs(y="Mean W")
even.g <- ggplot(even.df, aes(x=time, y=MeanW.pF, colour=st)) + geom_line() + ylim(0,30) + labs(y="Mean W")
outcomp.g <- ggplot(outcomp.df, aes(x=time, group=st, colour=st, y=MeanW.pF)) + geom_line() + ylim(0,30) + labs(y="Mean W")

grid.arrange(alog.g, prop.g, even.g, outcomp.g, ncol=2, nrow=2)


##########################################
### plots of w.mean with w.init varied ###
##########################################
w.mean.df <- subset(summary.df, is.na(summary.df$incComp) & summary.df$st==0.3 & summary.df$m==0.5 & summary.df$isl.r==5)
ggplot(w.mean.df, aes(x=w.mean, y=MeanW.pF, colour=feed.fn)) + geom_point(size=3) + labs(y="Mean W", x="w.init")

alog.df <- getSimsByTime(feedFunc="antilogit", param="w.mean", indices=1:8)
prop.df <- getSimsByTime(feedFunc="proportional", param="w.mean", indices=1:8)
even.df <- getSimsByTime(feedFunc="even", param="w.mean", indices=1:8)
outcomp.df <- getSimsByTime(feedFunc="outcompete", param="w.mean", indices=1:8)

alog.g <- ggplot(alog.df, aes(x=time, group=w.mean, colour=w.mean, y=MeanW.pF)) + geom_line() + ylim(0,20) + 
  labs(y="Mean W") + scale_colour_gradient(name="W.init")
prop.g <- ggplot(prop.df, aes(x=time, group=w.mean, colour=w.mean, y=MeanW.pF)) + geom_line() + ylim(0,20) + 
  labs(y="Mean W") + scale_colour_gradient(name="W.init")
even.g <- ggplot(even.df, aes(x=time, y=MeanW.pF, group=w.mean, colour=w.mean)) + geom_line() + ylim(0,20) + 
  labs(y="Mean W") + scale_colour_gradient(name="W.init")
outcomp.g <- ggplot(outcomp.df, aes(x=time, group=w.mean, colour=w.mean, y=MeanW.pF)) + geom_line() + ylim(0,20) + 
  labs(y="Mean W") + scale_colour_gradient(name="W.init")

grid.arrange(alog.g, prop.g, even.g, outcomp.g, ncol=2, nrow=2)



###########################################
### plots of w.mean with incComp varied ###
###########################################
incComp.df <- subset(summary.df, is.na(summary.df$incComp)==FALSE)
ggplot(incComp.df, aes(x=incComp, y=MeanW.pF, colour=feed.fn)) + geom_vline(x=20, linetype=2, colour="gray70") + 
  geom_point(size=3) + labs(y="Mean W", x="Productivity")
ggplot(incComp.df, aes(x=incComp, y=VarW.pF, colour=feed.fn)) + geom_vline(x=20, linetype=2, colour="gray70") +
  geom_point(size=3) + labs(y="Variance in w", x="Productivity")
ggplot(incComp.df, aes(x=MeanW.pF, y=sqrt(VarW.pF), colour=feed.fn)) + geom_point(aes(size=incComp), alpha=0.75) + 
  labs(y="Standard deviation in W", x="Mean W")

alog.df <- getSimsByTime(feedFunc="antilogit", param="incComp", indices=1:8)
prop.df <- getSimsByTime(feedFunc="proportional", param="incComp", indices=1:8)
even.df <- getSimsByTime(feedFunc="even", param="incComp", indices=1:8)
outcomp.df <- getSimsByTime(feedFunc="outcompete", param="incComp", indices=1:8)

alog.g <- ggplot(alog.df, aes(x=time, group=incComp, colour=incComp, y=MeanW.pF)) + geom_line() + ylim(0,40) + 
  labs(y="Mean W") + scale_colour_gradient(name="Productivity")
prop.g <- ggplot(prop.df, aes(x=time, group=incComp, colour=incComp, y=MeanW.pF)) + geom_line() + ylim(0,40) + 
  labs(y="Mean W") + scale_colour_gradient(name="Productivity")
even.g <- ggplot(even.df, aes(x=time, y=MeanW.pF, group=incComp, colour=incComp)) + geom_line() + ylim(0,40) + 
  labs(y="Mean W") + scale_colour_gradient(name="Productivity")
outcomp.g <- ggplot(outcomp.df, aes(x=time, group=incComp, colour=incComp, y=MeanW.pF)) + geom_line() + ylim(0,40) + 
  labs(y="Mean W") + scale_colour_gradient(name="Productivity")

grid.arrange(alog.g, prop.g, even.g, outcomp.g, ncol=2, nrow=2)

ggplot(alog.df, aes(x=time, group=incComp, colour=incComp, y=sqrt(VarW.pF))) + geom_line() + ylim(0,7) +
  labs(y="Standard deviation in W") + scale_colour_gradient(name="Productivity")
ggplot(prop.df, aes(x=time, group=incComp, colour=incComp, y=sqrt(VarW.pF))) + geom_line() + ylim(0,7) + 
  labs(y="Standard deviation in W") + scale_colour_gradient(name="Productivity")
ggplot(even.df, aes(x=time, y=sqrt(VarW.pF), group=incComp, colour=incComp)) + geom_line() + ylim(0,7) + 
  labs(y="Standard deviation in W") + scale_colour_gradient(name="Productivity")
ggplot(outcomp.df, aes(x=time, group=incComp, colour=incComp, y=sqrt(VarW.pF))) + geom_line() + ylim(0,7) + 
  labs(y="Standard deviation in W") + scale_colour_gradient(name="Productivity")



#####################################
### plots of w.mean with f varied ###
#####################################
f.df <- subset(summary.df, summary.df$f!=1 & summary.df$st==0.3 & summary.df$m==0.5 & summary.df$isl.r==5 & summary.df$w.mean==0.5)
ggplot(f.df, aes(x=f, y=MeanW.pF, colour=feed.fn)) + geom_point(size=3) + labs(y="Mean W", x="f")

alog.df <- getSimsByTime(feedFunc="antilogit", param="f", indices=1:8)
prop.df <- getSimsByTime(feedFunc="proportional", param="f", indices=1:8)
even.df <- getSimsByTime(feedFunc="even", param="f", indices=1:8)
outcomp.df <- getSimsByTime(feedFunc="outcompete", param="f", indices=1:8)

alog.g <- ggplot(alog.df, aes(x=time, group=f, colour=f, y=MeanW.pF)) + geom_line() + ylim(0,21) + 
  labs(y="Mean W") + scale_colour_gradient(name="f")
prop.g <- ggplot(prop.df, aes(x=time, group=f, colour=f, y=MeanW.pF)) + geom_line() + ylim(0,21) + 
  labs(y="Mean W") + scale_colour_gradient(name="f")
even.g <- ggplot(even.df, aes(x=time, y=MeanW.pF, group=f, colour=f)) + geom_line() + ylim(0,21) + 
  labs(y="Mean W") + scale_colour_gradient(name="f")
outcomp.g <- ggplot(outcomp.df, aes(x=time, group=f, colour=f, y=MeanW.pF)) + geom_line() + ylim(0,21) + 
  labs(y="Mean W") + scale_colour_gradient(name="f")

grid.arrange(alog.g, prop.g, even.g, outcomp.g, ncol=2, nrow=2)




################################
### plots of w distributions ###
################################
# w.mean
  wInit.distr <- subset(init.distr, is.na(init.distr$incComp) & init.distr$st==0.3 & init.distr$m==0.5 & init.distr$isl.r==5)
  wFinal.distr <- subset(final.distr, is.na(final.distr$incComp) & final.distr$st==0.3 & final.distr$m==0.5 & final.distr$isl.r==5)

  ggplot(wFinal.distr, aes(x=w, colour=factor(w.mean))) + geom_density(trim=TRUE) + facet_wrap(~feed.fn) + 
    scale_colour_brewer(name="Mean w.init", type="seq", palette=4) + labs(x="Final w")
  ggplot(wInit.distr, aes(x=w, colour=factor(w.mean))) + geom_density(trim=TRUE) + facet_wrap(~feed.fn) + 
    scale_colour_brewer(name="Mean w.init", type="seq", palette=4) + labs(x="Initial w")


# st
  stInit.distr <- subset(init.distr, is.na(init.distr$incComp) & init.distr$w.mean==1 & init.distr$m==0.5 & init.distr$isl.r==5)
  stFinal.distr <- subset(final.distr, is.na(final.distr$incComp) & final.distr$w.mean==1 & final.distr$m==0.5 & final.distr$isl.r==5)

  ggplot(stFinal.distr, aes(x=w, colour=factor(st))) + geom_density(trim=TRUE) + facet_wrap(~feed.fn) + 
    scale_colour_brewer(name="st", type="seq", palette=4) + labs(x="Final w")
  ggplot(stInit.distr, aes(x=w, colour=factor(st))) + geom_density(trim=TRUE) + facet_wrap(~feed.fn) + 
    scale_colour_brewer(name="st", type="seq", palette=4) + labs(x="Initial w")


# incComp
  incCompInit.distr <- subset(init.distr, is.na(init.distr$incComp)==FALSE)
  incCompFinal.distr <- subset(final.distr, is.na(final.distr$incComp)==FALSE)
  incCompInit.distr$compFac <- factor(incCompInit.distr$incComp)
  levels(incCompInit.distr$compFac)[3:4] <- "15"
  levels(incCompInit.distr$compFac)[4:5] <- "20"
  levels(incCompInit.distr$compFac)[5:6] <- "25"
  levels(incCompInit.distr$compFac)[6:7] <- "30"
  levels(incCompInit.distr$compFac)[7:8] <- "35"
  levels(incCompInit.distr$compFac)[8:9] <- "40"
  incCompFinal.distr$compFac <- factor(incCompFinal.distr$incComp)
  levels(incCompFinal.distr$compFac)[3:4] <- "15"
  levels(incCompFinal.distr$compFac)[4:5] <- "20"
  levels(incCompFinal.distr$compFac)[5:6] <- "25"
  levels(incCompFinal.distr$compFac)[6:7] <- "30"
  levels(incCompFinal.distr$compFac)[7:8] <- "35"
  levels(incCompFinal.distr$compFac)[8:9] <- "40"
  

  ggplot(incCompFinal.distr, aes(x=w, colour=compFac)) + geom_density(trim=TRUE) + facet_wrap(~feed.fn) + 
    scale_colour_brewer(name="Productivity", type="seq", palette=4) + labs(x="Final w")
  ggplot(incCompInit.distr, aes(x=w, colour=compFac)) + geom_density(trim=TRUE) + facet_wrap(~feed.fn) + 
    scale_colour_brewer(name="Productivity", type="seq", palette=4) + labs(x="Initial w")












