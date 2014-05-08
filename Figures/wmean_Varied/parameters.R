island.pars <- list(isl.r=5, isl.c=5, E.mean=20, E.sd=0, prod.mean=20, prod.sd=0)
pop.pars <- list(N.init=40, w.mean=0.3, w.CV=0.04)  # w.mean varied within a graph
move.pars <- list(m=0.5, move.fn='random')
pred.pars <- list(p=0.5, prP=0.2, predOpt=10, predEnd=300, pred.fn='none')
feed.pars <- list(st=0.3, feed.fn='antilogit')   # feed.fn varied across graphs
repro.pars <- list(f=1, babybump=1, repro.fn='log')
nsims <- 50
MeanW.init <- 0.2987469
MeanW.pF <- 3.300375
MeanW.pR <- 3.360419
SumW.init <- 11.94988
SumW.pF <- 500.7391
SumW.pR <- 1248.523
N.init <- 40
N.pF <- 153.8
N.pR <- 374.98
