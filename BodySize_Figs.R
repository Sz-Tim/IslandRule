library(ggplot2)

mn.mn.anti <- rowMeans(mns.pR.anti)
mn.mn.even <- rowMeans(mns.pR.even)
mn.mn.outC <- rowMeans(mns.pR.outC)
mn.mn.prop <- rowMeans(mns.pR.prop)

matplot(1:maxt, mns.pR.anti, lty=1, type="l", col=rgb(0,0,0,0.5), ylim=c(0,12), xlab="Generations", ylab="Mean w after reproduction")
lines(1:maxt, mn.mn.anti, lwd=3, col="black")
matlines(1:maxt, mns.pR.even, lty=1, type="l", col=rgb(1,0,0,0.5))
lines(1:maxt, mn.mn.even, lwd=3, col="darkred")
matlines(1:maxt, mns.pR.outC, lty=1, type="l", col=rgb(0,1,0,0.5))
lines(1:maxt, mn.mn.outC, lwd=3, col="darkgreen")
matlines(1:maxt, mns.pR.prop, lty=1, type="l", col=rgb(0,0,1,0.5))
lines(1:maxt, mn.mn.prop, lwd=3, col="darkblue")
abline(v=predEnd, lty=3)

legend("topleft", bty="n", c("Outcompete", "Even", "Antilogit", "Proportion"), col=c("green", "red", "black", "blue"), lty=1, lwd=3)
