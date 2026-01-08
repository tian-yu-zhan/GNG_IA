#Flat prior with mean 0 and large SD:
Pr.diff <- 0
Pr.sig <- 1e6
#Final N per arm (equal randomization)
N <-50
#TPPs and tuning parameters:
minTPP <- 1.5
baseTPP <- 3
P.min <- 0.8
P.base <- 0.2
#P.nogo <- 0.65
#Interim with M per arm (M < N)
M <- 30
#Interim accelerate if pred. probability greater than a threshold:
P.int <- c(1,0.6,0.8)
#Number of simulations in Monte Carlo:
N.sim <- 500000

#Known sample SD
#Assuming if we have a known true DIFF
Calc_Interim <- function(sig=1,DIFF=0,M=30)
{
  ###STUDY END###
  #Posterior SD
  sig.post <- sqrt(1/(1/(Pr.sig^2)+1/(2*(sig^2)/N)))
  #Minimal posterior mean needed for GO
  t.Go <- max(minTPP+sig.post*qnorm(P.min), baseTPP+sig.post*qnorm(P.base))
  #Minimal observed diff at study end for GO
  Go <- (t.Go*(1/(Pr.sig^2)+1/(2*(sig^2)/N))-Pr.diff/(Pr.sig^2))*(2*(sig^2)/N)
  #Maximal posterior mean needed for NO-GO
  #t.No <- min(minTPP+sig.post*qnorm(P.nogo), baseTPP+sig.post*qnorm(P.base))
  #Maximal observed diff at study end for NO-GO
  #No <- (t.No*(1/(Pr.sig^2)+1/(2*(sig^2)/N))-Pr.diff/(Pr.sig^2))*(2*(sig^2)/N)
  
  #### INTERIM ####
  #Posterior at interim with mean at avg(the first M diff) and SD
  sig.int <- sqrt(1/(1/(Pr.sig^2)+1/(2*(sig^2)/M)))
  #Predictive prob for the remainder (N-M) diff with SD
  sig.pred <- sqrt(2*(sig^2)+sig.int^2)
  #Minimal average of the first M diff needed for interim accelerate
  Go.int  <- Go+qnorm(P.int)*sqrt(N-M)*sig.pred/N
  
  #Probability of making a "Go" decision at final
  Prob.final <- 1-pnorm(Go,mean=DIFF,sd=2*(sig^2)/N)
  
  #Probability of interim accelerate
  Prob.int <- 1-pnorm(Go.int,mean=DIFF,sd=2*(sig^2)/M)
  
  #Probability of mismatch -- interim greater than Go.int but final less than Go
  ANYGO <- NULL
  INTGO <- NULL
  for (i in 1:length(P.int))
  {
    temp.int <- rnorm(N.sim, mean=DIFF, sd=2*(sig^2)/M)
    temp.rest <- rnorm(N.sim, mean=DIFF, sd=2*(sig^2)/(N-M))
    ANYGO[i] <- length(temp.int[temp.int>Go.int[i] | (temp.int*M+temp.rest*(N-M))/N>Go])/N.sim
    INTGO[i] <- length(temp.int[temp.int>Go.int[i]])/N.sim
  }
  
  out.data <- data.frame(Go,Prob.final,Go.int,Prob.int,ANYGO,INTGO)
  return(out.data)
}


Calc_Interim(sig=2.5,DIFF=1.5,M=30)
Calc_Interim(sig=2.5,DIFF=2.0,M=30)
Calc_Interim(sig=2.5,DIFF=2.5,M=30)
Calc_Interim(sig=2.5,DIFF=3.0,M=30)


###PLOT THE OUTPUT###
library(latex2exp)
DIFF_plot <- seq(1.5,3.0,0.02)
ANYGO_plot=matrix(NA,nrow=3,ncol=length(DIFF_plot))
INTGO_plot=matrix(NA,nrow=3,ncol=length(DIFF_plot))
for (i in 1:length(DIFF_plot))
{
  temp <- Calc_Interim(sig=2.5,DIFF=DIFF_plot[i],M=30)
  ANYGO_plot[,i] <- temp$ANYGO
  INTGO_plot[,i] <- temp$INTGO
}

  
plot(DIFF_plot,INTGO_plot[2,],type='l',col='green',ylim=c(0,1),xlim=c(1.5,3.0),
     xlab='True Treatment Effect',ylab='Probability of Go',axes=FALSE,
     main=TeX(r'(Sample SD, $\sigma$ = 2.5)'))
points(DIFF_plot,INTGO_plot[3,],type='l',col='green',lty=2)
points(DIFF_plot,INTGO_plot[1,],type='l',col='green',lty=3)
points(DIFF_plot,ANYGO_plot[2,],type='l',lty=1)
points(DIFF_plot,ANYGO_plot[3,],type='l',lty=2)
points(DIFF_plot,ANYGO_plot[1,],type='l',lty=3)
axis(1,at=c(1.5,1.8,2.1,2.4,2.7,3.0),tck=0.02)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),tck=0.02)
mtext("Interim with N=30 out of 50",side=4)
lines(x=c(1.5,3.0),y=c(0.2,0.2),lty=6,col=gray(0.8,0.5))
lines(x=c(1.5,3.0),y=c(0.4,0.4),lty=6,col=gray(0.8,0.5))
lines(x=c(1.5,3.0),y=c(0.6,0.6),lty=6,col=gray(0.8,0.5))
lines(x=c(1.5,3.0),y=c(0.8,0.8),lty=6,col=gray(0.8,0.5))
lines(x=c(1.5,3.0),y=c(1.0,1.0),lty=6,col=gray(0.8,0.5))

legend(1.6,0.75,lty=c(1,2,3),col=c(1,1,1),
       legend=c('Interim threshold = 0.6','Interim threshold = 0.8','No interim threshold'),box.lty=0)
legend(1.6,0.95,lty=c(1,1),col=c(1,3),
       legend=c('Probability of Go at any time','Probability of Go at interim'),box.lty=0)


