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
M <- seq(10,N-10,10)
#Interim accelerate if pred. probability greater than a threshold:
P.int <- seq(0.02,0.98,0.02)
#Number of simulations in Monte Carlo:
N.sim <- 100000

#Known sample SD
#Assuming if we have a known true DIFF
Calc_Interim <- function(sig=1,DIFF=2.7,M=10)
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
  MISMATCH <- NULL
  for (i in 1:length(P.int))
  {
    temp.int <- rnorm(N.sim, mean=DIFF, sd=2*(sig^2)/M)
    temp.rest <- rnorm(N.sim, mean=DIFF, sd=2*(sig^2)/(N-M))
    MISMATCH[i] <- length(temp.int[temp.int>Go.int[i] & (temp.int*M+temp.rest*(N-M))/N<Go])/N.sim
  }
  
  out.data <- data.frame(Go,Prob.final,Go.int,Prob.int,MISMATCH)
  return(out.data)
}

###PLOT THE OUTPUT###
library(latex2exp)

temp <- Calc_Interim(sig=2.5,DIFF=2.8,M=10)
plot(P.int,temp$Prob.int,type='l',col='green',ylim=c(0,1),xlim=c(0,1),
     xlab='',ylab='Probility of Go',axes=FALSE,
     main=TeX(r'(Sample SD, $\sigma$ = 2.5)'))
axis(1,at=c(0,0.2,0.4,0.6,0.8,1),tck=0.02)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),tck=0.02)
mtext("Interim with N=10 out of 50",side=4)
lines(x=c(0,1),y=c(temp$Prob.final[1],temp$Prob.final[1]),lty=2,col=4)
points(P.int,temp$MISMATCH,type='l')

temp <- Calc_Interim(sig=2.5,DIFF=2.8,M=20)
plot(P.int,temp$Prob.int,type='l',col='green',ylim=c(0,1),xlim=c(0,1),
     xlab='',ylab='',axes=FALSE,
     main='True Effect = 2.8')
axis(1,at=c(0,0.2,0.4,0.6,0.8,1),tck=0.02)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),tck=0.02)
axis(4,tick=FALSE,labels=FALSE)
mtext("Interim with N=20 out of 50",side=4)
lines(x=c(0,1),y=c(temp$Prob.final[1],temp$Prob.final[1]),lty=2,col=4)
points(P.int,temp$MISMATCH,type='l')

temp <- Calc_Interim(sig=2.5,DIFF=2.8,M=30)
plot(P.int,temp$Prob.int,type='l',col='green',ylim=c(0,1),xlim=c(0,1),
     xlab='Pred. Probability Threshold at Interim',ylab='Probility of Go',axes=FALSE,
     main='')
axis(1,at=c(0,0.2,0.4,0.6,0.8,1),tck=0.02)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),tck=0.02)
mtext("Interim with N=30 out of 50",side=4)
lines(x=c(0,1),y=c(temp$Prob.final[1],temp$Prob.final[1]),lty=2,col=4)
points(P.int,temp$MISMATCH,type='l')
legend(0.05,0.6,lty=c(2,1,1),col=c(4,3,1),
       legend=c('Final Go',
                'Interim Go',
                'Discordant Go'),y.intersp=2)

temp <- Calc_Interim(sig=2.5,DIFF=2.8,M=40)
plot(P.int,temp$Prob.int,type='l',col='green',ylim=c(0,1),xlim=c(0,1),
     xlab='Pred. Probability Threshold at Interim',ylab='Probility of Go',axes=FALSE,
     main='')
axis(1,at=c(0,0.2,0.4,0.6,0.8,1),tck=0.02)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),tck=0.02)
axis(4,tick=FALSE,labels=FALSE)
mtext("Interim with N=40 out of 50",side=4)
lines(x=c(0,1),y=c(temp$Prob.final[1],temp$Prob.final[1]),lty=2,col=4)
points(P.int,temp$MISMATCH,type='l')

