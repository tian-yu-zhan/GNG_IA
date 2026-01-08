library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(234984675)

#####################################
# function to calculate beta binomial predictive probability

pbetabinom_c <- function (q, size, m, s) {
  # if (any(q < 0))
  #     stop("q must contain non-negative values")
  if (any(size < 0))
    stop("size must contain non-negative values")
  if (any(m <= 0) || any(m >= 1))
    stop("m must lie between 0 and 1")
  if (any(s <= 0))
    stop("s must be positive")
  len <- max(length(q), length(m), length(s), length(size))
  if (length(q) != len) {
    if (length(q) == 1)
      q <- rep(q, len)
    else stop("length of q incorrect")
  }
  if (length(size) != len) {
    if (length(size) == 1)
      size <- rep(size, len)
    else stop("size must be the same length as q")
  }
  if (length(m) != len) {
    if (length(m) == 1)
      m <- rep(m, len)
    else stop("m and q must have the same length")
  }
  if (length(s) != len) {
    if (length(s) == 1)
      s <- rep(s, len)
    else stop("s and q must have the same length")
  }
  if (any(q > size)){
    # Updating to correctly deal with this siduation
    # stop("q must be <= size")
    if(all(q>size)) return(rep(1,len))
    Val <- q<= size
    out <- rep(1,len)
    out[Val] <- pbetabinom_c(q[Val],size[Val],m[Val],s[Val])
    return(out)
  }
  if (any(q < 0)){
    # stop("q must contain non-negative values")
    if(all(q < 0)) return(rep(0,len))
    Val <- q>=0
    out <- rep(0,len)
    out[Val] <- pbetabinom_c(q[Val],size[Val],m[Val],s[Val])
    return(out)
  }
  t <- s * m
  u <- s * (1 - m)
  res <- vector("numeric", length(q))
  for (i in 1:length(q)) {
    qq <- 0:q[i]
    res[i] <- sum(exp(lbeta(qq + t[i], size[i] - qq + u[i]) -
                        lbeta(t[i], u[i]) + lchoose(size[i], qq)))
  }
  res
}

#####################################################
# Setup Study design and calculate decision thresholds for final and interim analysis.

Params <- tibble(
  alpha = .5, # shape 1 parameter for the prior beta distribution.
  beta = .5, # shape 2 parameter for the prior beta distribution.
  Final.N = 40, # final sample size.
  Interm.Start = 1, # Point at which interim analysis will start.
  Delta.lrv = c(.1, .25, .5, .75, .9), #Lower reference value for final decision making
  Delta.tv = c(.2, .35, .6, .85, .95), # target value for final decision making
  tau.lrv = .8, # probability threshold associated with the lower reference value for final decision making.
  tau.tv = .2, # probability threshold associated with the target value for final decision making
  tau.ng = .65, # probability threshold associated with the no go decision for final decision
  tau.interim.go = .8, # Probability threshold for the interim go decision
  tau.interim.nogo = .8 # probability threshold for the interim no-go decision (Not advised or reported in these outputs)
) %>%
  mutate(Senario = 1:n())

# Simulation Parameters. Difference is used to define the true rate in the simulations and is the difference from the lower reference value (Delta.lrv)
Difference <- c(-.05,0,0.05,.1)
nSim <- 5000

# probability for all possible outcomes
Params.Detailed <- crossing(Senario = Params$Senario, Final.Resp = 1:max(Params$Final.N)) %>%
  left_join(Params,by = 'Senario') %>%
  filter(Final.Resp <= Final.N) %>%
  mutate(
    Prob.tv = 1 - pbeta(Delta.tv,alpha+Final.Resp,beta + (Final.N - Final.Resp)),
    Prob.lrv = 1 - pbeta(Delta.lrv,alpha+Final.Resp,beta + (Final.N - Final.Resp)),
    Decision = case_when(
      Prob.lrv > tau.lrv & Prob.tv > tau.tv ~ 'Go',
      Prob.lrv <= tau.ng & Prob.tv <= tau.tv ~ 'No Go',
      T ~ 'Consider'
    )
  ) %>%
  group_by(Senario) %>%
  mutate(
    Min.Go = min(Final.Resp[Decision == 'Go']),
    Max.NoGo = max(Final.Resp[Decision == 'No Go'])
    )

# Collapse to just the decision threshold points.
Thresholds <- Params.Detailed %>%
  select(-Prob.tv, - Prob.lrv, - Decision, - Final.Resp) %>%
  group_by(Senario) %>%
  slice(1)

# Interim predictive go probabilities for all outcomes.
Interim.Detailed <- crossing(
  Senario = Params$Senario,
  Interm.Resp = 0:max(Params$Final.N),
  Interm.N = min(Params$Interm.Start):max(Params$Final.N)
) %>%
  left_join(Thresholds,by = 'Senario') %>%
  filter(
    Interm.Resp <= Interm.N,
    Interm.N <= Final.N,
    Interm.N >= Interm.Start
    ) %>%
  rowwise() %>%
  mutate(
    Prob.Go = 1-pbetabinom_c(Min.Go-Interm.Resp-1,Final.N-Interm.N,(alpha+Interm.Resp)/(alpha+beta+Interm.N),alpha+beta+Interm.N),
    Prob.NoGo = pbetabinom_c(Max.NoGo-Interm.Resp,Final.N-Interm.N,(alpha+Interm.Resp)/(alpha+beta+Interm.N),alpha+beta+Interm.N)
  ) %>%
  ungroup() %>%
  mutate(
    Int.Go = Prob.Go >= tau.interim.go,
    Int.NoGo = Prob.NoGo >= tau.interim.nogo
  )

# Interim decision thresholds
Interim <- Interim.Detailed %>%
  group_by(Senario,Interm.N) %>%
  summarise(
  int.Min.Go = min(Interm.Resp[Int.Go]),
  int.Max.NoGo = max(Interm.Resp[Int.NoGo])
  )

#########################################################
# Simulations


# Runs nSim simulations each of Final.N subjects.
Sim.Raw <- crossing(Senario = 1:nrow(Params),Difference = Difference,Sim = 1:nSim,Subject = 1:max(Params$Final.N)) %>%
  left_join(Thresholds,by = 'Senario') %>%
  filter(Subject <= Final.N) %>%
  mutate(Response = rbinom(n(),1,Delta.lrv+Difference)) %>%
  left_join(Interim,by = c('Senario',Subject = 'Interm.N')) %>%
  group_by(Senario,Sim,Delta.lrv,Difference) %>%
  mutate(CumulativeResponses = cumsum(Response)) %>%
  ungroup() %>%
  mutate(
    Interim.Decision = case_when(
      CumulativeResponses >= int.Min.Go ~ "Go",
      CumulativeResponses <= int.Max.NoGo ~ "No Go",
      T ~ "Contenue"
    )
  )

# Calculates the final decision for each simulation.
FinalDecison <- Sim.Raw %>%
  filter(Final.N == Subject) %>%
  mutate(
    Final.Decision = case_when(
      CumulativeResponses >= Min.Go ~ 'Go',
      CumulativeResponses <= Max.NoGo ~ 'No Go',
      T ~ "Consider"
    )) %>%
  select(Senario,Sim,Delta.lrv,Difference,Final.Decision)

# Aggregates the simulations into interim go probability for each interim point and the associated discordant go probability.
Sim.Results <- Sim.Raw %>%
  filter(!is.na(int.Min.Go)) %>%
  left_join(FinalDecison, by = c('Senario','Sim','Delta.lrv','Difference')) %>%
  group_by(Senario,Subject,Delta.lrv,Difference) %>%
  summarise(
    Go = sum(Interim.Decision == 'Go')/nSim,
    `Discordant Go` = sum(Interim.Decision == 'Go' & Final.Decision != 'Go')/nSim,
    .groups = 'drop'
  )

# Generates table shown in manuscript.
Sim.Results %>%
  mutate(
    DifferenceLab = paste('LRV ',ifelse(Difference<0,'-','+'),' ',abs(Difference)),
    DifferenceLab = ifelse(Difference == Difference[1],paste('True Rate =',DifferenceLab),DifferenceLab),
    DifferenceLab = factor(DifferenceLab,levels = unique(DifferenceLab)[order(unique(Difference))]),
    lrvLab = paste('LRV =',Delta.lrv),
    lrvLab = factor(lrvLab,levels = unique(lrvLab)[order(unique(Delta.lrv))])
    ) %>%
  pivot_longer(c(Go, `Discordant Go`),names_to = 'Decision',values_to = 'Probability') %>%
  ggplot(aes(x=Subject,y=Probability,color = Decision)) +
  geom_line()+
  geom_point()+
  scale_color_manual(values = c('red','green'))+
  facet_grid(lrvLab ~ DifferenceLab)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  ggtitle(
    'Probability of Decision vs Interim Sample Size',
    subtitle = 'grouped by true response rate and LRV'
  )+
  ylab('Probability')+
  xlab('Observations at Interim Decision')
