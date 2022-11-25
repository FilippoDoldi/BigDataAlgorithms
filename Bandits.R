#Setup
rm(list = ls())
dir.work <- "./"
setwd(dir.work)
library(tidyverse) 

sdtart <- -1L # setting starting seed as 
nArmsEachGame <- 4L # to be set
numTrials <- 10000L # to be set

# name strategies = name objects
ArmEpsGreedy <- "ArmEpsGr"
ArmUCB1 <- "ArmUpperConfidence"
ArmThompson <- "ArmThompsonSampling"

#Modified function pullarm
  pullArm <- function(namePull = "Pull1",nArm = 1L,seed_start = 0L, 
                      nPlayGame = 1L, Complications=FALSE){
    if(!exists(namePull,envir = .GlobalEnv)){
      assign(namePull, list(sequence = tibble(numArm = integer(),
                                              rewards = numeric()),
                            stats = tibble(
                              numArm = integer(),
                              meanArmSoFar = numeric(),
                              numPlayed = integer()
                            )) , envir = .GlobalEnv)
    }
    
    if(nPlayGame>0){
      Pull <- get(namePull,envir = .GlobalEnv)
      HistTN <- filter(Pull$stats,numArm == nArm)
      
      if(HistTN %>% nrow()==0 ){
        Pull$stats <- add_row(Pull$stats,numArm=nArm,
                              meanArmSoFar=0.0,numPlayed=0L)
        HistTN <- filter(Pull$stats,numArm == nArm)
      }
      
      # mean rewards: .4, .5, .5, .6, .75, .75, .7, .6, .5, 
      #      .9, (20/11)*1/2, (20/12)*1/2, ...
      
      if(nArm>10){
        a <- 1
        b <- nArm/10
      } else {
        switch (nArm,
                {a=8
                b=12},
                {a=1
                b=1},
                {a=0.05
                b=0.05},
                {a=.6
                b=.4},
                {a=.75
                b=.25},
                {a=3
                b=1},
                {a=.7
                b=.3},
                {a=1.2
                b=.8},
                {a=.5
                b=.5},
                {a=1.8
                b=.2}
        )
      }
      
      set.seed(seed = 1e7*nArm + HistTN$numPlayed + seed_start)
      rewards <- rbeta(n = nPlayGame,shape1 = a,shape2 = b)
      
      #If TRUE we take the arms as bernoulli with success 
      # if there aren't complications
      #If FALSE (default) we do nothing
      if(Complications)
      {
        if(rewards<0.2) rewards<-0
        else rewards<-1
      }
      
      Pull$sequence <- bind_rows(Pull$sequence,
                                 tibble(numArm = nArm,
                                        rewards = rewards))
      
      sumRewSofar <- sum(rewards,HistTN$meanArmSoFar*HistTN$numPlayed)
      
      #  Pull
      
      Pull$stats <- Pull$stats %>% mutate(
        numPlayed = case_when(
          numArm == nArm ~ numPlayed + nPlayGame,
          TRUE ~ numPlayed + 0),
        meanArmSoFar = case_when(
          numArm == nArm ~ sumRewSofar/numPlayed,
          TRUE ~ meanArmSoFar))
      
      assign(namePull, Pull , envir = .GlobalEnv)
      return(rewards)
    }
    
    return()
  }
##########################################################################
#PART 1

if(exists(ArmEpsGreedy)){
  rm(list = c(ArmEpsGreedy))
}
if(exists(ArmUCB1)){
  rm(list = c(ArmUCB1))
}
if(exists(ArmThompson)){
  rm(list = c(ArmThompson))
}

#UCB1 Starting Parameters
muHat<-rep(0,nArmsEachGame) # estimate of payoff per arm so far
m<-rep(0,nArmsEachGame)     # number of pulls per arm so far
UCB<-rep(0,nArmsEachGame)   # estimated upper confidence bound per arm

#Thompson Starting Parameters
success<-rep(0,nArmsEachGame)  #number of successes per arm
failure<-rep(0,nArmsEachGame)  #number of failures per arm
theta<-rep(0,nArmsEachGame)    #will contain the beta distribution

#Number of complications
ComplGreedy_Before<-0
ComplUCB1_Before<-0
ComplThompson_Before<-0

# initialization: one bandit for each arm
for(j in seq_len(nArmsEachGame)){
  
  #EpsGreedy
  tmp <- pullArm(namePull = ArmEpsGreedy,
                 nArm = j,
                 nPlayGame = 1,
                 seed_start = sdtart)
  #if complication count it
  if(tmp<0.2) ComplGreedy_Before<-ComplGreedy_Before+1
  
  #UCB1
  tmp <- pullArm(namePull = ArmUCB1,
                 nArm = j,
                 nPlayGame = 1,
                 seed_start = sdtart)
  muHat[j] <- tmp
  m[j] <- m[j] +1
  #if complication count it
  if(tmp<0.2) ComplUCB1_Before<-ComplUCB1_Before+1
  
  #Thompson
  tmp <- pullArm(namePull = ArmThompson,
                 nArm = j,
                 nPlayGame = 1,
                 seed_start = sdtart)
  # we compute a bernoulli with probabilty of success equal the reward
  # if it's 1 it's a success else is a failure (for the pulled arm)
  if(rbernoulli(1,tmp)==TRUE) success[j] <-success[j] + 1
  else failure[j] <- failure[j] + 1
  #if complication count it
  if(tmp<0.2) ComplThompson_Before<-ComplThompson_Before+1
}

#Starting point
get(ArmEpsGreedy)#=get(ArmUCB1)=get(ArmThompson)

#Cycle on the trials
for(j in seq_len(numTrials)){
  if(j %% 100 == 0){
    cat("Num trials ",j,"out of",numTrials,"\n")
  }

  # EpsGreedyAlgorithm
  exploration <- runif(1) < 4/j
  if(exploration){
    tmp <- pullArm(namePull = ArmEpsGreedy,
                   nArm = sample.int(nArmsEachGame,size = 1),
                   nPlayGame = 1,
                   seed_start = sdtart) 
  } else {
    bestSoFar <- get(ArmEpsGreedy)$stats %>% 
      filter(meanArmSoFar == max(meanArmSoFar)) 
    
    tmp <- pullArm(namePull = ArmEpsGreedy,
                   nArm = bestSoFar$numArm,
                   nPlayGame = 1,
                   seed_start = sdtart)
  }
  #if complication count it
  if(tmp<0.2) ComplGreedy_Before<-ComplGreedy_Before+1
  
  # UCB1
  # Computing the upper bound for each arm by the formula seen in lesson
  for(i in 1:nArmsEachGame)
    UCB[i] <- muHat[i] + (1/j)*sqrt(2*log(j,exp(1))/m[i])
  # i chose to take the parameter alpha that trades off exploration
  # and exploitation as 1/j with the idea that we grew confident in our 
  # esteem as time grows 
  
  ArmToPull <- which(UCB==max(UCB))
  if(length(ArmToPull)!=1)  ArmToPull<-sample(ArmToPull,1)
  #We randomly choose one if there are more arm to pull by this method
  
  tmp <- pullArm(namePull = ArmUCB1,
                 nArm = ArmToPull,
                 nPlayGame = 1,
                 seed_start = sdtart)
  
  # Updating parameters m and muHat by the formulas seen in lesson
  m[ArmToPull] <- m[ArmToPull] + 1
  muHat[ArmToPull] <- 1/m[ArmToPull]*(get(ArmUCB1)$stats$meanArmSoFar[ArmToPull] 
                                      + (m[ArmToPull] - 1)*muHat[ArmToPull])
  
  #if complication count it
  if(tmp<0.2) ComplUCB1_Before<-ComplUCB1_Before+1
  
  # Thompson Sampling
  # computing the beta distribution for each arm based 
  # on previous successes and failures
  for(i in 1:nArmsEachGame)
       theta[i]<-rbeta(1,1+success[i],1+failure[i])

  ArmToPull <- which(theta==max(theta))
  if(length(ArmToPull)!=1)  ArmToPull<-sample(ArmToPull,1)
  #We randomly choose one if there are more arm to pull by this method
  
  tmp <- pullArm(namePull = ArmThompson,
                nArm = ArmToPull,
                nPlayGame = 1,
                seed_start = sdtart)
  
  # After founding the maximum theta we pull the arm obtaining his reward and 
  # we compute a bernoulli with probabilty of success equal to the reward
  # if it's 1 it's a success else is a failure
  if(rbernoulli(1,tmp)==TRUE) success[ArmToPull] <-success[ArmToPull] + 1
  else failure[ArmToPull] <- failure[ArmToPull] + 1
  
  #if complication count it
  if(tmp<0.2) ComplThompson_Before<-ComplThompson_Before+1
}

# plot time series
seqEpdGreedy <- get(ArmEpsGreedy)$sequence %>%
  rowid_to_column("numGame") %>%
  mutate(meanRewards = cumsum(rewards)/numGame) %>%
  add_column(strategy = ArmEpsGreedy)

seqEpdUCB1 <- get(ArmUCB1)$sequence %>%
  rowid_to_column("numGame") %>%
  mutate(meanRewards = cumsum(rewards)/numGame) %>%
  add_column(strategy = ArmUCB1)

seqEpdThompson <- get(ArmThompson)$sequence %>%
  rowid_to_column("numGame") %>%
  mutate(meanRewards = cumsum(rewards)/numGame) %>%
  add_column(strategy = ArmThompson)


dataTimeSeries <- bind_rows(seqEpdGreedy,seqEpdUCB1,seqEpdThompson) %>%
  mutate(numGame = numGame -nArmsEachGame) %>%
  filter(numGame > 0) 

ggplot(data = dataTimeSeries,mapping = aes(x = numGame,
                                           y = meanRewards)) + 
  geom_line(aes(color = strategy, linetype = strategy)) +
  ylab("mean of Rewards")+
  xlab("number of games")

#Final results
GreedyStats1<-get(ArmEpsGreedy)$stats
UCB1Stats1<-get(ArmUCB1)$stats
ThompsonStats1<-get(ArmThompson)$stats
  #Show them
  GreedyStats1
  UCB1Stats1
  ThompsonStats1
# END PART 1#################################################################

rm(bestSoFar,dataTimeSeries,seqEpdGreedy,seqEpdThompson,seqEpdUCB1)

#PART 2
if(exists(ArmEpsGreedy)){
  rm(list = c(ArmEpsGreedy))
}
if(exists(ArmUCB1)){
  rm(list = c(ArmUCB1))
}
if(exists(ArmThompson)){
  rm(list = c(ArmThompson))
}
#UCB1 Starting Parameters
muHat<-rep(0,nArmsEachGame)
m<-rep(1,nArmsEachGame)
UCB<-rep(0,nArmsEachGame)

#Thompson Starting Parameters
success<-rep(0,nArmsEachGame)
failure<-rep(0,nArmsEachGame)
theta<-rep(0,nArmsEachGame)

#Number of complications in the second case
ComplGreedy_After<-0
ComplUCB1_After<-0
ComplThompson_After<-0

# initialization: one bandit for each arm
for(j in seq_len(nArmsEachGame)){
  tmp <- pullArm(namePull = ArmEpsGreedy,
                 nArm = j,
                 nPlayGame = 1,
                 seed_start = sdtart, Complications = TRUE)
  #if complication count it
  if(tmp<0.2) ComplGreedy_After <- ComplGreedy_After +1
  
  tmp <- pullArm(namePull = ArmUCB1,
                 nArm = j,
                 nPlayGame = 1,
                 seed_start = sdtart, Complications = TRUE)
  muHat[j] <- tmp
  #if complication count it
  if(tmp<0.2) ComplUCB1_After <- ComplUCB1_After +1
  
  tmp <- pullArm(namePull = ArmThompson,
                 nArm = j,
                 nPlayGame = 1,
                 seed_start = sdtart, Complications = TRUE)
  if(rbernoulli(1,tmp)==TRUE) success[j] <-success[j] + 1
  else failure[j] <- failure[j] + 1
  #if complication count it
  if(tmp<0.2) ComplThompson_After <- ComplThompson_After +1
}
#Starting point
get(ArmEpsGreedy)#=get(ArmUCB1)=get(ArmThompson)

for(j in seq_len(numTrials)){
  if(j %% 100 == 0){
    cat("Num trials ",j,"out of",numTrials,"\n")
  }
  
  # EpsGreedyAlgorithm
  exploration <- runif(1) < 4/j
  if(exploration){
    tmp <- pullArm(namePull = ArmEpsGreedy,
                   nArm = sample.int(nArmsEachGame,size = 1),
                   nPlayGame = 1,
                   seed_start = sdtart, Complications = TRUE) 
  } else {
    bestSoFar <- get(ArmEpsGreedy)$stats %>% 
      filter(meanArmSoFar == max(meanArmSoFar))
    
    ArmToPull<-bestSoFar$numArm
    if(length(ArmToPull)!=1) ArmToPull<-sample(ArmToPull,1)
    #We randomly choose one if there are more arm to pull by this method
    
    tmp <- pullArm(namePull = ArmEpsGreedy,
                   nArm = ArmToPull,
                   nPlayGame = 1,
                   seed_start = sdtart, Complications = TRUE)
  }
  #if complication count it
  if(tmp<0.2) ComplGreedy_After <- ComplGreedy_After +1
  
  # UCB1
  for(i in 1:nArmsEachGame)
    UCB[i] <- muHat[i] + (1/j)*sqrt(2*log(j,exp(1))/m[i])
  
  ArmToPull <- which(UCB==max(UCB))
  if(length(ArmToPull)!=1) ArmToPull<-sample(ArmToPull,1)
  #We randomly choose one if there are more arm to pull by this method
  
  tmp <- pullArm(namePull = ArmUCB1,
                 nArm = ArmToPull,
                 nPlayGame = 1,
                 seed_start = sdtart, Complications = TRUE)
  
  m[ArmToPull] <- m[ArmToPull] + 1
  muHat[ArmToPull] <- 1/m[ArmToPull]*(get(ArmUCB1)$stats$meanArmSoFar[ArmToPull] 
                                      + (m[ArmToPull] - 1)*muHat[ArmToPull])
  #if complication count it
  if(tmp<0.2) ComplUCB1_After <- ComplUCB1_After +1
  
  # Thompson Sampling
  for(i in 1:nArmsEachGame)
    theta[i]<-rbeta(1,1+success[i],1+failure[i])
  
  ArmToPull <- which(theta==max(theta))
  if(length(ArmToPull)!=1)  ArmToPull<-sample(ArmToPull,1)
  #We randomly choose one if there are more arm to pull by this method
  
  tmp <- pullArm(namePull = ArmThompson,
                 nArm = ArmToPull,
                 nPlayGame = 1,
                 seed_start = sdtart, Complications = TRUE)
  
  if(tmp) success[ArmToPull] <-success[ArmToPull] + 1 
  if(!tmp) failure[ArmToPull] <-failure[ArmToPull] + 1
  
  #if complication count it
  if(tmp<0.2) ComplThompson_After <- ComplThompson_After +1
}

# plot time series, now we need to show the number of complications vs
# the number of games hence we change the following mutate by counting
# the cumulative complications 
seqEpdGreedy <- get(ArmEpsGreedy)$sequence %>%
  rowid_to_column("numGame") %>%
  mutate(Complications = cumsum(1-rewards)) %>%
  add_column(strategy = ArmEpsGreedy)

seqEpdUCB1 <- get(ArmUCB1)$sequence %>%
  rowid_to_column("numGame") %>%
  mutate(Complications = cumsum(1-rewards)) %>%
  add_column(strategy = ArmUCB1)

seqEpdThompson <- get(ArmThompson)$sequence %>%
  rowid_to_column("numGame") %>%
  mutate(Complications = cumsum(1-rewards)) %>%
  add_column(strategy = ArmThompson)

dataTimeSeries <- bind_rows(seqEpdGreedy,seqEpdUCB1,seqEpdThompson) %>%
  mutate(numGame = numGame -nArmsEachGame) %>%
  filter(numGame > 0) 

ggplot(data = dataTimeSeries,mapping = aes(x = numGame,
                                           y = Complications)) + 
  geom_line(aes(color = strategy, linetype = strategy)) +
  ylab("number of complications")+
  xlab("number of games")

#Final results
GreedyStats2<-get(ArmEpsGreedy)$stats
UCB1Stats2<-get(ArmUCB1)$stats
ThompsonStats2<-get(ArmThompson)$stats
  #Show them
  GreedyStats2
  UCB1Stats2
  ThompsonStats2

#Complications' comparison before and after
paste("Greedy prima:",ComplGreedy_Before,"Greedy dopo:",ComplGreedy_After)
paste("UCB1 prima:",ComplUCB1_Before,"UCB1 dopo:",ComplUCB1_After)
paste("Thompson prima:",ComplThompson_Before,
      "Thompson dopo:",ComplThompson_After)
