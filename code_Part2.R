set.seed(42)
n = 1000

#Expected values
labda1 <- 16.95 # Part1
mu1 <- 0.429
sd1 <- 0.097
shapeArr <- 7.25
rateArr <- 8.06

shapedur <- 11.39
ratedur <- 16.78

#95% Pessimistic
#labda1 <- 18.81
#mu1 <- 0.439
#sd1 <- 0.104
#shapeArr <- 8.98
#rateArr <- 6.84
#shapedur <- 13.81
#ratedur <- 14.09


gridsize1 <- 12 # We work in scheduling appointments per 5 minutes 12
gridsize2 <- 18 #18
score1 <- matrix(data = 0, nrow = gridsize1+1, ncol = gridsize2+1)
score2 <- matrix(data = 0, nrow = gridsize1+1, ncol = gridsize2+1)

otsum1 <- matrix(data = 0, nrow = gridsize1+1, ncol = gridsize2+1)
otsum2 <- matrix(data = 0, nrow = gridsize1+1, ncol = gridsize2+1)
fsum1 <-  matrix(data = 0, nrow = gridsize1+1, ncol = gridsize2+1)
fsum2 <- matrix(data = 0, nrow = gridsize1+1, ncol = gridsize2+1)
dsum1 <- matrix(data = 0, nrow = gridsize1+1, ncol = gridsize2+1)
dsum2 <- matrix(data = 0, nrow = gridsize1+1, ncol = gridsize2+1)
#Loops
for(k in 1:n){
      
  #DataGeneratingProcesses
  maxT <- 9 # 8:00-17:00 is 9h
  jobsType1times <- c()
  t1 <- 0
      
  while(t1 < maxT){
    q <- rexp(1, rate = labda1/9) # Rate is labda per h from poisson in 1.
    t1 <- t1 + q
        
    if(t1 < maxT){
      jobsType1times <- c(jobsType1times,t1)
    }
  }
  
  t2 <- 0
  jobsType2times <- c()
  while(t2 < maxT){
    q <- rgamma(1,shapeArr, rate = rateArr) # Implement distribution from part1
    t2 <- t2 + q
    if(t2 < maxT){
      jobsType2times <- c(jobsType2times,t2)
    }
  }
  jobsType1lengths <- rnorm(length(jobsType1times), mu1, sd = sd1) 
  jobsType2lengths <- rgamma(length(jobsType2times),shapedur, rate = ratedur) #placeholder
  #GRIDSEARCH
  for(i in 0:gridsize1 ){
    scheduledTime1 <- i/12
    for(j in 0:gridsize2 ){
      scheduledTime2 <- j/12
      
      ot11 = 0
      ot12 = 0
      ot21 = 0
      ot22 = 0
      delay12 = 0
      delay11 = 0
      delay21 = 0
      delay22 = 0
      
      # OLD METHOD
      schedule11 <- c()
      for(l in 1:length(jobsType1times)){
        schedule11 <- c(schedule11,(l-1)*scheduledTime1)
      }
      schedule12 <- c()
      for(l in 1:length(jobsType2times)){
        schedule12 <- c(schedule12,(l-1)*scheduledTime2)
      }
      
      
      if(length(jobsType1lengths) > 1){
        trueschedule11 <- c(0)
      for(l in 2:length(jobsType1times)){
        trueschedule11 <- c(trueschedule11,max((trueschedule11[l-1]+jobsType1lengths[l-1]),schedule11[l]))
      }
      }

      trueschedule12 <- c(0)
      if(length(jobsType2lengths) > 1){
      for(l in 2:length(jobsType2times)){
        trueschedule12 <- c(trueschedule12,max((trueschedule12[l-1]+jobsType2lengths[l-1]), schedule12[l]))
      }
      }

      
      if(length(jobsType1lengths)>0){
        delay11 <- trueschedule11 - schedule11
        
        for(l in 1:length(delay11)){
          if(delay11[l] < 0){
            delay11[l] <- 0
          }
        }
        delay11 <- sum(delay11)
        ot11 <- max(0,trueschedule11[length(jobsType1lengths)] + jobsType1lengths[length(jobsType1lengths)] -9)
      }
      else{
        ot11 <- 0
        delay11 <- 0
      }
      if(length(jobsType2lengths)>0){
        delay12 <- trueschedule12 - schedule12
        for(l in 1:length(delay12)){
          if(delay12[l] < 0){
            delay12[l] <- 0
          }
        }
        delay12 <- sum(delay12)
        ot12 <- max(0,trueschedule12[length(jobsType2lengths)] + jobsType2lengths[length(jobsType2lengths)] -9)
          
      }
      else{
        ot12 <- 0
        delay12<- 0
      }
      infeasible1 <- 0
      if(ot11 > 1 || ot12 > 1){
        infeasible1 <- 1
      }
      
      
      #NEW METHOD
      
      jobsM1 <- c()
      jobsM2 <- c()
      TM1 <- 0
      TM2 <- 0
      typeM2 <- c()
      typeM1 <- c()
      
      schedule21 <- c()
      schedule22 <- c()
      lastend21 <- 0
      lastend22 <- 0
      
      c1 <- 1
      c2 <- 1
      stop <- F
      while(!stop){
        bestT <- Inf
        bestC <- 0
        if(c1 <= length(jobsType1lengths)){
          bestT <- jobsType1times[c1]
          bestC <- 1
        }
        if(c2 <= length(jobsType2lengths) && jobsType2times[c2] < bestT){
          bestC <- 2
        }
        
        bestM <- 1
        if(lastend22 < lastend21){
          bestM <- 2
        }
        
        if(bestM == 1){
          if(bestC == 1){
            jobsM1 <- c(jobsM1, c1)
            c1 <- c1+1
            schedule21 <- c(schedule21,lastend21)
            lastend21 <- lastend21 + scheduledTime1
            typeM1 <- c(typeM1,1)
          }
          if(bestC == 2){
            jobsM1 <- c(jobsM1, c2)
            c2 <- c2+1
            schedule21 <- c(schedule21,lastend21)
            lastend21 <- lastend21 + scheduledTime2
            typeM1 <- c(typeM1,2)
          }
        }
        if(bestM == 2){
          if(bestC == 1){
            jobsM2 <- c(jobsM2, c1)
            c1 <- c1+1
            schedule22 <- c(schedule22,lastend22)
            lastend22 <- lastend22 + scheduledTime1
            typeM2 <- c(typeM2,1)
          }
          if(bestC == 2){
            jobsM2 <- c(jobsM2, c2)
            c2 <- c2+1
            schedule22 <- c(schedule22,lastend22)
            lastend22 <- lastend22 + scheduledTime2
            typeM2 <- c(typeM2,2)
          }
        }
        if(c1 > length(jobsType1lengths) && c2 > length(jobsType2lengths)){
          stop = T
        }
      }
      lastTrueend21 <- 0
      trueschedule21 <- c()
      if(length(schedule21) > 0){
        for(l in 1:length(schedule21)){
          trueschedule21 <- c(trueschedule21,lastTrueend21)
          if(typeM1[l] == 1){
            lastTrueend21 <- max(schedule21[l+1],lastTrueend21 + jobsType1lengths[jobsM1[l]])
          }
          if(typeM1[l] == 2){
            lastTrueend21 <- max(schedule21[l+1],lastTrueend21 + jobsType2lengths[jobsM1[l]])
          }
        }
      }
      lastTrueend22 <- 0
      trueschedule22 <- c()
      if(length(schedule22) > 0){
        for(l in 1:length(schedule22)){
          trueschedule22 <- c(trueschedule22,lastTrueend22)
          if(typeM2[l] == 1){
            lastTrueend22 <- max(schedule22[l+1],lastTrueend22 + jobsType1lengths[jobsM2[l]])
          }
          if(typeM2[l] == 2){
            lastTrueend22 <- max(schedule22[l+1],lastTrueend22 + jobsType2lengths[jobsM2[l]])
          }
        }
      }
      
      
      
      if(length(typeM1)>0){
        delay21 <- trueschedule21 - schedule21

        for(l in 1:length(delay21)){
          if(delay21[l] < 0){
            delay21[l] <- 0
          }
        }
        delay21 <- sum(delay21)
        if(typeM1[length(typeM1)] == 1){
        ot21 <- max(0,trueschedule21[length(typeM1)] + jobsType1lengths[jobsM1[length(typeM1)]] -9)
        }
        if(typeM1[length(typeM1)] == 2){
          ot21 <- max(0,trueschedule21[length(typeM1)] + jobsType2lengths[jobsM1[length(typeM1)]] -9)
        }
      }
      else{
        ot21 <- 0
        delay21 <- 0
      }
      if(length(typeM2)>0){
        delay22 <- trueschedule22 - schedule22
        
        for(l in 1:length(delay22)){
          if(delay22[l] < 0){
            delay22[l] <- 0
          }
        }
        delay22 <- sum(delay22)
        if(typeM2[length(typeM2)] == 1){
          ot22 <- max(0,trueschedule22[length(typeM2)] + jobsType1lengths[jobsM2[length(typeM2)]] -9)
        }
        if(typeM2[length(typeM2)] == 2){
          ot22 <- max(0,trueschedule22[length(typeM2)] + jobsType2lengths[jobsM2[length(typeM2)]] -9)
        }
      }
      else{
        ot22 <- 0
        delay22 <- 0
      }
      infeasible2 <- 0
      if(ot21 > 1 || ot22 > 1){
        infeasible2 <- 1
        
      }
      
      score2[i+1,j+1] <- score2[i+1,j+1] + (10*infeasible2 + (delay21 + delay22)/(length(jobsType1lengths)+length(jobsType2lengths)) + (ot21 + ot22))/n
      score1[i+1,j+1] <- score1[i+1,j+1] + (10*infeasible1 + (delay11 + delay12)/(length(jobsType1lengths)+length(jobsType2lengths)) + (ot11 + ot12))/n
      
      otsum1[i+1,j+1] <- otsum1[i+1,j+1] + (ot11 + ot12)/n
      otsum2[i+1,j+1] <- otsum2[i+1,j+1] + (ot21 + ot22)/n
      fsum1[i+1,j+1] <- fsum1[i+1,j+1] + infeasible1/n
      fsum2[i+1,j+1] <- fsum2[i+1,j+1] + infeasible2/n
      dsum1[i+1,j+1] <- dsum1[i+1,j+1] + (delay11 + delay12)/(length(jobsType1lengths)+length(jobsType2lengths))/n
      dsum2[i+1,j+1] <- dsum2[i+1,j+1] + (delay12 + delay22)/(length(jobsType1lengths)+length(jobsType2lengths))/n
    }
  }
}








scoremaxi <- c()
for(i in 0:gridsize1){
  scoremaxi <- c(scoremaxi, max(score1[i+1,]))
}
seqi <- seq(from = 0, to = gridsize1, by = 1)
plot(seqi,scoremaxi, type = "b")

scoremaxj <- c()
for(j in 0:gridsize2){
  scoremaxj <- c(scoremaxj, max(score1[,j+1]))
}
seqj <- seq(from = 0, to = gridsize2, by = 1)
plot(seqj,scoremaxj, type = "b")




besti <- 0
bestj <- 0
bestscore <- Inf
bestscore2 <- Inf
bestij2 <- c(0,0)
for(i in 0:gridsize1){
  for(j in 0:gridsize2){
    if(score1[i+1,j+1] < bestscore){
      bestscore <- score1[i+1,j+1]
      besti <- i
      bestj <- j
    }
    if(score2[i+1,j+1] < bestscore2){
      bestscore2 <- score2[i+1,j+1]
      bestij2 <- c(i,j)
    }
  }
}
print(cbind(besti,bestj)*5)
print(bestij2*5)
print(bestscore)
print(bestscore2)

print(dsum1[besti+1,bestj+1]*60)
print(dsum2[bestij2[1]+1,bestij2[2]+1]*60)
print(fsum1[besti+1,bestj+1])
print(fsum2[bestij2[1]+1,bestij2[2]+1])
print(otsum1[besti+1,bestj+1]*60)
print(otsum2[bestij2[1]+1,bestij2[2]+1]*60)
