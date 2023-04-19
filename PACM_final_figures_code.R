library(tidyverse)
library (ggplot2)
library(reshape)
#####CONSTANTS####


#number of tasks m

m<-2

#Probability exponent
eta<-10 # in paper range from 1 to 30
#all tasks same demand rate
demandRate<-0.6 
# if all demand differtn stimulusRate<-vector(mode="numeric", length=m)

#reduction in stimulus if all individuals 
#in colony do the task
alpha<-50

#quit probability
tau<-0.2

#### FUNCTIONS####

#stimulus update function
stimFunc<-function(stimVec,activityMat,demandRate,alpha,n){
  output<-stimVec+demandRate-alpha*(colSums(activityMat)/n)
  output[output<0]<-0
  return(output)
}

probabilityfunc<-function(threshMat,stimVec,n,m,eta){
  probMat<-matrix(data = 0, nrow=n,ncol=m)
  for (i in 1:m) {
    s<-stimVec[i]
    denominator<-s^eta+threshMat[,i]^eta
    P<-(s^eta)/denominator
    probMat[,i]<-P
  }
  return(probMat)
}


pop_thresh_mean<-10
means_total<-2*pop_thresh_mean

ant1_thresh1_vec<-1:(means_total-1)

ant1_thresh2_vec<-means_total-ant1_thresh1_vec


n_vec<-1:30

sd_vec<-seq(from=0,to=0.3,by=0.01)

reps<-20
timesteps<-10000
########################################
roleDeterminismFunc<-function(n,m,threshMean,threshSD,individualThresh,reps,timesteps){
  timeProportion<-matrix(data=0,nrow=reps,ncol=m)
  for (r in 1:reps) {
    thresholdMat<-matrix(data=rnorm(n*m,mean = threshMean,sd=threshSD*threshMean),nrow=n,ncol=m)
    thresholdMat[1,]<-c(individualThresh,2*threshMean-individualThresh)
    
    activityMat<-matrix(data=0,nrow=n,ncol=m)
    
    #active(1) or inactive (0)
    
    activity<-vector(mode="numeric", length=n)
    
    #stimulus
    stimulusVec<-vector(mode="numeric", length=m)
    
    antTrack<-vector(mode="numeric", length=timesteps)
    for (t in 1:timesteps) {
      #update probabilities
      probabilities<-probabilityfunc(threshMat=thresholdMat, stimVec = stimulusVec,n=n,m=m, eta=eta)
      
      #update individuals
      activity<-rowSums(activityMat)
      quits<-rbinom(n=n,size=1,prob=tau)
      for (i in 1:n) {
        #active individuals may become inactive
        if (activity[i]==1&&quits[i]==1) {
          activityMat[i,]<-matrix(data=0,nrow=1,ncol=m)
        }
        #inactive individuals may become active 
        if (activity[i]==0) {
          seq<-sample(1:m,m)
          for (j in 1:m) {
            activity[i]<-sum(activityMat[i,])
            if (activity[i]==0) {
              taskNum<-seq[j]
              taskProb<-probabilities[i,taskNum]
              activityMat[i,taskNum]<-rbinom(n=1,size=1,prob=taskProb)
            }
          }
        }
      }
      
      #change stimulus from activity in the past timestep
      stimulusVec<-stimFunc(stimulusVec,activityMat,demandRate,alpha,n) 
      
      #record data
      antTrack[t]<-sum(activityMat[1,]*1:m)
    }
    task1<-sum(antTrack==1)
    task2<-sum(antTrack==2)
    
    timeProportion[r,]<-c(sum(task1)/timesteps,sum(task2)/timesteps)
  }
  return(timeProportion)
}
###################
roleDeterminismFunc(n=10,m=2,threshMean=10, threshSD=0.1, individualThresh = 7,reps=3,timesteps=100) %>% colMeans()



## vary group size####

task1_mat_mean_groupsize_thresh<-matrix(data=0,nrow=length(ant1_thresh1_vec),ncol=length(n_vec))
colnames(task1_mat_mean_groupsize_thresh)<-n_vec
rownames(task1_mat_mean_groupsize_thresh)<-ant1_thresh1_vec
task1_mat_sd_groupsize_thresh<-matrix(data=0,nrow=length(ant1_thresh1_vec),ncol=length(n_vec))
colnames(task1_mat_sd_groupsize_thresh)<-n_vec
rownames(task1_mat_sd_groupsize_thresh)<-ant1_thresh1_vec


ColSD_func<-function(data){
  apply(data,2, sd)
}


for (i in 1:length(ant1_thresh1_vec)) {
  thresVar<-ant1_thresh1_vec[i]
  antlist<-lapply(X=n_vec, FUN = roleDeterminismFunc,m=2,threshMean=10, threshSD=0.1, individualThresh = thresVar,reps=100,timesteps=10000)
  
  #time means
  timeMeans<-lapply(X=antlist, FUN=colMeans)
  timeMeans<-matrix(unlist(timeMeans), ncol=2,byrow=TRUE)
  task1_time<-(timeMeans[,1])/(timeMeans[,1]+timeMeans[,2])
  task1_mat_mean_groupsize_thresh[i,]<-task1_time
  
  
  #time sds
  timeSDs<-lapply(X=antlist,FUN=ColSD_func)
  timeSDs<-matrix(unlist(timeSDs), ncol=2,byrow=TRUE)
  task1_sd<-timeSDs[,1]
  task1_mat_sd_groupsize_thresh[i,]<-task1_sd
}

mean_melt<-melt(task1_mat_mean_groupsize_thresh)
colnames(mean_melt)<-c("Task 1 threshold","Colony Size","Proportion of Active Time on Task 1")
ggplot(mean_melt,aes(X2,X1))+
  geom_tile(aes(fill=value))

SD_melt<-melt(task1_mat_sd_groupsize_thresh)
colnames(SD_melt)<-c("Task 1 threshold","Colony Size","Std Deviation of Active Time on Task 1")

ggplot(SD_melt,aes(X2,X1))+
  geom_tile(aes(fill=value))



#########change group SD ##########


task1_mat_mean_sd_thresh<-matrix(data=0,nrow=length(ant1_thresh1_vec),ncol=length(sd_vec))
colnames(task1_mat_mean_sd_thresh)<-sd_vec
rownames(task1_mat_mean_sd_thresh)<-ant1_thresh1_vec
task1_mat_sd_sd_thresh<-matrix(data=0,nrow=length(ant1_thresh1_vec),ncol=length(sd_vec))
colnames(task1_mat_sd_sd_thresh)<-sd_vec
rownames(task1_mat_sd_sd_thresh)<-ant1_thresh1_vec


ColSD_func<-function(data){
  apply(data,2, sd)
}


for (i in 1:length(ant1_thresh1_vec)) {
  thresVar<-ant1_thresh1_vec[i]
  antlist<-lapply(X=sd_vec, FUN = roleDeterminismFunc,m=2,threshMean=10, n=10, individualThresh = thresVar,reps=100,timesteps=10000)
  
  #time means
  timeMeans<-lapply(X=antlist, FUN=colMeans)
  timeMeans<-matrix(unlist(timeMeans), ncol=2,byrow=TRUE)
  task1_time<-(timeMeans[,1])/(timeMeans[,1]+timeMeans[,2])
  task1_mat_mean_sd_thresh[i,]<-task1_time
  
  
  #time sds
  timeSDs<-lapply(X=antlist,FUN=ColSD_func)
  timeSDs<-matrix(unlist(timeSDs), ncol=2,byrow=TRUE)
  task1_sd<-timeSDs[,1]
  task1_mat_sd_sd_thresh[i,]<-task1_sd
}

mean_melt<-melt(task1_mat_mean_sd_thresh)

colnames(mean_melt)<-c("Task 1 threshold","Colony Threshold Variation","Proportion of Active Time on Task 1")
ggplot(mean_melt,aes(y=`Task 1 threshold`,x=`Colony Threshold Variation`))+
  geom_tile(aes(fill=`Proportion of Active Time on Task 1`))

SD_melt<-melt(task1_mat_sd_sd_thresh)
colnames(SD_melt)<-c("Task 1 threshold","Colony Threshold Variation","Std Deviation of Active Time on Task 1")

ggplot(SD_melt,aes(y=`Task 1 threshold`,x=`Colony Threshold Variation`))+
  geom_tile(aes(fill=`Std Deviation of Active Time on Task 1`))

