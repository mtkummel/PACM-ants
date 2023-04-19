#####CONSTANTS####

ColSD_func<-function(data){
  apply(data,2, sd)
}


# number of individuals n

n<-30

#number of tasks m

m<-2

#Probability exponent
eta<-7 # in paper range from 1 to 30
#all tasks same demand rate
demandRate<-0.6 
# if all demand differtn stimulusRate<-vector(mode="numeric", length=m)

#reduction in stimulus if all individuals 
#in colony do the task
alpha<-50

#quit probability
tau<-0.2

corrStep       <- 200 #number of time steps for calculation of correlation 


threshMean<-10 #forr both tasks
threshSD<-0.1
#time
timesteps<-10000
ThreshM        <- c(10, 10) #population threshold means 
ThreshSD       <- ThreshM * 0 #population threshold standard deviations 
InitialStim    <- c(0, 0) #intital vector of stimuli
StimRates      <- c(0.6, 0.6) #vector of stimuli increase rates  
threshSlope    <- 30 #exponent parameter for threshold curve shape  
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active


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

############################### #################
nvec<-c(2,4,6,8,12,16)

DataCollectionMat<-matrix(data=0,nrow=length(nvec),ncol=3)
colnames(DataCollectionMat)<-c("n","Specialization","sem")
reps<-100
replicateMat<-matrix(data=0,nrow = reps,ncol=2)
m=2
threshMean=10 
threshSD=0.1

timesteps=10000

for (q in 1:length(nvec)) {
  
n=nvec[q]

for (r in 1:reps) {
  

thresholdMat<-matrix(data=rnorm(n*m,mean = threshMean,sd=threshSD*threshMean),nrow=n,ncol=m)


activityMat<-matrix(data=0,nrow=n,ncol=m)

AntTaskTrack_1<-matrix(data=0,nrow=timesteps,ncol=n)
AntTaskTrack_2<-matrix(data=0,nrow=timesteps,ncol=n)


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
  AntTaskTrack_1[t,]<-activityMat[,1]
  AntTaskTrack_2[t,]<-activityMat[,2]
}

windownum<-floor(timesteps)/200
spear_corr_mat<-matrix(data = 0,nrow=(windownum-2),ncol=2)


for (i in 1:(windownum-2) ) {
  window_T1_1<-AntTaskTrack_1[(200*i):(200*i+200),]
  window_T2_1<-AntTaskTrack_2[(200*i):(200*i+200),]
  task1_active_1<-colSums(window_T1_1)
  task2_active_1<-colSums(window_T2_1)
  
  window_T1_2<-AntTaskTrack_1[(200*i+200):(200*i+400),]
  window_T2_2<-AntTaskTrack_2[(200*i+200):(200*i+400),]
  task1_active_2<-colSums(window_T1_2)
  task2_active_2<-colSums(window_T2_2)
  
  spear_corr_1<-cor.test(task1_active_1,task1_active_2,method = 'spearman')
  spear_corr_2<-cor.test(task2_active_1,task2_active_2,method = 'spearman')
  spear_corr_mat[i,]<-c(spear_corr_1$estimate,spear_corr_2$estimate)
}
replicateMat[r,]<-c(mean(colMeans(spear_corr_mat, na.rm=T)),sd(spear_corr_mat,na.rm=T))
}
DataCollectionMat[q,]<-c(n, mean(replicateMat[,1], na.rm=T), mean(replicateMat[,2], na.rm=T)/sqrt(length(spear_corr_mat)))
}
DataCollectionMat<-DataCollectionMat %>% as.data.frame()
ggplot(DataCollectionMat,aes(x=n,y=Specialization))+
  geom_line(color="darkorange")+
  geom_point(color="darkorange")+
  ylim(0,0.8)+
  geom_pointrange(aes(ymin=Specialization-sem,ymax=Specialization+sem),color="orange3")+
  ggtitle("Group specialization as a function of colony size\nThreshold SD=0.1")+
  xlab("Number of ants in colony") +
  ylab("Specialization")
