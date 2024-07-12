
meanmeta<-function(mu,sigma1,sigma2,n,nmin,nmax){
  D<-data.frame()
  for(i in 1:n){
    n1<-sample(nmin:nmax,1)
    mu1<-rnorm(1,0,sigma1)
    mu2<-rnorm(n1,mu1,sigma2)
    mu2<-mu2+mu
    D<-rbind(D,data.frame(mu=mu2,metastudy=i,n=n1))
  }
  
  return(D)
}

library(nlme)
X<-read.csv("../data/Cochrane/alzheimer.csv",sep=";",dec=",")
D<-X[,1:12]
D$Alter[which(D$Alter=="not reported")]<-NA
D$Alter<-as.numeric(gsub(",",".",D$Alter))
DDa<-D[-which(is.na(D$Alter)),]
lmres<-lm(Alter~1,data=DDa)
lmresA<-lme(Alter~1,random=~1|Metastudie, data=DDa)
if(AIC(lmresA)>AIC(lmres))
  sigma1<-intervals(lmresA)$re$M[2][1,1] else sigma1=0
mu1<-lmresA$coefficients$fixed
sigma2<-lmresA$sigma
sigma1<-intervals(lmresA)$re$M[2][1,1]
Nmeta<-length(unique(DDa$Metastudie))
Nstudies<-unique(DDa$Anzahl.Studien)

print("################################")
print("Alzheimer")
print(paste("mean age",mu1))
print(paste("standard deviation of age within a metastudy sigma2",sigma2))
print(paste("standard deviation of age between metastudies sigma1",sigma1)) 
print(paste("Numbers of studies vary between",min(Nstudies),"and",max(Nstudies)))
print("################################")

################################
set.seed(232313)
CI<-data.frame()

N<-1000
L<-vector()
for(i in 1:N){
  D<-meanmeta(mu=mu1,sigma1=sigma1,sigma2=sigma2,n=10,nmin = min(Nstudies), nmax = max(Nstudies))
  CItmp<-intervals(lme(mu~1,random=~1|metastudy,data=D),level=0.95,which="fixed")$fixed
  CI<-rbind(CI,CItmp)
}
######################
print(paste("mean of the median",mean(CI[,2])))
print(paste("Mean of the lower bound",mean(CI[,1]),"mean of the upper bound",mean(CI[,3])))
print(paste("Mean of the difference between upper and lower bound ",mean(CI[,3]-CI[,1])))

###########################################################
X<-read.csv("../data//Cochrane/diabetis.csv",sep=";")

D<-X[,1:13]
D$Alter[which(D$Alter=="not reported")]<-NA
D$Alter<-as.numeric(gsub(",",".",D$Alter))
I.na<-which(is.na(D$Anzahl.Studien))
DDa<-D[-c(which(is.na(D$Alter)),I.na),]
DDa<-DDa[-106,]
lmres<-lm(Alter~1,data=DDa)
lmresA<-lme(Alter~1,random=~1|Metastudie, data=DDa)
if(AIC(lmresA)>AIC(lmres))
  sigma1<-intervals(lmresA)$re$M[2][1,1] else sigma1=0
mu1<-lmresA$coefficients$fixed
sigma2<-lmresA$sigma
sigma1<-intervals(lmresA)$re$M[2][1,1]
Nmeta<-length(unique(DDa$Metastudie))
Nstudies<-unique(DDa$Anzahl.Studien)

print("################################")
print("Diabetis")
print(paste("mean age",mu1))
print(paste("standard deviation of age within a metastudy sigma2",sigma2))
print(paste("standard deviation of age between metastudies sigma1",sigma1)) 
print(paste("Numbers of studies vary between",min(Nstudies),"and",max(Nstudies)))
print("################################")

set.seed(232313)
CI<-data.frame()

N<-1000
L<-vector()
for(i in 1:N){
  D<-meanmeta(mu=mu1,sigma1=sigma1,sigma2=sigma2,n=10,nmin = min(Nstudies), nmax = max(Nstudies))
  CItmp<-intervals(lme(mu~1,random=~1|metastudy,data=D),level=0.95,which="fixed")$fixed
  CI<-rbind(CI,CItmp)
}
######################
print(paste("mean of the median",mean(CI[,2])))
print(paste("Mean of the lower bound",mean(CI[,1]),"mean of the upper bound",mean(CI[,3])))
print(paste("Mean of the difference between upper and lower bound ",mean(CI[,3]-CI[,1])))
############################################################
X<-read.csv("../data/Cochrane/stroke.csv",sep=";")

D<-X[,1:13]
D$Metastudie<-as.factor(D$Metastudie)
D$Alter[which(D$Alter=="not reported")]<-NA
D$Alter<-as.numeric(gsub(",",".",D$Alter))
DDa<-D[-which(is.na(D$Alter)),]
lmres<-lm(Alter~1,data=DDa)
lmresA<-lme(Alter~1,random=~1|Metastudie, data=DDa)
mu1<-lmresA$coefficients$fixed
sigma2<-lmresA$sigma
if(AIC(lmresA)>AIC(lmres))
  sigma1<-intervals(lmresA)$re$M[2][1,1] else sigma1=0
Nmeta<-length(unique(DDa$Metastudie))
Nstudies<-unique(DDa$Anzahl.Studien)

print("################################")
print("Stroke")
print(paste("mean age",mu1))
print(paste("standard deviation of age within a metastudy sigma2",sigma2))
print(paste("standard deviation of age between metastudies sigma1",sigma1)) 
print(paste("Numbers of studies vary between",min(Nstudies),"and",max(Nstudies)))
print("################################")
set.seed(232313)
CI<-data.frame()

N<-1000
L<-vector()
for(i in 1:N){
  D<-meanmeta(mu=mu1,sigma1=sigma1,sigma2=sigma2,n=10,nmin = min(Nstudies), nmax = max(Nstudies))
  CItmp<-intervals(lme(mu~1,random=~1|metastudy,data=D),level=0.95,which="fixed")$fixed
  CI<-rbind(CI,CItmp)
}
######################
print(paste("mean of the median",mean(CI[,2])))
print(paste("Mean of the lower bound",mean(CI[,1]),"mean of the upper bound",mean(CI[,3])))
print(paste("Mean of the difference between upper and lower bound ",mean(CI[,3]-CI[,1])))



