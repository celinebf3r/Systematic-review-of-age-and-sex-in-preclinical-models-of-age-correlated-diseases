######## test of multmodes ##################
library(multimode)
set.seed(1234126)
X<-read.csv("../data//preclinical/Pilotstudie_ finale preclinical Tabelle.csv",sep=";",dec=",")
D<-X[3:22,-1]
A1<-D[,5]
A2<-as.numeric(D[,11])
A3<-D[,17]
A4<-D[,23]
A5<-D[,29]
modres1<-modetest(A1,B=10000,BMC=5000)
modres2<-modetest(A2,B=10000,BMC=5000)
modres3<-modetest(A3,B=10000,BMC=5000)
modres4<-modetest(A4,B=10000,BMC=5000)
modres5<-modetest(A5,B=10000,BMC=5000)
par(mfrow=c(2,3))
hist(A1,main = "COPD",xlab="age")
hist(A2,main ="IHD",xlab="age")
hist(A3,main="Alzheimer's disease",xlab="age")
hist(A4,main="Stroke",xlab="age")
hist(A5,main="diabetes mellitus type 2",xlab="age")

print(paste("p value for COPD",modres1$p.value))
print(paste("p value for IHD",modres2$p.value))
print(paste("p value for Alzheimer's disease",modres3$p.value))
print(paste("p value for stroke",modres4$p.value))
print(paste("p value for diabetes mellitus type 2",modres5$p.value))
            