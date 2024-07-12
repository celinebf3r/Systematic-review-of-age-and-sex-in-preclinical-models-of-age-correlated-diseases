samplesize<-function(D,k,DELTA){
  K<-c(5,11,17,23,29)
  DISEASE<-c("COPD","Heart disease","Alzheimer","Stroke","Diabetis melitus")
  A1<-as.numeric(D[,K[k]])

  B1<-function(n){sample(A1,size=n,replace=TRUE)}
  B2<-function(n){sample(A1+DELTA,size=n,replace=TRUE)}

  boot.sim.size1<-sim.ssize.wilcox.test(rx = B1, ry = B2, n.max = 300,
                                      iter = 1000,sig.level = 0.05/5,alternative = "less")
  i.p<-min(which(boot.sim.size1$emp.power>0.8))
  if(i.p==Inf){print("n.max>300")}else{
  n.max<-boot.sim.size1$n[i.p]
  n.min<-n.max-10
  }
  boot.sim.size2<-sim.ssize.wilcox.test(rx = B1, ry = B2, n.max = n.max,n.min=n.min,step.size=1,
                                        iter = 10000,sig.level = 0.05/5,alternative = "less")
  
  i.p<-min(which(boot.sim.size2$emp.power>0.8))
  n<-boot.sim.size2$n[i.p]
  
 # print(paste("number of articles needed for ",DISEASE[k],n,"DELTA=",DELTA))
  return(n)
}

samplepower<-function(D,k,DELTA,n){
  K<-c(5,11,17,23,29)
  DISEASE<-c("COPD","IHD","Alzheimer’s disease","stroke","diabetes mellitus type 2")
  A1<-as.numeric(D[,K[k]])
  
  B1<-function(n){sample(A1,size=n,replace=TRUE)}
  B2<-function(n){sample(A1+DELTA,size=n,replace=TRUE)}
  
  res<-sim.power.wilcox.test(nx=n,rx = B1, ny=n,ry = B2,
                                        iter = 100000,sig.level = 0.05/5,alternative = "less")
  
  
 # print(paste("power for ",DISEASE[k],"(effect size in weeks=",DELTA,")","sample size=",n))
  POWER<-length(which(res$Asymptotic$H1$pvalue<0.01))/length(res$Asymptotic$H1$pvalue)
#  print(paste("==>POWER",POWER))
  return(POWER)
}

samplemode<-function(D,k,n,delta=0,alpha=0.05){
  K<-c(5,11,17,23,29)
  DISEASE<-c("COPD","IHD","Alzheimer’s disease","stroke","diabetes mellitus type 2")
  A1<-as.numeric(D[,K[k]])
  
  B1<-function(n,delta){sample(A1+delta,size=n,replace=TRUE)}
  x<-B1(n,0)
  y<-B1(n,delta)
  
  modres1<-modetest(x)
  modres2<-modetest(y)
  if(modres1$p.value<alpha|modres2$p.value<alpha){
    
    i.vert<-which.max(abs(diff(density(y)$y[which(diff(sign(diff(density(y)$y)))==-2)+1])))
    i1<-(which(diff(sign(diff(density(y)$y)))==-2)+1)[i.vert]
    i2<-(which(diff(sign(diff(density(y)$y)))==-2)+1)[i.vert+1]
    THRESH<-density(y)$x[which.min(density(y)$y[i1:i2])+i1]
  
  
  p1<-length(which(x<THRESH))/n
  p2<-length(which(x>THRESH))/n
  p12<-length(which(y<THRESH))/n
  p22<-length(which(y>THRESH))/n
  n1<-length(which(x<THRESH))
  n2<-length(which(x>THRESH))
  n12<-length(which(y<THRESH))
  n22<-length(which(y>THRESH))
  n_new<-power.prop.test(p1=p2,p2=p2+0.2,sig.level = 0.01,power=0.8)$n
  if(n_new<n){
    higherg<-c(n2,n22)
    PV<-prop.test(higherg,c(n,n))$p.value
  }else{
    xx<-B1(n_new-n,delta)
    yy<-B1(n_new-n,delta)
    x<-c(x,xx)
    y<-c(y,yy)
   
    n1<-length(which(x<THRESH))
    n2<-length(which(x>THRESH))
    n12<-length(which(y<THRESH))
    n22<-length(which(y>THRESH))
    higherg<-c(n2,n22)
    PV<-prop.test(higherg,c(n_new,n_new))$p.value
  }
  }else{
    PV<-wilcox.test(x,y)$p.value
  }
  
  return(PV)
}

meanmeta<-function(mu,sigma1,sigma2,n){
  D<-data.frame()
  for(i in 1:n){
    n1<-sample(5:70,1)
    mu1<-rnorm(1,0,sigma1)
    mu2<-rnorm(n1,mu1,sigma2)
    mu2<-mu2+mu
    D<-rbind(D,data.frame(mu=mu2,metastudy=i,n=n1))
  }
  
  return(D)
}