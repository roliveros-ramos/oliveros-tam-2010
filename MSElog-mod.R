##############################################
# MANAGMENT STRATEGY EVALUATION
# Ricardo Oliveros-Ramos @ CIMOBP - IMARPE.
# e-mail: roliveros@imarpe.gob.pe
# Versión 0.9 - 24/06/2009
##############################################
##################################
#### Simulation control ##########
##################################
t.ini=proc.time()
# strategy="E2" 	# to be chosen beetween E1 (min), E2 (precautionary ICES)
# E3 (max viable), E4 (precautionary viable), E5 (viable)
obs=0 				# 1 for observed catches
N=1000 				# number of simulations (Montecarlo)
T=50 # time horizon
spinup=20
colapse=FALSE	# Scenarium
sce=if(colapse) "post" else "pre"
uncertainty=TRUE	# catch uncertainty on/off
uncertainty.b=TRUE # biomass uncertainty
est=c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10")
#est="mor" # moratorium
SSBs=c()
n=1; t=1
##################################
#### Parameters ##################
##################################
# population growth
R=1.5
k=1 # carrying capacity
K=R*k/(R-1)
# Reference points
Bmin=0.25
Bpa=1.5*Bmin
Ymin=0.05 # variar
# Ymin.viab=0.068 # verificar, MVY
# Ymin.viab=Bmin-K*(1-sqrt(1-(4*Bmin/(R*K))))/2
Fpa=round(2*(R-1)/6,2) # 0.33=2F_MSY/3
# uncertainty
sigma.mu=0.1
sigma.sd=sqrt(0.09)
omega.mu=0.05
omega.sd=sqrt(0.04)
delta=0.1
sigma=rlnorm(N*T,mean=log(sigma.mu),sd=sigma.sd)
omega=rlnorm(N*T,mean=log(omega.mu),sd=omega.sd)
delta.b=rnorm(N*T,mean=0,sd=delta)
factor=if(uncertainty) (1+sigma)*(1+omega) else double(N*T)+1
factor=matrix(factor,ncol=T,nrow=N)
factor.b=if(uncertainty.b) delta.b else double(N*T)
factor.b=matrix(factor.b,ncol=T,nrow=N)
##################################
#### Initial conditions ##########
##################################
# Y.health=as.matrix(read.table("yield.txt"))
# Y.colapse=as.matrix(read.table("yield2.txt"))
# Y.obs=if(colapse) Y.colapse else Y.health
# Initial biomasses
B.health=2*Bmin # healthy stock
B.colapse=0.6*Bmin # post-colapse stock
B.ini=if(colapse) B.colapse else B.health
# ENSO
# ENSO=if(colapse) ENSO.colapse else ENSO.health
# T=length(ENSO) # time horizon
times.health=seq(from=0,by=1,len=T+1)
times.colapse=seq(from=0,by=1,len=T+1)
times=if(colapse) times.colapse else times.health
##################################
#### Dynamics ####################
##################################
dynamics.old=function(B,Y) {
  stock=B-Y
  stock=R*stock*(1-stock/K)
  stock=max(stock,0)
  return(stock)
}
dynamics=function(B,Y) {
  stock=R*B*(1-B/K)-Y
  stock=max(stock,0)
  return(stock)
}
##################################
#### Viability results ###########
##################################
#Ymax= function(B) B-(0.5*R*k/(R-1))*(1-sqrt(1-(4*(R-1)*Bmin)/(R*R*k)))
Ymax= function(B) R*B-R*B*B/K-Bmin
test=function(B) {
  TEST=if(B>=Bmin) TRUE else FALSE
  return(TEST)
}
MVY=Ymax(Bmin)
MSY=(R-1)*k/4
FMSY=(R-1)/2
##################################
#### Management strategies #######
##################################
yield=function(B,strategy) {
  Y=if(B<Bmin) 0 else switch(strategy,
                             E1=Ymin,
                             E2=Ymax(B),
                             E3=0.5*(Ymax(B)+Ymin),
                             E4=if(Ymax(B)>=Ymin) runif(1,Ymin,Ymax(B)) else Ymin,
                             E5=max(FMSY*B,Ymin),
                             E6=if(B>=Bpa) Fpa*B else Ymin,
                             E7=if(B>=Bpa) 0.5*(Fpa*B+Ymin) else Ymin,
                             E8=if(B>=Bpa) runif(1,Ymin,Fpa*B) else Ymin,
                             E9=MVY,
                             E10=if(B>=Bpa) 0.5*(Fpa*B+MVY) else Ymin,
                             mor=0)
  return(Y)
}
##################################
#### Analysis ####################
##################################
ref.1=function(x) {
  lim1=c(mean(x),min(x),max(x))
  return(lim1)
}
##################################
#### Simulation ##################
##################################
if(obs==0) {
  B.est=NULL
  Prob.lt=NULL
  Prob.min=NULL
  Prob.pa=NULL
  Prob.mor=NULL
  Prob.mor2=NULL
  B.ind=NULL
  Y.ind=NULL
  B.serie=array(0,dim=c(length(times),length(ref.1(1)),length(est)))
  Y.serie=array(0,dim=c(length(times)-1,length(ref.1(1)),length(est)))
  for(strategy in est) {
    B.t=NULL
    Y.t=NULL
    prob.lt=0
    prob.min=0
    prob.pa=0
    prob.mor=0
    prob.mor2=0
    for(n in 1:N) {
      pop=rbind(B.ini)
      yields=c()
      for(t in 1:T) {
        Y=yield(pop[t,]*(1+factor.b[n,t]),strategy)
        pop.new=dynamics(pop[t,],Y*factor[n,t])
        pop=rbind(pop,pop.new)
        yields=c(yields,Y)
      }
      B.t=cbind(B.t,pop)
      Y.t=cbind(Y.t,yields)
      prob.lt=prob.lt+min(prod(floor(pop[1:spinup]/Bmin)),1)
      prob.min=prob.min+min(prod(floor(pop[-(1:spinup)]/Bmin)),1)
      prob.pa=prob.pa+min(prod(floor(pop[-(1:spinup)]/Bpa)),1)
      prob.mor=prob.mor+sum(yields==0)/length(yields)
      prob.mor2=prob.mor2+sum(yields<=Ymin)/length(yields)
    }
    #########################
    # Scenarium variability #
    #########################
    bio.m=t(apply(B.t,1,ref.1)) # Biomass variability
    y.m=t(apply(Y.t,1,ref.1)) # Yield variability
    B.serie[,,which(est==strategy)]=bio.m
    Y.serie[,,which(est==strategy)]=y.m
    B.est=cbind(B.est,bio.m[,1]) # Spawners mean serie
    #########################
    # Indicators ############
    #########################
    Prob.lt=rbind(Prob.lt,1-prob.lt/N)
    Prob.min=rbind(Prob.min,1-prob.min/N)
    Prob.pa=rbind(Prob.pa,1-prob.pa/N)
    Prob.mor=rbind(Prob.mor,prob.mor/N)
    Prob.mor2=rbind(Prob.mor2,prob.mor2/N)
    B.ind.1=mean(apply(B.t[-(1:spinup),],2,mean))/Bmin
    B.ind.2=mean(apply(B.t[-(1:spinup),],2,sd))/mean(apply(B.t[-(1:spinup),],2,mean))
    B.ind=rbind(B.ind,c(B.ind.1,B.ind.2))
    Y.ind.1=mean(apply(Y.t[-(1:spinup),],2,mean))/Ymin
    Y.ind.2=mean(apply(Y.t[-(1:spinup),],2,sd))/mean(apply(Y.t[-(1:spinup),],2,mean))
    Y.ind=rbind(Y.ind,c(Y.ind.1,Y.ind.2))
  }
  B.est=data.frame(B.est)
  dimnames(B.est)=list(times,est)
  write.table(B.est,file=paste(sce,"F_",Fpa,"R_",R,"_B_mean.csv",sep=""),quote = FALSE, sep = ",", col.names = NA, row.names = TRUE)
  boxplot(B.est[-(1:spinup),])
  ind=cbind(Prob.lt,Prob.min,Prob.pa,Prob.mor,Prob.mor2,B.ind,Y.ind)
  vars=c("B","Y")
  inds=c("Risk.lt","Risk.min","Risk.pa","Prob.mor","Prob.cm")
  for(v in vars) {
    indss=paste(v,".ind.",1:2,sep="")
    inds=c(inds,indss)
  }
  dimnames(ind)=list(est,inds)
  write.table(ind,file=paste(sce,"F_",Fpa,"R_",R,"_indicadores.csv",sep=""),quote = FALSE, sep = ",", col.names = NA, row.names = TRUE)
  names.1=list(as.character(times),c("Mean","Min","Max"),est)
  names.2=list(as.character(times[-length(times)]),c("Mean","Min","Max"),est)
  dimnames(B.serie)=names.1
  dimnames(Y.serie)=names.2
  B.serie=data.frame(B.serie)
  write.table(B.serie,file=paste(sce,"F_",Fpa,"R_",R,"_B.csv",sep=""),quote = FALSE, sep = ",", col.names = NA, row.names = TRUE)
  Y.serie=data.frame(Y.serie)
  write.table(Y.serie,file=paste(sce,"F_",Fpa,"R_",R,"_Y.csv",sep=""),quote = FALSE, sep = ",", col.names = NA, row.names = TRUE)
} # End if MSE
if(obs==1) {
  pop=rbind(B.ini)
  for(t in 1:T) {
    Y=Y.obs[t]
    yields=c(yields,Y)
    pop.new=dynamics(pop[t,],Y)
    pop=rbind(pop,pop.new)
  }
  biomass=pop
  plot(times,biomass,type="b",ylim=c(0,1.1*k))
}
t.end=proc.time()
proceso=t.end-t.ini
proceso=proceso[3]
horas=floor(proceso/3600)
proceso=proceso-3600*horas
minutos=floor(proceso/60)
proceso=proceso-60*minutos
segundos=round(proceso,0)
print(paste("tiempo de proceso=",horas,"h",minutos,"m",segundos,"s",sep=""))
