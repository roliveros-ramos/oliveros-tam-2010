##############################################
# MANAGMENT STRATEGY EVALUATION
# Probability map
# Ricardo Oliveros-Ramos @ CIMOBP - IMARPE.
# e-mail: roliveros@imarpe.gob.pe
# Versi?n 1.1 - 24/06/2009
##############################################
##################################
#### Simulation control ##########
##################################
t.ini=proc.time()
# strategy="E2" 	# to be chosen beetween E1 (min), E2 (precautionary ICES)
# E3 (max viable), E4 (precautionary viable), E5 (viable)
obs=0 				# 1 for observed catches
N=100 				# number of simulations (Montecarlo)
T=50 # time horizon
spinup=20
colapse=TRUE 		# Scenarium
sce=if(colapse) "post" else "pre"
uncertainty=TRUE	# model uncertainty on/off
fpa=seq(from=0.2,to=1.2,by=0.02) # 0.02 
bpa=seq(form=1,to=4,by=0.05) # 0.05
estra="E6"
est=rep(estra,length(fpa)*length(bpa))
#est="mor" # moratorium
SSBs=c()
refpoint.1=rep(bpa,times=rep(length(fpa),length(bpa))) # Bpa
refpoint.2=rep(fpa,length(bpa)) # Fpa
refpoints=cbind(refpoint.1,refpoint.2)
n=1; t=1
##################################
#### Parameters ##################
##################################
# population growth
R=1.8
k=1 # carrying capacity
K=R*k/(R-1)
# Reference points
Bmin=0.2
Bpa=1.5*Bmin
Ymin=0.05 # variar
Ymin.viab=0.068 # verificar, MVY
Fpa=0.30
Fviab=0.20 #no se usa
# uncertainty
sigma.mu=0.1
sigma.sd=sqrt(0.09)
omega.mu=0.05
omega.sd=sqrt(0.04)
sigma=rlnorm(N*T,mean=log(sigma.mu),sd=sigma.sd)
omega=rlnorm(N*T,mean=log(omega.mu),sd=omega.sd)
factor=if(uncertainty) (1+sigma)*(1+omega) else double(N*T)+1
factor=matrix(factor,ncol=T,nrow=N)
##################################
#### Initial conditions ##########
##################################
# Y.health=as.matrix(read.table("yield.txt"))
# Y.colapse=as.matrix(read.table("yield2.txt"))
# Y.obs=if(colapse) Y.colapse else Y.health
# Initial biomasses
B.health=2*Bmin # healthy stock
B.colapse=0.9*Bmin # post-colapse stock
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
dynamics=function(B,Y) {
  stock=B-Y
  stock=R*stock*(1-stock/K)
  stock=max(stock,0)
  return(stock)
}
##################################
#### Viability results ###########
##################################
Ymax= function(B) B-(0.5*R*k/(R-1))*(1-sqrt(1-(4*(R-1)*Bmin)/(R*R*k)))
test=function(B) {
  TEST=if(B>=Bmin) TRUE else FALSE
  return(TEST)
}
MVY=Ymax(Bmin)
##################################
#### Management strategies #######
##################################
yield=function(B,strategy) {
  Y=if(B<Bmin) 0 else switch(strategy,
                             E1=Ymin,
                             E2=Ymax(B),
                             E3=MVY,
                             E4=0.5*(Ymax(B)+Ymin),
                             E5=runif(1,Ymin,Ymax(B)),
                             E6=if(B>=Bpa) Fpa*B else Ymin,
                             E7=if(B>=Bpa) 0.5*(Fpa*B+Ymin) else Ymin,
                             E8=if(B>=Bpa) runif(1,Ymin,Fpa*B) else Ymin,
                             mor=0)
  return(Y)
}
##################################
#### Analysis ####################
##################################
ref.1=function(x) {
  lim1=c(mean(x))
  return(lim1)
}
##################################
#### Simulation ##################
##################################
if(obs==0) {
  B.est=cbind()
  Prob.lt=rbind()
  Prob.min=rbind()
  Prob.pa=rbind()
  B.ind=rbind()
  Y.ind=rbind()
  B.serie=array(0,dim=c(length(times),length(ref.1(1)),length(est)))
  Y.serie=array(0,dim=c(length(times)-1,length(ref.1(1)),length(est)))
  i=0
  for(strategy in est) {
    i=i+1
    Bpa=refpoints[i,1]*Bmin
    Fpa=refpoints[i,2]
    Fviab=refpoints[i,2]
    B.t=cbind()
    Y.t=cbind()
    prob.lt=0
    prob.min=0
    prob.pa=0
    for(n in 1:N) {
      pop=rbind(B.ini)
      yields=c()
      for(t in 1:T) {
        Y=yield(pop[t,],strategy)
        pop.new=dynamics(pop[t,],Y*factor[n,t])
        pop=rbind(pop,pop.new)
        yields=c(yields,Y)
      }
      B.t=cbind(B.t,pop)
      Y.t=cbind(Y.t,yields)
      prob.lt=prob.lt+min(prod(floor(pop[1:spinup]/Bmin)),1)
      prob.min=prob.min+min(prod(floor(pop[-(1:spinup)]/Bmin)),1)
      prob.pa=prob.pa+min(prod(floor(pop[-(1:spinup)]/Bpa)),1)
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
    B.ind.1=mean(apply(B.t[-(1:spinup),],2,mean))/Bmin
    B.ind.2=mean(apply(B.t[-(1:spinup),],2,sd))/mean(apply(B.t[-(1:spinup),],2,mean))
    B.ind=rbind(B.ind,c(B.ind.1,B.ind.2))
    Y.ind.1=mean(apply(Y.t[-(1:spinup),],2,mean))/Ymin
    Y.ind.2=mean(apply(Y.t[-(1:spinup),],2,sd))/mean(apply(Y.t[-(1:spinup),],2,mean))
    Y.ind=rbind(Y.ind,c(Y.ind.1,Y.ind.2))
  }
  ind=cbind(refpoints,Prob.lt,Prob.min,Prob.pa,B.ind,Y.ind)
  vars=c("B","Y")
  inds=c("Bpa","Fpa","Prob.lt","Prob.min","Prob.pa")
  for(v in vars) {
    indss=paste(v,".ind.",1:2,sep="")
    inds=c(inds,indss)
  }
  dimnames(ind)=list(est,inds)
  write.table(ind,file=paste(sce,"_indicadores-BF-N",N,"-",estra,".csv",sep=""), quote = FALSE, sep = ",", col.names = TRUE, row.names = TRUE)
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
