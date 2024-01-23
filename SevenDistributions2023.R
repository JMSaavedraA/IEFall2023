#   ____________________________________________________________________________
#   LIKELIHOOD INFERENCES FOR THE GAMMA DISTRIBUTION                        ####
#   Exponential, Normal, Lognormal, Weibull, Gumbel and 
#   Inverse Gaussian distributions are also considered.
#   Note that fitting a Weibull to data is equivalent to fitting 
#   a Gumbel to logdata. But a Gumbel is also fitted here to the data
#   Last revised by Eloísa on October 2023

##  .........................................................................
##       Cleaning workspace area:                                ####
rm(list=ls(all=TRUE))

##  ............................................................................
##  Required functions:                                                     ####
#   -------------------------------- --
#   -------------------------------- --
#  The following functions will be used in the main program:

# This function is used to find the root of the first derivative of lp(alpha)=0:
Lgampsi=function(ALPHA,T1,T3,N){
  if (ALPHA>0)
   { sali=log(ALPHA*N/T1)-digamma(ALPHA) +(T3/N)
   } else  {
    sali=55555
   }
  return(sali)
}
#------------------------- --
# This function calculates the logprofile likelihood of the 
#  Gamma shape parameter:
lpgamalfa=function(alfas,T1,T3,N) {
 A=sum(alfas<=0)
    if (A>1) {'Warning alpha must be positive'
      lp= -14
    } else {
      lp = -N*lgamma(alfas) + (N*alfas*log(alfas)) -
            N*alfas*log(T1/N) + alfas*(T3-N)
    }
return(lp)
}

#------------------------- --
# This function calculates Gamma's rp(alpha)-lnc (for profile like. intervals)
rpgamalfalnc=function(alfas,T1,T3,N,ALPHAMLE,CC){
  rsali=lpgamalfa(alfas,T1,T3,N)-lpgamalfa(ALPHAMLE,T1,T3,N)-log(CC)
  return(rsali)
} 

#------------------------- --
# Gamma loglikelihood of shape parameter alpha and Gamma mean mu
lvgamalfamu=function(vec2,T1,T3,N){
   ALFA=vec2[1]
   MU=vec2[2]
   if (ALFA>0 & MU>0) {
   lv=-N*lgamma(ALFA)+ N*ALFA*log(ALFA)+(ALFA*T3)-
       (ALFA*T1/MU)-(N*ALFA*log(MU))
   } else {
     lv=-50;
   }
}
  
#------------------------- --
# This function is used for
# computing numerically profile likelihood of Gamma mean:
lpgamamu=function(vec1,T1,T3,N,MU){
  ALFA=vec1
    if (ALFA>0 & MU>0) {
    lv=-N*lgamma(ALFA)+ N*ALFA*log(ALFA)+(ALFA*T3)-
      (ALFA*T1/MU)-(N*ALFA*log(MU))
  } else {
    lv=-50;
  }
}

#------------------------- --
# This function is used for obtaining numerically
# profile like. intervals of Gamma mean mu by calculating rp(mu)-lnc:
rpgamamulnc=function(MU,T1,T3,N,ALPHAMLE,MUMLE,CC){
  #rp(MU) at MU is calculated first:
  interv=c(0.001,ALPHAMLE*15)  #search interval for alpha
  auxsali=optimize(lpgamamu,interval=interv,T1,T3,N,MU,maximum=TRUE,
                   lower=min(interv),upper=max(interv),tol=0.0000000000001)             
  lpgmu=auxsali$objective   
  sali=lpgmu-lpgamamu(ALPHAMLE,T1,T3,N,MUMLE)-log(CC)
  return(sali)
}

#--------------------------------- -
# Computing likelihood levels for usual confidence levels:
# THese are for estimating Gamma parameters:
cgamlevels=function(N){
  # Obtained jointly with Elías Mercado for his BA Thesis 
  #for a single Gamma parameter:
  cg90=0.2585-(1/(-0.5171+(1.6667*n)))
  cg95=0.1465-(1/(2.029+(2.084*n)))
  cg99=0.0362-(1/(14.609+(4.886*n)))
  cgamint=c(cg90,cg95,cg99)
  names(cgamint)=c("90%","95%","99%")
  #for both Gamma parameters:
  cgr90=0.1-(1/(1.060+(4.617*n)))
  cgr95=0.05-(1/(4.932+(7.030*n)))
  cgr99=0.01-(1/(51.82+(21.66*n)))
  cgamr=c(cgr90,cgr95,cgr99)
  names(cgamr)=c("90%","95%","99%")
  cesgam=list("cgamint"=cgamint,"cgamr"=cgamr)
  return(cesgam)
}

#These are likelihood levels for normal parameters, for mean mu, 
# scale sigma and for joint estimation (mu,sigma):
cnorlevels=function(n){
  #para mu:
  cc90=0.2585-(1/(-0.4743+(1.9029*n)))
  cc95=0.1465-(1/(0.8469+(2.3544*n)))
  cc99=0.0362-(1/(8.753+(5.551*n)))
  ccesmu=c(cc90,cc95,cc99)
  names(ccesmu)=c("90%","95%","99%")
  #para sigma:
  cs90=0.2585-(1/(0.775+(1.641*n)))
  cs95=0.1465-(1/(2.757+(2.082*n)))
  cs99=0.0362-(1/(18.294+(5.229*n)))
  ccessig=c(cs90,cs95,cs99)
  names(ccessig)=c("90%","95%","99%")
  #para (mu,sigma)
  ctous90=0.10-(1/(1.11+(4.655*n)))
  ctous95=0.05-(1/(1.61+(7.27*n)))
  ctous99=0.01-(1/(19.51+(23.59*n)))
  ctous=c(ctous90,ctous95,ctous99)
  names(ctous)=c("90%","95%","99%")
  cesn=list("cmu"=ccesmu,"csig"=ccessig,"cmusig"=ctous)
  return(cesn)
}



#------------------------- --
#The following 3 functions are used for inferences about the Weibull 
# Gumbel distributions:
K1bet=function(BETA,XX){
  if (BETA>0)
  { k1=sum(exp(BETA*XX))
  } else {
    cat("Error in K1: Beta must be positive")
    k1=0}
  return(k1)}
#  -------------------- --
K2bet=function(BETA,XX){
  if (BETA>0)
  { k2=sum(XX*exp(BETA*XX))
  } else {
    cat("Error in K2: Beta must be positive")
    k2=0}
  return(k2)}
#  --------------------- --
K3bet=function(BETA,XX){
  if (BETA>0)
  { k3=sum((XX**2)*exp(BETA*XX))
  } else {
    cat("Error in K3: Beta must be positive")
    k3=0}
  return(k3)}
#  ----------------------- -- 
#  This function is used to find the mle of Weibull shape parameter beta,
#   as the root of the first derivative of the log profile like. lp(beta)
#   Its arguments are beta, logdata XX, their sum TT, and sample size NN.
Weibetmle=function(BETA,XX,TT,NN){
  if (BETA>0){
    k1=K1bet(BETA,XX)
    k2=K2bet(BETA,XX)
    sali=(NN/BETA)+TT-(NN*k2/k1)
  }else {
    sali=5555}
  return(sali)}
#  ----------------------- -- 
#  This function calculates the Weibull log likelihood of (a,beta).
#  Beta is the Weibull shape parameter and a is the Weibull log scale.
#  Its arguments are vector VEC2=(a,beta), logdata XX, TT, and NN.
lvabetaWei=function(VEC2,XX,TT,NN){
  AA=VEC2[1]
  BETA=VEC2[2]
  if (BETA>0){
    k1=K1bet(BETA,XX)
    sali=NN*log(BETA)+(BETA*TT)-(NN*AA*BETA)-exp(-AA*BETA)*k1
  }else {
    sali=-999999999999}
  return(sali)
}

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
# Normal  rp(mu)                                                        
rpmu=function(mu,t1,t2,n){
  HH=t2-(t1^2)/n
  muhat=t1/n
  Hmu=t2-(2*t1*mu)+(n*(mu^2))
  rp=-(n/2)*(log(Hmu)-log(HH))
  return(rp)}  


#   ____________________________________________________________________________
#                 MAIN PROGRAM                                              ####
#  ............................. ........

#   Either simulate or input Gamma Data as Y here:

##  ............................................................................
#          Simulated Gamma Data                                           ####
n=15 #Sample size
alpha=2 #Gamma shape parameter
mugam=15 #Gamma mean
beta=mugam/alpha  #Gamma scale parameter
# Gamma data are simulated here:
X = rgamma(n,shape=alpha,scale=beta)
xsort = sort(X) #ordered data

##  ............................................................................
#   Gamma Data: Grice and Bain 20 mice subject to radiation                 ####
# This data set has been considered as Gamma and Weibull    
  X =  c(40,62,69,77,83,88,94,101,109,115,123,125,128,136,
         137,152,152,153,160,165)
  n=length(X)
  xsort = sort(X)
  
##  ............................................................................
#  Holes in Bread "Pan Kruasan Bimbo"                                     ####
#  Mean length of all holes per loaf in millimeters:
  X=c(9.0,11.0,30.3,16.5,42.3,26.8,21.0,1.0,21.0,38.7,12.5,20.5,6.0,14.0,19.5,
      20.5,34.6,36.7,29.0,16.5,31.2,12.7,19.5,26.2,10.3,12.7,16.8,15.9,14.6,
      10.8,23.4,22.3,18.8,17.3,14.6,10.0,13.5,17.6,25.5,10.8,9.8,17.6,16.2,
      14.7,16.9,17.3,13.7,11.8,15.5,13.3,13.5,11.7,10.0,19.4,15.9,17.5,19.0,
      12.3,13.6,8.8,14.3,14.6,14.5,18.4,15.0,16.3,12.6,15.1,15.3,16.5,14.0,
      21.2,20.1,11.5,14.3,17.0,17.0,9.5,20.5,20.6,12.0,18.6,24.0,13.5,18.0,
      11.6,14.4,10.3,17.0,17.0,29.0,24.4,12.6,19.5,17.5,13.3,10.3,26.4,11.5,
      17.1,20.0,19.9,11.3,13.7,19.8,13.6,16.5,13.5,15.9,14.8,15.0,12.5,26.0)
  xsort=sort(X)
  n=length(X)

  #---------- --
  # Datos Hirose
  X=c(2.1,3.4,2.7,3.2,2.5,2.7,2.7,3.4,3.3,3.9,2.9,3,3.4,3.5,3.5,3.4,3.4,3.6)
  xsort=sort(X)
  n=length(X)
  ##  ............................................................................
  # Trilobitas Counts                                                       ####    
  #    
  DATMAT=read.csv(file="C:/Users/Eloisa Diaz Frances/Google Drive/RprogsEloisa2022/GammaInferences2023/all_counts.csv",
                  header=TRUE)
  Contis=DATMAT[1:1835,2]
  X=log(Contis)+.5
  xsort=sort(X)
  n=length(X)
  
  
  ##  ............................................................................
  #  Maximum length of Holes per loaf                                    ####    
# Maximum length of hole per loaf of bread "Pan Croissant Bimbo" :   
DATMAT=read.csv(file="C:/Users/Eloisa Diaz Frances/Google Drive/RprogsEloisa2022/GammaInferences2023/Agujeros.csv",
                header=TRUE)
  Maxis=DATMAT[,3]
  X=Maxis
  xsort=sort(X)
  n=length(X)

# ____________________________________ ####
#        Statistics and MLES                                              ####
  t1=sum(X) # Sufficient statistic for Exp, Normal, Lognormal, Gamma and 
            # InvGaussian distributions
  t2=sum(X**2) # Sufficient statistic for Normal and Lognormal
  t3=sum(log(X)) # Second sufficient Statistic for Gamma distribution
  t4=sum(1/X)  #Second sufficient statistic for INverse Gaussian distrib

    
##  ............................................................................
#   Gamma MLEs                                                             ####
# Method of moments estimator of Gamma shape parameter:
alphamom = (t1**2)/((n*t2)-(t1**2))
#mle of Gamma, Exponential, and Normal mean:
mumle=t1/n
#mles of Gamma shape and scale (alpha,betagam) parameters:
alphamle=uniroot(Lgampsi,c(alphamom/10,alphamom*10),t1,t3,n,
                lower=alphamom/10,upper=alphamom*10,tol=0.000001)$root
betagamle=mumle/alphamle

##  ............................................................................
#   Weibull mles                                                 ####
# Log data:
Z=log(X)
tz=sum(Z)
zsort=sort(Z)
zbar=mean(Z)
#  The Gumbel method of moments estimate of Weibull shape beta is:
betaWeimom=pi*sqrt(n)/sqrt(6*sum((Z-zbar)**2))
# The mle of beta is calculated here:
betaWei=uniroot(Weibetmle,c(betaWeimom/10,betaWeimom*10),Z,tz,n,
                lower=betaWeimom/10,upper=betaWeimom*10,tol=0.000001)$root
#  Calculating the mle of "a", Gumbel location parameter:
K1=K1bet(betaWei,Z)
K2=K2bet(betaWei,Z)
K3=K3bet(betaWei,Z)
amle=log(K1/n)/betaWei
##  ............................................................................
#  Gumbel mles                                                 ####
# Original data turn to be the logdata:
ZG=X
tzG=sum(ZG)
zsortG=sort(ZG)
#  The Gumbel method of moments for beta for GUMBEL data X:
zbarG=mean(ZG)
betaGmom=pi*sqrt(n)/sqrt(6*sum((ZG-zbarG)**2))
# The mle of beta is calculated here:
betaGum=uniroot(Weibetmle,c(betaGmom/10,betaGmom*10),ZG,tzG,n,
                lower=betaGmom/10,upper=betaGmom*10,tol=0.000001)$root
#  Calculating the mle of a location Gumbel:
K1G=K1bet(betaGum,Z)
K2G=K2bet(betaGum,Z)
K3G=K3bet(betaGum,Z)
amleGum=log(K1G/n)/betaGum

##  ............................................................................
#  Inverse Gaussian mles                                                 ####
# The mean mle is also mumle
lambmle=1/((t4/n)-(n/t1))


# ____________________________________ ####
#                 MODEL VALIDATION                                          ####
##  ............................................................................
#  AICs:----------                                                          ####
##  ............................................................................
#   AIC Gamma                                                              ####
LnfgamY = -n*lgamma(alphamle) - (n*alphamle*log(betagamle)) +
         (alphamle-1)*t3 - (t1/betagamle);
AICgam=-2*LnfgamY + 4
# ---------------------- ---
#   AIC Exponential  ----
AICexp = (2*n*log(t1/n)) + 2*n+2
# ---------------------- ---
#   AIC Normal                                                              ####
H=t2 - ((t1^2)/n)
sigmle=sqrt(H/n)
AICnor=n*log(2*pi*H/n)+n+4
# ---------------------- ---
#   AIC LogNormal                                                           ####
Z=log(xsort)
t1LN=sum(Z)
t2LN=sum(Z^2)
HLN=t2LN-(t1LN^2)/n
# Lognormal mles:
muLN=t1LN/n
sigLN=sqrt(HLN/n)
AICLognor = 2*(t1LN) + n*log(2*pi*HLN/n) + n + 4 

# ---------------------- ---
#   AIC Weibull for Data X, Gumbel for LogData Y=lnX                             ####
# The estimated Weibull logdensity at (amle,betaWei) is:
# LnfWeiX=lvabetaWei(c(amle,betaWei),Z,tz,n) - tz
# AICWei= -2*LnfWeiX + 4

AICWei= 2*t3*(1-betaWei) -2*n*log(betaWei) +
        2*n*amle*betaWei + 2*exp(-amle*betaWei)*K1 +4
# ---------------------- ---
#   AIC Gumbel for Data X,     ####
# amleGum y betaGum calculados para datos X que son Gumbel
AICGum=-2*n*log(betaGum)-2*betaGum*tzG +2*n*amleGum*betaGum + 4 +
        2*exp(-amleGum*betaGum)*sum(exp(betaGum*ZG))

# ---------------------- ---
#   AIC Inverse Gaussian for Data X    ####
AICGauInv= n*log(2*pi) + log(lambmle) + 3*t3 -
         (2*n)/mumle + lambmle*((t1/mumle)+t4)+4


# ---------------------- ---
"Akaike Information Criteria for estimated distributions:"
rbind(AICexp,AICgam,AICnor,AICLognor,AICWei,AICGum,AICGauInv)

##  ............................................................................
#  AVD Average Vertical Distances in Probability Plots                     ####

# Transformed data for Exp, Gamma, Normal, Lognormal, and Weibull distributions:
uexp=1-exp(-xsort/mumle)
ugam=pgamma(xsort,shape=alphamle,scale=betagamle)
unor=pnorm(xsort,mean=mumle,sd=sigmle)
ulognor=plnorm(xsort,meanlog=muLN,sdlog=sigLN) 
#Weibull distrib for data X is equivalent to Gumbel for logdata Z
Zsort=sort(Z)
ugum=1-exp(-exp(betaWei*(Zsort-amle)))
#GUmbel distrib directly for data X=ZG:
ZGsort=sort(ZG)
ugumZG=1-exp(-exp(betaGum*(ZGsort-amleGum)))
#Inverse Gaussian:
aux1=sqrt(lambmle/xsort)*((xsort/mumle)-1)
aux2=-sqrt(lambmle/xsort)*((xsort/mumle)+1)
uGauInv=pnorm(aux1,mean=0,sd=1)+exp(2*lambmle/mumle)*pnorm(aux2,mean=0,sd=1)

# Ebetas are the theoretical quantiles of Uniform(0,1) distribution
# of same probabilities:
enes=(1:n)
Ebetas=enes/(n+1)
varian=(enes*(n+1-enes))/(((n+1)^2)*(n+2))
## Joint 95% confidence band for N prediction intervals of uniform quantiles:
tau=1-(.05/n)
tau1=(1-tau)/2
tau2=(1+tau)/2
aenes= n+1-enes
Ban1=qbeta(tau1,enes,aenes)
Ban2=qbeta(tau2,enes,aenes)
# AVD (average Vertical distances of points to identity line are calculated:
AVDExp=sum(abs(uexp-Ebetas)/sqrt(varian))/n
AVDGam=sum(abs(ugam-Ebetas)/sqrt(varian))/n
AVDNor=sum(abs(unor-Ebetas)/sqrt(varian))/n
AVDLognor=sum(abs(ulognor-Ebetas)/sqrt(varian))/n
AVDWei=sum(abs(ugum-Ebetas)/sqrt(varian))/n

"AVD for Exp, Gamma, Normal, Lognormal, and Weibull distributions:"
rbind(AVDExp,AVDGam,AVDNor,AVDLognor,AVDWei)


# Proposed probability plots for 7 distributions:
X11()
plot(Ebetas,uexp,pch=19,cex=.5,ylab="F(xi;mumv)",
     xlab="Cuantiles uniformes (0,1)", xlim=c(0,1),ylim=c(0,1),
     main="Gráfica de probabilidad Exponencial",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)

X11()
plot(Ebetas,ugam,pch=19,cex=.5,ylab="F(xi;mumv,alfamv)",
     xlab="Cuantiles uniformes (0,1)", xlim=c(0,1),ylim=c(0,1),
     main="Gráfica de probabilidad Gamma",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)

X11()
plot(Ebetas,unor,pch=19,cex=.5,ylab="F(xi;mumv,sigmamv)",
     xlab="Cuantiles uniformes (0,1)", xlim=c(0,1),ylim=c(0,1),
     main="Gráfica de probabilidad Normal",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)

X11()
plot(Ebetas,ulognor,pch=19,cex=.5,ylab="F(xi;muLNmv,sigmaLNmv)",
     xlab="Cuantiles uniformes (0,1)", xlim=c(0,1),ylim=c(0,1),
     main="Gráfica de probabilidad Lognormal",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)

X11()
plot(Ebetas,ugum,pch=19,cex=.5,ylab="F(xi;betaWei,aGum)",
     xlab="Cuantiles uniformes (0,1)", xlim=c(0,1),ylim=c(0,1),
     main="Gráfica Prob. Weibull",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)



X11()
plot(Ebetas,uGauInv,pch=19,cex=.5,ylab="F(xi;betaWei,aGum)",
     xlab="Cuantiles uniformes (0,1)", xlim=c(0,1),ylim=c(0,1),
     main="Gráfica Prob. Gaussiana Inversa",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)


X11()
plot(Ebetas,ugumZG,pch=19,cex=.5,ylab="F(xi;betaWei,aGum)",
     xlab="Cuantiles uniformes (0,1)", xlim=c(0,1),ylim=c(0,1),
     main="Gráfica Prob. Gumbel para datos X",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)



#_________________________________________ ####



t1 <- sum(X)
t2 <- sum(X^2)
H <- t2 - t1^2/113
qchisq(0.025,112)
qchisq(0.975,112)

pt(2.1584)