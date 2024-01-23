# Inferencias para datos con censura para 7 distribuciones estadísticas
# Revisado por Eloísa Díaz Francés en octubre de 2023

##  .........................................................................
##   Limpieza espacio de trabajo:                                ####
rm(list=ls(all=TRUE))

# Instalar antes elpaquete "actuar" para poder usar comandos de la Gaussiana Inv.
library(actuar)
##  ............................................................................
##  Funciones que se usarán:                                                     ####
#   -------------------------------- --
# Estadística Anderson Darling para comparación de modelos. 
# Se le alimenta el vector ZZ de datos ordenados transformados con la
# distribución estimada: ZZ=Fhat(xsort).
AndersonDarling=function(ZZ){
  nn=length(ZZ)
  ns=seq(1,nn,by=1)
  aux1=(2*ns-1)
  aux2=(2*nn+1-(2*ns))
  AD= - nn - (1/nn)*sum((aux1*log(ZZ))+aux2*log(1-ZZ))
  return(AD)
}
#--------------- --
# Estadística propuesta por Eloísa del promedio ponderado de distancias
# verticales en la gráfica de probabilidad.
WAVD=function(ZZ){
  nn=length(ZZ)
  ns=seq(1,nn,by=1)
  aens= nn+1-ns
  Ebetas=ns/(n+1)
  varian=(ns*aens)/(((nn+1)^2)*(nn+2))  
  sali = sum(abs(ZZ-Ebetas)/sqrt(varian))/nn
  return(sali)
}
#--------------- --
# Función para hallar la raíz de primera derivada de lp(alfa)=0 para el
# parámetro de forma alfa de la distribución Gamma.
Lgampsi=function(ALPHA,T1,T3,N){
  if (ALPHA>0)
  { sali=log(ALPHA*N/T1)-digamma(ALPHA) +(T3/N)
  } else  {
    sali=55555
  }
  return(sali)
}
#------------------------- --
# Log verosimilitud perfil del parámetro forma Gamma:
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
# Para encontrar numericamente intervalos verosim perfil de nivel c 
# para parámetro alfa de la Gamma se calcula rp(alpha)-lnc 
rpgamalfalnc=function(alfas,T1,T3,N,ALPHAMLE,CC){
  rsali=lpgamalfa(alfas,T1,T3,N)-lpgamalfa(ALPHAMLE,T1,T3,N)-log(CC)
  return(rsali)
} 

#------------------------- --
# Logverosimilitud (alpha,mu) para la distribución Gamma 
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
# Para calcular numéricamente la logperfil de la media Gamma:
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
# Para obtener intervalos perfiles de la media GAmma, se calcula
#  rp(mu)-lnc:
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
# Se calculan niveles de verosimilitud de intervalos y regiones para 
# inferencias con la distribución Gamma (uno y dos parámetros)
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

# Se calculan niveles de verosimilitud de intervalos y regiones para 
# inferencias con la distribución Normal (uno y dos parámetros)
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
# Para inferencias con distribuciones Weibull y Gumbel de mínimos: 
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
#  Para encontrar el emv del parámetro de forma beta Weibull dando con la 
# raíz de la primera derivada de la perfil de beta igual a cero:
# Los argumentos son beta, los logdatos XX, s suma TT y el tamaño muestral NN.
Weibetmle=function(BETA,XX,TT,NN){
  if (BETA>0){
    k1=K1bet(BETA,XX)
    k2=K2bet(BETA,XX)
    sali=(NN/BETA)+TT-(NN*k2/k1)
  }else {
    sali=5555}
  return(sali)}
#  ----------------------- -- 
#  Se calcula logverosimilitud Weibull (a,beta), beta es forma y aa logescala.
#  Los argumentos son VEC2=(a,beta), logdata XX, sumaTT, y tamaño mtra NN.
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
#  Funciones para datos censurados por intervalo ####
# Normal loglikelihood for interval censored data in matrix DATMAT (nx2)
# (to be maximized wiht OPTIM)
lvnorCEN=function(vec2,DATMAT){
  epsi=0.000001
  mu=vec2[1]
  sigma=vec2[2]
  if (sigma >0){
    # Prob. of each obs. of falling in its corresponding interval:
    difi = pnorm(DATMAT[,2],mu,sigma)-pnorm(DATMAT[,1],mu,sigma)
    # very small values are replaced before taking logarithm:
    difnew = (difi<epsi)*epsi + (difi>=epsi)*difi
    lv = sum(log(difnew))
  }else{
    lv=-99999999999999999 # discouraging impossible values for parameters
  }
  return(lv)
} 
#---------------------- --
# Lognormal loglikelihood for interval censored data in matrix DATMAT (nx2)
lvLognorCEN=function(vec2,DATMAT){
  epsi=0.000001
  muLN=vec2[1]
  sigmaLN=vec2[2]
  if (sigmaLN >0){
    # Prob. of each obs. of falling in its corresponding interval:
    difi = plnorm(DATMAT[,2],meanlog=muLN,sdlog=sigmaLN) - 
      plnorm(DATMAT[,1],meanlog=muLN,sdlog=sigmaLN)
    # very small values are replaced before taking logarithm:
    difnew = (difi<epsi)*epsi + (difi>=epsi)*difi
    lv = sum(log(difnew))
  }else{
    lv=-99999999999999999 # discouraging impossible values for parameters
  }
  return(lv)
} 
#----------------------- --
# Gamma loglikelihood for interval censored data in matrix DATMAT (nx2)
lvGamCEN=function(vec2,DATMAT){
  epsi=0.000001
  alfa=vec2[1]
  mugam=vec2[2]
  betagam=mugam/alfa
  if ((alfa>0) & (mugam>0)) {
    # Prob. of each obs. of falling in its corresponding interval:
    difi = pgamma(DATMAT[,2],shape=alfa,scale=betagam) - 
      pgamma(DATMAT[,1],shape=alfa,scale=betagam)
    # very small values are replaced before taking logarithm:
    difnew = (difi<epsi)*epsi + (difi>=epsi)*difi
    lv = sum(log(difnew))
  }else{
    lv=-99999999999999999 # discouraging impossible values for parameters
  }
  return(lv)
} 
#---------------------- --
# Exponential loglikelihood for interval censored data in matrix DATMAT (nx2)
# This distribution is parametrized in terms of its mean theta
lvExpCEN=function(vec1,DATMAT){
  epsi=0.000001
  theta=vec1
  if (theta>0) {
    # Prob. of each obs. of falling in its corresponding interval:
    difi = pexp(DATMAT[,2],rate=1/theta) - pexp(DATMAT[,1],rate=1/theta)
    # very small values are replaced before taking logarithm:
    difnew = (difi<epsi)*epsi + (difi>=epsi)*difi
    lv = sum(log(difnew))
  }else{
    lv=-99999999999999999 # discouraging impossible values for parameters
  }
  return(lv)
} 
#---------------------- --
lvWeiCEN=function(vec2,DATMAT){
  epsi=0.000001
  betaWei=vec2[1] #shape parameter
  sigWei=vec2[2]  #scale parameter
  if ((betaWei>0) & (sigWei>0)) {
    # Prob. of each obs. of falling in its corresponding interval:
    difi = pweibull(DATMAT[,2],shape=betaWei,scale=sigWei) - 
      pweibull(DATMAT[,1],shape=betaWei,scale=sigWei)
    # very small values are replaced before taking logarithm:
    difnew = (difi<epsi)*epsi + (difi>=epsi)*difi
    lv = sum(log(difnew))
  }else{
    lv=-99999999999999999 # discouraging impossible values for parameters
  }
  return(lv)
} 

#---------------------- --
#Gumbel distribution of minima:
pGumbelmin=function(xx,vec2){
  agum=vec2[1] #Gumbel location
  bgum=vec2[2] #Gumbel scale
  if (bgum>0){
    Fgum=1-exp(-exp((xx-agum)/bgum))
  } else { 
    Fgum=999
  }
  return(Fgum)  
}
#Int censored loglikelihood for Gumbel of minima:
lvGuminCEN=function(vec2,DATMAT){
  epsi=0.000001
  agum=vec2[1] #location parameter
  bgum=vec2[2]  #scale parameter
  if (bgum>0) {
    # Prob. of each obs. of falling in its corresponding interval:
    difi = pGumbelmin(DATMAT[,2],c(agum,bgum)) - 
      pGumbelmin(DATMAT[,1],c(agum,bgum))
    # very small values are replaced before taking logarithm:
    difnew = (difi<epsi)*epsi + (difi>=epsi)*difi
    lv = sum(log(difnew))
  }else{
    lv=-99999999999999999 # discouraging impossible values for parameters
  }
  return(lv)
} 

#---------------------- --
#--------Inverse Gaussian log likelihood with interval censoring: --- --
lvGauInvCEN=function(vec2,DATMAT){
  epsi=0.000001
  muGI=vec2[1] #mean parameter
  lamGI=vec2[2]  #shape parameter
  if ((lamGI>0) & (muGI>0)) {
    # Prob. of each obs. of falling in its corresponding interval:
    difi = pinvgauss(DATMAT[,2],mean=muGI,shape=lamGI) - 
      pinvgauss(DATMAT[,1],mean=muGI,shape=lamGI)
    # very small values are replaced before taking logarithm:
    difnew = (difi<epsi)*epsi + (difi>=epsi)*difi
    lv = sum(log(difnew))
  }else{
    lv=-99999999999999999 # discouraging impossible values for parameters
  }
  return(lv)
} 
# ____________________________________ ####
#               PROGRAMA PRINCIPAL                                    ####
#  ............................. ........
# En todos los casos dar un vector X con datos redondeados, aunque
# sea una aproximación
#---------------------------- --
# __________________________________________________________________ ####
#     A) Simulating Interval Censored Data                               ####
# Gamma rounded data are simulated, assuming that the precision of the 
# measuring instrument is W=2h
# The ratio of W to the median, Ratio=W/Q50gam, is relevant for the number
# of repeated observations of continuous random variables. When RR>0.2 there 
# are usually several repetitions of some observations.
n=60 #Sample size
Alphagam=1.72 #Gamma shape parameter
Mugam=1   5 #Gamma mean
Betagam=Mugam/Alphagam  #Gamma scale parameter
#Median of this Gamma distribution:
Q50gam=qgamma(0.50,shape=Alphagam,scale=Betagam)
#Either the Ratio or W must be specified. Here the Ratio is specified as:
Ratio=0.5;
# Precision of the measuring instrument is calculated as:
W=(Ratio*Q50gam)
h=W/2
# Exact Gamma data are simulated here:
Xexa = rgamma(n,shape=Alphagam,scale=Betagam)
xsort = sort(Xexa) #ordered data
xmax=max(Xexa)
cotamax=floor(xmax/W)*W + W
#Observed data will be rounded to one of these values:
# {0,2h,4h,...,(2k)h,...} (an even multiple of h)
VecH=seq(from=0,to=cotamax,by=W)
#Xround is the vector of registered data and may contain zeros
Xround=rep(0,n)
for(j in 1:n){
  indi=which.min(abs(xsort[j]-VecH))
  Xround[j]=VecH[indi]
}
#For calculating t1,t2,t3 zeros are replaced with h:
X=(Xround==0)*h +(Xround>0)*Xround #Vector without zeros
# Matrix of intervals is provided here:
XMAT=cbind(Xround-h,Xround+h)
#Any interval[-h,h] is changed to [0.000001,h]:
epsi=0.000001
XMAT[,1]=(XMAT[,1]==(-h))*(XMAT[,1]+h+epsi) + (XMAT[,1]>0)*XMAT[,1]
Xlogmat=log(XMAT)
#Different values are counted to acknowledge the ammount of repeated values:
Z=unique(X)
nuni=length(Z)
"For n interval cens. obs:"
n
"There were only unrepeated values:"
nuni
"The maximum number of repeated values was:"
max(table(X))
#  ................
##  ............................................................................
#   Gamma Data: Grice and Bain 20 mice subject to radiation                 ####
# This data set has been considered as Gamma and Weibull    
  X =  c(40,62,69,77,83,88,94,101,109,115,123,125,128,136,
         137,152,152,153,160,165)
  n=length(X)
  xsort = sort(X)
  epsi=1/2
  XMAT=cbind((xsort-epsi),(xsort+epsi))
  
  #---------- --
  # Datos Hirose
  X=c(2.1,3.4,2.7,3.2,2.5,2.7,2.7,3.4,3.3,3.9,2.9,3,3.4,3.5,3.5,3.4,3.4,3.6)
  xsort=sort(X)
  n=length(X)
  epsi=0.1/2
  XMAT=cbind((xsort-epsi),(xsort+epsi))
##  ............................................................................
  #  Longitudes máximas por hogaza de pan Cruapán de Bimbo                                   ####    
X=c(1,6,12,12,14,14,15,15,15,15,15,15,15,16,16,16,16,16,17,17,18,18,18,18,18,18,
    18,18,18,18,19,19,19,19,19,21,21,21,21,21,22,22,22,22,23,23,23,23,23,24,24,
    24,24,24,24,25,25,25,25,25,25,25,25,25,27,27,27,27,27,27,28,29,29,29,29,29,
    29,30,30,31,32,32,33,33,34,34,34,34,35,35,36,37,38,38,39,39,40,40,42,43,45,
    45,48,48,50,51,54,54,59,71,73,77,78)
  xsort=sort(X)
  n=length(X)
    Xuni=sort(unique(X))
  nuni=length(Xuni)
  difis=Xuni[2:nuni]-Xuni[1:(nuni-1)]
  difimin=min(difis)
  epsi=difimin/2
  XMAT=cbind((xsort-epsi),(xsort+epsi))
  
#----------------------------------------------------- --  
  
  # ____________________________________ ####
###  Emvs SIN CENSURA para valores iniciales   ####
  #    Estadísticas suficientes datos exactos:                     ####
  t1=sum(X) # Para distribuciones Exp, Normal, Lognormal, Gamma e Inv Gaussiana
  t2=sum(X**2) # Para la Normal y Lognormal
  t3=sum(log(X)) # Para la Gamma
  t4=sum(1/X)  # Para la Inversa Gaussiana

##  ............................................................................
  # EMV media Exp,Gamma y Normal: ####
  mumle=t1/n
  # EMVs Gamma                                                           ####
# Estimador momentos forma Gamma:
alphamom = (t1**2)/((n*t2)-(t1**2))

# emv de parámetros de forma Gamma, alphamle,  y escala, betagamle:
alphamle=uniroot(Lgampsi,c(alphamom/10,alphamom*10),t1,t3,n,
                lower=alphamom/10,upper=alphamom*10,tol=0.000001)$root
betagamle=mumle/alphamle
#------------------------ --
# EMVs Normales            ####
H=t2 - ((t1^2)/n)
sigmle=sqrt(H/n)

# EMVs Log Normales            ####
muLN=tz/n
tz2=sum(Z^2)
sigLN=sqrt((tz2/n)-(muLN^2))

##  ............................................................................
# EMVs Gaussiana Inversa:                                       ####
# La media mumle es el emv de mu. Este es el emv del parámetro de forma lambda:
lambmle=1/((t4/n)-(n/t1))

##  ............................................................................
# EMVs  Weibull                                             ####
# Log data:
Z=log(X)
tz=sum(Z)
zsort=sort(Z)
zbar=mean(Z)
#  Estimador de momentos Gumbel para forma beta de la Weibull:
betaWeimom=pi*sqrt(n)/sqrt(6*sum((Z-zbar)**2))
# EMV de forma beta aquí:
betaWei=uniroot(Weibetmle,c(betaWeimom/10,betaWeimom*10),Z,tz,n,
                lower=betaWeimom/10,upper=betaWeimom*10,tol=0.000001)$root
# Se calcula el emv de a, localización de los logdatos Gumbel.
# Notar que la escala WEibull sigma=exp(a).
K1=K1bet(betaWei,Z)
K2=K2bet(betaWei,Z)
K3=K3bet(betaWei,Z)
amle=log(K1/n)/betaWei
##  ............................................................................
#  Emvs Gumbel para datos X (escala original)                          ####
# Ahora Z=exp(X) son Weibull de parámetros amlegum,betaGum.
ZG=xsort
tzG=sum(ZG)
zsortG=ZG
# Estimador de momentos Gumbel para beta WEibull de Z=exp(X):
zbarG=mean(ZG)
betaGmom=pi*sqrt(n)/sqrt(6*sum((ZG-zbarG)**2))
# Emv de beta:
betaGum=uniroot(Weibetmle,c(betaGmom/10,betaGmom*10),ZG,tzG,n,
                lower=betaGmom/10,upper=betaGmom*10,tol=0.000001)$root
# Emv de localización a Gumbel:
K1G=K1bet(betaGum,ZG)
K2G=K2bet(betaGum,ZG)
K3G=K3bet(betaGum,ZG)
amleGum=log(K1G/n)/betaGum
tetamle=amleGum
deltamle=1/betaGum
#_________________________________________ ####
#  AICs SIN CENSURA para comparar           ####
##  ............................................................................
# ---------------------- ---
#   AIC Exponencial  
AICexp = (2*n*log(t1/n)) + 2*n+2
# ---------------------- ---
#   AIC Gamma                                       
LnfgamY = -n*lgamma(alphamle) - (n*alphamle*log(betagamle)) +
  (alphamle-1)*t3 - (t1/betagamle);
AICgam=-2*LnfgamY + 4
# ---------------------- ---
#   AIC Normal        
AICnor=n*log(2*pi*H/n)+n+4
# ---------------------- ---
#   AIC LogNormal     
Z=log(xsort)
t1LN=sum(Z)
t2LN=sum(Z^2)
HLN=t2LN-(t1LN^2)/n
# Lognormal mles:
muLN=t1LN/n
sigLN=sqrt(HLN/n)
AICLognor = 2*(t1LN) + n*log(2*pi*HLN/n) + n + 4 

# ---------------------- ---
#   AIC Weibull para datos X o Gumbel para logdatos Y=lnX 
# La logdensidad estiamda Weibull en (amle,betaWei) es:
# LnfWeiX=lvabetaWei(c(amle,betaWei),Z,tz,n) - tz
# AICWei2= -2*LnfWeiX + 4
# se calcula directo aquí el AIC Weibull:
AICWei= 2*t3*(1-betaWei) -2*n*log(betaWei) +
  2*n*amle*betaWei + 2*exp(-amle*betaWei)*K1 +4
# ---------------------- ---
#   AIC Gumbel para datos X:     
# amleGum y betaGum calculados para datos X que son Gumbel
AICGum=2*n*log(deltamle)-2*(t1-(n*tetamle))/deltamle  +
  2*exp(-tetamle/deltamle)*sum(exp(xsort/deltamle)) +4 
# ---------------------- ---
#   AIC Gaussiana Inversa para datos X    
AICGauInv= n*log(2*pi) - n*log(lambmle) + 3*t3 -
  (2*n*lambmle)/mumle + lambmle*((t1/(mumle^2))+t4) + 4
# ---------------------- ---
"Criterio de Info.de Akaike para datos exactos y distribuciones estimadas:"
AICSIN=c(AICexp,AICgam,AICnor,AICLognor,AICGauInv,AICWei,AICGum)
names(AICSIN)=c("Exp","Gam","Nor","Lnor","GauInv","Wei","Gum")
cbind(sort(AICSIN))


##_______________________________________________________ ####

### EMVS CON CENSURA:    ####
# Exponencial CEN    ####
SALIS=optimize(lvExpCEN,c(mumle/10,mumle*10),XMAT,
               lower=mumle/10,upper=mumle*10,tol=0.000001,
               maximum=TRUE)
muExpC=SALIS$maximum
#-------------------------- --
# Gamma CEN    ####
vecini=c(alphamom,mumle) #initial values of Gamma shape and mean 

maxversali=optim(vecini,lvGamCEN,gr="Null",XMAT,method="Nelder-Mead",
                 control=list(fnscale=-2) )
# if fnscale is negative, optim MAXIMIZES
#mles:
emvs=maxversali$par
alfaGamC=emvs[1]
mugamC=emvs[2]
betagamC=mugamC/alfaGamC
#loglikelihood at the mle:
LVmaxGamC=maxversali$value

#----------------- --
# Normal CEN    ####
vecini=c(mumle,sigmle) #initial values for mu and sigma
maxversali=optim(vecini,lvnorCEN,gr="Null",XMAT,method="Nelder-Mead",
                 control=list(fnscale=-2) )
# if fnscale is negative, optim MAXIMIZES
#mles:
emvs=maxversali$par
munorC=emvs[1]
signorC=emvs[2]
#Loglikelihood at the mle:
LVmaxnorC=maxversali$value
#-------------------- --
# Lognormal CEN    ####
Xlogmat=log(XMAT)
vecini=c(muLN,sigLN) #initial values are exact mles
maxversali=optim(vecini,lvnorCEN,gr="Null",Xlogmat,method="Nelder-Mead",
                 control=list(fnscale=-2) )
# if fnscale is negative, optim MAXIMIZES instead of minimizing the function
emvs=maxversali$par
muLNC=emvs[1]
sigLNC=emvs[2]
#loglikelihood at the mle:
LVmaxLNC=maxversali$value
#----------------- --
# Gaussiana Inversa CEN    ####
#convenient initial values InvGaussian parameters:
lambmle=1/((t4/n)-(n/t1))
vecini=c(mumle,lambmle)  
maxversali=optim(vecini,lvGauInvCEN,gr="Null",XMAT,method="Nelder-Mead",
                 control=list(fnscale=-2) )
# if fnscale is negative, optim MAXIMIZES
#mles:
emvs=maxversali$par
muGIC=emvs[1]
lambGIC=emvs[2]
#loglikelihood at the mle:
LVmaxGauInvC=maxversali$value

#----------------- --
# Weibull CEN    ####
#convenient initial values for Weibull parameters (beta shape and sigma scale):
betaWeini=pi*sqrt(n)/sqrt(6*sum((log(X)-mean(log(X)))**2))
aini=log(K1bet(betaWeini,log(X))/n)/betaWeini
sigWeini=exp(aini)
vecini=c(betaWeini,sigWeini) #initial values 
maxversali=optim(vecini,lvWeiCEN,gr="Null",XMAT,method="Nelder-Mead",
                 control=list(fnscale=-2) )
# if fnscale is negative, optim MAXIMIZES
#mles:
emvs=maxversali$par
betaWeiC=emvs[1]
sigWeiC=emvs[2]
#loglikelihood at the mle:
LVmaxWeiC=maxversali$value

#----------------- --
# Gumbel CEN    ####
#initial values are exact mles:
vecini=c(amleGum,deltamle) 
maxversali=optim(vecini,lvGuminCEN,gr="Null",XMAT,method="Nelder-Mead",
                 control=list(fnscale=-2) )
# if fnscale is negative, optim MAXIMIZES
#mles:
emvs=maxversali$par
amleGumC=emvs[1]
bGumC=emvs[2]
#loglikelihood at the mle:
LVmaxGuminC=maxversali$value
#----------------- --

#-------------------------------------- --

# Van EMVS sin y con CENSURA: ####
"Normal mles exact and cens for mu and sigma are:"
c(mumle,sigmle)
c(munorC,signorC)
"Logmormal mles exact and cens for mu and sigma are:"
c(muLN,sigLN)
c(muLNC,sigLNC)
"Exponential mean mles exact and cens:"
c(mumle,muExpC)
"Gamma shape,mean mles exact and cens:"
c(alphamle,mumle)
c(alfaGamC,mugamC)
"Gamma scale beta mle exact and cens:"
c(mumle/alphamle,mugamC/alfaGamC)
"Weibull shapebeta, scalesigma mles exact and cens:"
c(betaWei,exp(amle))
c(betaWeiC,sigWeiC)
"Inv Gaussian mean,shape mles exact and cens:"
c(mumle,lambmle)
c(muGIC,lambGIC)
##_______________________________________________________ ####
##  ......................................................  ..................
#   AIC CON CENSURA                                            ####
# Normal:
PrSamNorC=lvnorCEN(c(munorC,signorC),XMAT)
AICenNorC=-2*(PrSamNorC)+4
# Lognormal:
PrSamLN=lvnorCEN(c(muLNC,sigLNC),Xlogmat)
AICenLNC=-2*(PrSamLN)+4
# Gamma:
PrSamGamC=lvGamCEN(c(alfaGamC,mugamC),XMAT)
AICenGamC=-2*(PrSamGamC)+4
# Exponential:
PrSamExpC=lvExpCEN(muExpC,XMAT)
AICenExpC=-2*(PrSamExpC)+2
# Weibull:
PrSamWeiC=lvWeiCEN(c(betaWeiC,sigWeiC),XMAT)
AICenWeiC=-2*(PrSamWeiC)+4
# Gumbel min:
PrSamGumC=lvGuminCEN(c(amleGumC,bGumC ),XMAT)
AICenGumC=-2*(PrSamGumC)+4
# Inverse Gaussian:
PrSamIGC=lvGauInvCEN(c(muGIC,lambGIC),XMAT)
AICenGauInvC=-2*(PrSamIGC)+4

"Criterio Inform Akaike CON CENSURA:"
AICens=c(AICenExpC,AICenGamC,AICenNorC,AICenLNC,AICenGauInvC,AICenWeiC,AICenGumC)
names(AICens)=c("Exp","Gam","Nor","Lnor","GauInv","Wei","Gumbel")
cbind(sort(AICens))
##_______________________________________________________ ####
#---------------------------------------------- --
##  ............................................................................
# Gráficas Probab. CENSURA                     ####
# Matrices ordenadas de intervalos de censura de los datos:
XMatord=XMAT[order(XMAT[,1],decreasing=FALSE),]
XmatLNord=log(XMatord)
#Bandas 95% confianza para PP:
Ebetas=enes/(n+1)
varian=(enes*(n+1-enes))/(((n+1)^2)*(n+2))
## Joint 95% confidence band for N prediction intervals of uniform quantiles:
tau=1-(.05/n)
tau1=(1-tau)/2
tau2=(1+tau)/2
aenes= n+1-enes
Ban1=qbeta(tau1,enes,aenes)
Ban2=qbeta(tau2,enes,aenes)

#Para las estadísticas Anderson Darling se usa esto:
enes=seq(1,n,by=1)
aux1=(2*enes-1)
aux2=(2*n+1-(2*enes))

#--------------------------------- --
# Matrices datos transf para PP y QQ ####
unormat=pnorm(XMatord,mean=munorC,sd=signorC)
uLNmat=pnorm(XmatLNord,mean=muLNC,sd=sigLNC)
uexpmat=pexp(XMatord,rate=1/muExpC)
ugamat=pgamma(XMatord,shape=alfaGamC,scale=betagamC )
uWeimat=pweibull(XMatord,shape=betaWeiC,scale=sigWeiC)
uGuminmat=pGumbelmin(XMatord,c(amleGumC,bGumC))
uGauInvmat=pinvgauss(XMatord,mean=muGIC,shape=lambGIC)

#------------------------------------------------- --  
# PPCens Exponential distribution:              ####
# Simulating uniform values falling in each transformed interval:
# Three simulated samples with exact obs from estimated INv. GAussia distribution:
anchos=uexpmat[,2]-uexpmat[,1]
#first simulated sample:
uis=uexpmat[,1]+anchos*runif(n)
uexp1=sort(uis)
uis=uexpmat[,1]+anchos*runif(n)
uexp2=sort(uis)
uis=uexpmat[,1]+anchos*runif(n)
uexp3=sort(uis)
#Average weighted vertical distances to identity line are calculated for the 3 samples:
AVD1=sum(abs(uexp1-Ebetas)/sqrt(varian))/n
AVD2=sum(abs(uexp2-Ebetas)/sqrt(varian))/n
AVD3=sum(abs(uexp3-Ebetas)/sqrt(varian))/n
WAVDCExp=(AVD1+AVD2+AVD3)/3
#Anderson Darling:
AD1=-n-(1/n)*sum((aux1*log(uexp1))+aux2*log(1-uexp1))
AD2=-n-(1/n)*sum((aux1*log(uexp2))+aux2*log(1-uexp2))
AD3=-n-(1/n)*sum((aux1*log(uexp3))+aux2*log(1-uexp3))
ANDExp=(AD1+AD2+AD3)/3
#---------------------------------- --  
# Exponential Prob Plot for 3 simulated exact samples:
X11()
plot(Ebetas,uexp1,type="l",pch=19,cex=.5,ylab="F(xi)",
     xlab="U(0,1) Quantiles", xlim=c(0,1),ylim=c(0,1),
     main="Exponential Cens PP w 3 sim exact samples",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)
lines(Ebetas,uexp2,lty=1,col=3)
lines(Ebetas,uexp3,lty=1,col=4)
#------------------------------------ --
# PPCens Gamma distribution:              ####
# Simulating uniform values falling in each transformed interval:
# Three simulated samples with exact obs from estimated Normal distribution:
anchos=ugamat[,2]-ugamat[,1]
#first simulated sample:
uis=ugamat[,1]+anchos*runif(n)
ugam1=sort(uis)
uis=ugamat[,1]+anchos*runif(n)
ugam2=sort(uis)
uis=ugamat[,1]+anchos*runif(n)
ugam3=sort(uis)
#Average weighted vertical distances to identity line are calculated for the 3 samples:
AVD1=sum(abs(u1-Ebetas)/sqrt(varian))/n
AVD2=sum(abs(u2-Ebetas)/sqrt(varian))/n
AVD3=sum(abs(u3-Ebetas)/sqrt(varian))/n
WAVDCgam=(AVD1+AVD2+AVD3)/3
#Anderson Darling:
AD1=-n-(1/n)*sum((aux1*log(ugam1))+aux2*log(1-ugam1))
AD2=-n-(1/n)*sum((aux1*log(ugam2))+aux2*log(1-ugam2))
AD3=-n-(1/n)*sum((aux1*log(ugam3))+aux2*log(1-ugam3))
ANDGam=(AD1+AD2+AD3)/3


#---------------------------------- --  
# Gamma Prob Plot for 3 simulated exact samples:
X11()
plot(Ebetas,ugam1,type="l",pch=19,cex=.5,ylab="F(xi)",
     xlab="U(0,1) Quantiles", xlim=c(0,1),ylim=c(0,1),
     main="Gamma Cens PP w 3 sim exact samples",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)
lines(Ebetas,ugam2,lty=1,col=3)
lines(Ebetas,ugam3,lty=1,col=4)
#------------------------------------- --
# PPCens Normal distribution:       ####
# Simulating uniform values falling in each transformed interval:
# Three simulated samples with exact obs from estimated Normal distribution:
anchos=unormat[,2]-unormat[,1]
#first simulated sample:
uis=unormat[,1]+anchos*runif(n)
unor1=sort(uis)
uis=unormat[,1]+anchos*runif(n)
unor2=sort(uis)
uis=unormat[,1]+anchos*runif(n)
unor3=sort(uis)
#Average weighted vertical distances to identity line are calculated for the 3 samples:
AVD1=sum(abs(unor1-Ebetas)/sqrt(varian))/n
AVD2=sum(abs(unor2-Ebetas)/sqrt(varian))/n
AVD3=sum(abs(unor3-Ebetas)/sqrt(varian))/n
WAVDCnor=(AVD1+AVD2+AVD3)/3
#Anderson Darling:
AD1=-n-(1/n)*sum((aux1*log(unor1))+aux2*log(1-unor1))
AD2=-n-(1/n)*sum((aux1*log(unor2))+aux2*log(1-unor2))
AD3=-n-(1/n)*sum((aux1*log(unor3))+aux2*log(1-unor3))
ANDnor=(AD1+AD2+AD3)/3
#---------------------------------- --  
# Normal Probability Plot for 3 simulated exact samples:
X11()
plot(Ebetas,unor1,type="l",pch=19,cex=.5,ylab="F(xi)",
     xlab="U(0,1) Quantiles", xlim=c(0,1),ylim=c(0,1),
     main="Normal Cens PP w 3 sim exact samples",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)
lines(Ebetas,unor2,lty=1,col=3)
lines(Ebetas,unor3,lty=1,col=4)

#------------------------------------------------- --  
# PPCens Lognormal distribution:              ####
# Simulating uniform values falling in each transformed interval:
# Three simulated samples with exact obs from estimated Normal distribution:
anchos=uLNmat[,2]-uLNmat[,1]
#first simulated sample:
uis=uLNmat[,1]+anchos*runif(n)
uLN1=sort(uis)
uis=uLNmat[,1]+anchos*runif(n)
uLN2=sort(uis)
uis=uLNmat[,1]+anchos*runif(n)
uLN3=sort(uis)
#Average weighted vertical distances to identity line are calculated for the 3 samples:
AVD1=sum(abs(uLN1-Ebetas)/sqrt(varian))/n
AVD2=sum(abs(uLN2-Ebetas)/sqrt(varian))/n
AVD3=sum(abs(uLN3-Ebetas)/sqrt(varian))/n
WAVDCLN=(AVD1+AVD2+AVD3)/3
#Anderson Darling:
AD1=-n-(1/n)*sum((aux1*log(uLN1))+aux2*log(1-uLN1))
AD2=-n-(1/n)*sum((aux1*log(uLN2))+aux2*log(1-uLN2))
AD3=-n-(1/n)*sum((aux1*log(uLN3))+aux2*log(1-uLN3))
ANDLN=(AD1+AD2+AD3)/3


#---------------------------------- --  
# LogNormal Prob Plot for 3 simulated exact samples:
X11()
plot(Ebetas,uLN1,type="l",pch=19,cex=.5,ylab="F(xi)",
     xlab="U(0,1) Quantiles", xlim=c(0,1),ylim=c(0,1),
     main="LogNormal Cens PP w 3 sim exact samples",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)
lines(Ebetas,uLN2,lty=1,col=3)
lines(Ebetas,uLN3,lty=1,col=4)

#------------------------------------------------- --  
# PPCens Inverse Gaussian distribution:              ####
# Simulating uniform values falling in each transformed interval:
# Three simulated samples with exact obs from estimated INv. GAussia distribution:
anchos=uGauInvmat[,2]-uGauInvmat[,1]
#first simulated sample:
uis=uGauInvmat[,1]+anchos*runif(n)
uig1=sort(uis)
uis=uGauInvmat[,1]+anchos*runif(n)
uig2=sort(uis)
uis=uGauInvmat[,1]+anchos*runif(n)
uig3=sort(uis)
#Average weighted vertical distances to identity line are calculated for the 3 samples:
AVD1=sum(abs(uig1-Ebetas)/sqrt(varian))/n
AVD2=sum(abs(uig2-Ebetas)/sqrt(varian))/n
AVD3=sum(abs(uig3-Ebetas)/sqrt(varian))/n
WAVDCInvGau=(AVD1+AVD2+AVD3)/3
#Anderson Darling:
AD1=-n-(1/n)*sum((aux1*log(uig1))+aux2*log(1-uig1))
AD2=-n-(1/n)*sum((aux1*log(uig2))+aux2*log(1-uig2))
AD3=-n-(1/n)*sum((aux1*log(uig3))+aux2*log(1-uig3))
ANDGI=(AD1+AD2+AD3)/3
#---------------------------------- --  
# InvGaussian Prob Plot for 3 simulated exact samples:
X11()
plot(Ebetas,uig1,type="l",pch=19,cex=.5,ylab="F(xi)",
     xlab="U(0,1) Quantiles", xlim=c(0,1),ylim=c(0,1),
     main="Inverse Gaussian Cens PP w 3 sim exact samples",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)
lines(Ebetas,uig2,lty=1,col=3)
lines(Ebetas,uig3,lty=1,col=4)

#------------------------------------------------- --  
# PPCens Weibull distribution:              ####
# Simulating uniform values falling in each transformed interval:
# Three simulated samples with exact obs from estimated Normal distribution:
anchos=uWeimat[,2]-uWeimat[,1]
#first simulated sample:
uis=uWeimat[,1]+anchos*runif(n)
uwei1=sort(uis)
uis=uWeimat[,1]+anchos*runif(n)
uwei2=sort(uis)
uis=uWeimat[,1]+anchos*runif(n)
uwei3=sort(uis)
#Average weighted vertical distances to identity line are calculated for the 3 samples:
AVD1=sum(abs(uwei1-Ebetas)/sqrt(varian))/n
AVD2=sum(abs(uwei2-Ebetas)/sqrt(varian))/n
AVD3=sum(abs(uwei3-Ebetas)/sqrt(varian))/n
WAVDCWei=(AVD1+AVD2+AVD3)/3
#Anderson Darling:
AD1=-n-(1/n)*sum((aux1*log(uwei1))+aux2*log(1-uwei1))
AD2=-n-(1/n)*sum((aux1*log(uwei2))+aux2*log(1-uwei2))
AD3=-n-(1/n)*sum((aux1*log(uwei3))+aux2*log(1-uwei3))
ANDWei=(AD1+AD2+AD3)/3
#---------------------------------- --  
# Weibull Prob Plot for 3 simulated exact samples:
X11()
plot(Ebetas,uwei1,type="l",pch=19,cex=.5,ylab="F(xi)",
     xlab="U(0,1) Quantiles", xlim=c(0,1),ylim=c(0,1),
     main="Weibull Cens PP w 3 sim exact samples",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)
lines(Ebetas,uwei2,lty=1,col=3)
lines(Ebetas,uwei3,lty=1,col=4)

#------------------------------------------------- --  
# PPCens Gumbel minima distribution:              ####
# Simulating uniform values falling in each transformed interval:
# Three simulated samples with exact obs from estimated Normal distribution:
anchos=uGuminmat[,2]-uGuminmat[,1]
#first simulated sample:
uis=uGuminmat[,1]+anchos*runif(n)
ugum1=sort(uis)
uis=uGuminmat[,1]+anchos*runif(n)
ugum2=sort(uis)
uis=uGuminmat[,1]+anchos*runif(n)
ugum3=sort(uis)
#Average weighted vertical distances to identity line are calculated for the 3 samples:
AVD1=sum(abs(ugum1-Ebetas)/sqrt(varian))/n
AVD2=sum(abs(ugum2-Ebetas)/sqrt(varian))/n
AVD3=sum(abs(ugum3-Ebetas)/sqrt(varian))/n
WAVDCGumin=(AVD1+AVD2+AVD3)/3
#Anderson Darling:
AD1=-n-(1/n)*sum((aux1*log(ugum1))+aux2*log(1-ugum1))
AD2=-n-(1/n)*sum((aux1*log(ugum2))+aux2*log(1-ugum2))
AD3=-n-(1/n)*sum((aux1*log(ugum3))+aux2*log(1-ugum3))
ANDGum=(AD1+AD2+AD3)/3

#---------------------------------- --  
# Gumbel minCEn Prob Plot for 3 simulated exact samples:
X11()
plot(Ebetas,ugum1,type="l",pch=19,cex=.5,ylab="F(xi)",
     xlab="U(0,1) Quantiles", xlim=c(0,1),ylim=c(0,1),
     main="Gumbel Min Cens PP w 3 sim exact samples",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(Ebetas,Ban1,lty=2,col=2)
lines(Ebetas,Ban2,lty=2,col=2)
lines(Ebetas,ugum2,lty=1,col=3)
lines(Ebetas,ugum3,lty=1,col=4)
#------------------------------------------------- --  
##_______________________________________________________ ####
# WAVD CON CENSURA   ####
# WAVD= average weighted Vertical distances of u(k) points to E(uk)in PPlot:
WAVDCens=c(WAVDCExp,WAVDCgam,WAVDCnor,WAVDCLN,WAVDCInvGau,WAVDCWei,WAVDCGumin)
names(WAVDCens)=c("Exp","Gam","Nor","Lnor","GausInv","Wei","Gumbel")
cbind(sort(WAVDCens))
# Anderson Darling CENSURA   ####
ANDCens=c(ANDExp,ANDGam,ANDnor,ANDLN,ANDGI,ANDWei,ANDGum)
names(ANDCens)=c("Exp","Gam","Nor","Lnor","GausInv","Wei","Gumbel")
"Anderson Darling Prom 3 sim samples:"
cbind(sort(ANDCens))

##_______________________________________________________ ####
# Gráficas QQ CON CENSURA:         ####
# En eje horizontal van los cuantiles estimados y en el vertical los 
# datos ordenados que son los cuantiles empíricos correspondientes.
# Se transforman de escala las gráficas de Probab. con censura obtenidas antes.
#--Los sig. son los cuantiles teóricos bajo la distribución estimada
# van en el eje horizontal y se comparan contra los cuantiles empíricos 
# en el eje vertical.
# Para los Qp teóricos se usa la distribución inversa estimada
# de las 7 distribuciones consideradas aplicada a Ebetas:
FinvExpC=-muExpC*log(1-Ebetas)
FinvGamC=qgamma(Ebetas,shape=alfaGamC,scale=mugamC/alfaGamC)
FinvNorC=qnorm(Ebetas,mean=munorC,sd=signorC)
FinvLNC=qlnorm(Ebetas,meanlog=muLNC,sdlog=sigLNC)
FinvGIC=qinvgauss(Ebetas,mean=muGIC,shape=lambGIC)
FinvWeiC=exp(log(sigWeiC)+log(-log(1-Ebetas))/betaWeiC)
#La dist. inv. Gumbel de mínimos para los datos es x=a+ln(-ln(1-u))/beta
FinvGumC=amleGumC+log(-log(1-Ebetas))*bGumC
#------------------------------- --

miniX=XMAT[1,1]
#maxiX=XMAT[n,2]
maxiX=XMAT[n,1]+10
# Se transforman las bandas de confianza Beta a la escala
# de cada una de las 7 distribuciones estim con censura:
Fban1ExpC=-muExpC*log(1-Ban1)
Fban2ExpC=-muExpC*log(1-Ban2)
Fban1GamC=qgamma(Ban1,shape=alfaGamC,scale=mugamC/alfaGamC)
Fban2GamC=qgamma(Ban2,shape=alfaGamC,scale=mugamC/alfaGamC)
Fban1NorC=qnorm(Ban1,mean=munorC,sd=signorC)
Fban2NorC=qnorm(Ban2,mean=munorC,sd=signorC)
Fban1LNC=qlnorm(Ban1,meanlog=muLNC,sdlog=sigLNC)
Fban2LNC=qlnorm(Ban2,meanlog=muLNC,sdlog=sigLNC)
Fban1GIC=qinvgauss(Ban1,mean=muGIC,shape=lambGIC)
Fban2GIC=qinvgauss(Ban2,mean=muGIC,shape=lambGIC)
Fban1WeiC=exp(log(sigWeiC)+log(-log(1-Ban1))/betaWeiC)
Fban2WeiC=exp(log(sigWeiC)+log(-log(1-Ban2))/betaWeiC)
Fban1GumC=amleGumC+log(-log(1-Ban1))*bGumC
Fban2GumC=amleGumC+log(-log(1-Ban2))*bGumC

#------------------------------------------- --
#QQ Exponencial con Censura:  ####
FinvExpC1=-muExpC*log(1-uexp1)
FinvExpC2=-muExpC*log(1-uexp2)
FinvExpC3=-muExpC*log(1-uexp3)
X11() 
plot(FinvExpC,FinvExpC1,pch=19,cex=.5,ylab="Empirical Quantiles",
     xlab="Estim.Quantiles", xlim=c(miniX,maxiX),ylim=c(miniX,maxiX),
     main="Exponencial gráfica QQCen",pty="m")
lines(c(miniX,maxiX),c(miniX,maxiX),lty=2,col=2)
lines(FinvExpC,Fban1ExpC,lty=2,col=2)
lines(FinvExpC,Fban2ExpC,lty=2,col=2)
lines(FinvExpC,FinvExpC2,lty=1,col=3)
lines(FinvExpC,FinvExpC3,lty=1,col=4)
#-------------------------- --
#QQ Gamma con Censura:  ####
FinvGamC1=-muExpC*log(1-uexp1)
FinvGamC2=-muExpC*log(1-uexp2)
FinvGamC3=-muExpC*log(1-uexp3)
X11() 
plot(FinvGamC,FinvGamC1,pch=19,cex=.5,ylab="Empirical Quantiles",
     xlab="Estim.Quantiles", xlim=c(miniX,maxiX),ylim=c(miniX,maxiX),
     main="Gamma gráfica QQCen",pty="m")
lines(c(miniX,maxiX),c(miniX,maxiX),lty=2,col=2)
lines(FinvGamC,Fban1GamC,lty=2,col=2)
lines(FinvGamC,Fban2GamC,lty=2,col=2)
lines(FinvGamC,FinvGamC2,lty=1,col=3)
lines(FinvGamC,FinvGamC3,lty=1,col=4)
#-------------------------- --
#QQ Normal con Censura:  ####
FinvNorC1=qnorm(unor1,mean=munorC,sd=signorC)
FinvNorC2=qnorm(unor2,mean=munorC,sd=signorC)
FinvNorC3=qnorm(unor3,mean=munorC,sd=signorC)
X11() 
plot(FinvNorC,FinvNorC1,pch=19,cex=.5,ylab="Empirical Quantiles",
     xlab="Estim.Quantiles", xlim=c(miniX,maxiX),ylim=c(miniX,maxiX),
     main="Normal gráfica QQCen",pty="m")
lines(c(miniX,maxiX),c(miniX,maxiX),lty=2,col=2)
lines(FinvNorC,Fban1NorC,lty=2,col=2)
lines(FinvNorC,Fban2NorC,lty=2,col=2)
lines(FinvNorC,FinvNorC2,lty=1,col=3)
lines(FinvNorC,FinvNorC3,lty=1,col=4)

#-------------------------- --
#QQ Lognormal con Censura:  ####
FinvLNC1=qlnorm(uLN1,mean=muLNC,sd=sigLNC)
FinvLNC2=qlnorm(uLN2,mean=muLNC,sd=sigLNC)
FinvLNC3=qlnorm(uLN3,mean=muLNC,sd=sigLNC)
X11() 
plot(FinvLNC,FinvLNC1,pch=19,cex=.5,ylab="Empirical Quantiles",
     xlab="Estim.Quantiles", xlim=c(miniX,maxiX),ylim=c(miniX,maxiX),
     main="Lognormal gráfica QQCen",pty="m")
lines(c(miniX,maxiX),c(miniX,maxiX),lty=2,col=2)
lines(FinvLNC,Fban1LNC,lty=2,col=2)
lines(FinvLNC,Fban2LNC,lty=2,col=2)
lines(FinvLNC,FinvLNC2,lty=1,col=3)
lines(FinvLNC,FinvLNC3,lty=1,col=4)

#-------------------------- --
#QQ Gaussiana Inversa con Censura:  ####
FinvGIC1=qinvgauss(uig1,mean=muGIC,shape=lambGIC)
FinvGIC2=qinvgauss(uig2,mean=muGIC,shape=lambGIC)
FinvGIC3=qinvgauss(uig3,mean=muGIC,shape=lambGIC)
X11() 
plot(FinvGIC,FinvGIC1,pch=19,cex=.5,ylab="Empirical Quantiles",
     xlab="Estim.Quantiles", xlim=c(miniX,maxiX),ylim=c(miniX,maxiX),
     main="Gaussiana Inversa gráfica QQCen",pty="m")
lines(c(miniX,maxiX),c(miniX,maxiX),lty=2,col=2)
lines(FinvGIC,Fban1GIC,lty=2,col=2)
lines(FinvGIC,Fban2GIC,lty=2,col=2)
lines(FinvGIC,FinvGIC2,lty=1,col=3)
lines(FinvGIC,FinvGIC3,lty=1,col=4)

#-------------------------- --
#QQ Weibull con Censura:  ####
FinvWeiC1=exp(log(sigWeiC)+log(-log(1-uwei1))/betaWeiC)
FinvWeiC2=exp(log(sigWeiC)+log(-log(1-uwei2))/betaWeiC)
FinvWeiC3=exp(log(sigWeiC)+log(-log(1-uwei3))/betaWeiC)
X11() 
plot(FinvWeiC,FinvWeiC1,pch=19,cex=.5,ylab="Empirical Quantiles",
     xlab="Estim.Quantiles", xlim=c(miniX,maxiX),ylim=c(miniX,maxiX),
     main="Weibull gráfica QQCen",pty="m")
lines(c(miniX,maxiX),c(miniX,maxiX),lty=2,col=2)
lines(FinvWeiC,Fban1WeiC,lty=2,col=2)
lines(FinvWeiC,Fban2WeiC,lty=2,col=2)
lines(FinvWeiC,FinvWeiC2,lty=1,col=3)
lines(FinvWeiC,FinvWeiC3,lty=1,col=4)

#-------------------------- --
#QQ Gumbel con Censura:  ####
FinvGumC1=amleGumC+log(-log(1-ugum1))*bGumC
FinvGumC2=amleGumC+log(-log(1-ugum2))*bGumC
FinvGumC3=amleGumC+log(-log(1-ugum3))*bGumC
X11() 
plot(FinvGumC,FinvGumC1,pch=19,cex=.5,ylab="Empirical Quantiles",
     xlab="Estim.Quantiles", xlim=c(miniX,maxiX),ylim=c(miniX,maxiX),
     main="Gumbel gráfica QQCen",pty="m")
lines(c(miniX,maxiX),c(miniX,maxiX),lty=2,col=2)
lines(FinvGumC,Fban1GumC,lty=2,col=2)
lines(FinvGumC,Fban2GumC,lty=2,col=2)
lines(FinvGumC,FinvGumC2,lty=1,col=3)
lines(FinvGumC,FinvGumC3,lty=1,col=4)


#_________________________________________ ####

# Graficar densidades estimadas                                           ####
# Se grafican las mejores densidades con y sin censura:
miniX=XMAT[1,1]
#maxiX=XMAT[n,2]
maxiX=XMAT[n,1]+10
AA=miniX
BB=maxiX
AA=0.01
BB=220
xs=seq(AA,BB,length.out=200)
#---- SIN CENSURA: ####
#Exponencial
yexp=exp(-xs/mumle)/mumle
#Gamma:
ygam=dgamma(xs,shape=alphamle,scale=betagamle)
#Normal
ynor=dnorm(xs,mumle,sigmle)
#Lognormal con censura:
yLN=dlnorm(xs,meanlog=muLN,sdlog=sigLN)
#Gaussiana Inversa:
#tras bajar paquete {actuar} esta es también desidad Gaussiana INversa:
yGI= dinvgauss(xs,mean=mumle,shape=lambmle)
#Weibull:
aux=betaWei*(log(xs)-amle)
ywei=(betaWei/xs)*exp(aux-exp(aux))
#Gumbel:
auxG=betaGum*((xs)-amleGum)
ygum=betaGum*exp(auxG-exp(auxG))
#---- CON CENSURA:  ####
#Exponencial
yExpC=exp(-xs/muExpC)/muExpC
#Gamma:
ygamC=dgamma(xs,shape=alfaGamC,scale=mugamC/alfaGamC)
#Normal
ynorC=dnorm(xs,mean=munorC,sd=signorC)
#Lognormal con censura:
yLN=dlnorm(xs,meanlog=muLNC,sdlog=sigLNC)
#Gaussiana Inversa:
#tras bajar paquete {actuar} esta es también desidad Gaussiana INversa:
yGI= dinvgauss(xs,mean=muGIC,shape=lambGIC)
c(muGIC,lambGIC)
#Weibull:
aux=betaWeiC*(log(xs)-log(sigWeiC))
yweiC=(betaWeiC/xs)*exp(aux-exp(aux))
#Gumbel:
auxG=((xs)-amleGumC)/bGumC
ygumC=exp(auxG-exp(auxG))/bGumC

#------------------- --
#tope altura gráfica:
densitop=max(c(yexp,ygam,ynor,yLN,yGI,ywei,ygum))

#--------------------- --
# PARA RATONES DE GRICE BAIN
#Para datos ratones de Grice-Bain fueron Gumbel Weibull y Normal
X11() 
plot(xs,ygum,ylab="f(x)",yaxs="i",
     xlab="x", xlim=c(xs[1],xs[200]),
     ylim=c(0,densitop),
     main="f(x) cens:Gumbel negro, Weibull rojo, Norm azul",type="l")
points(X,rep(0,n),type="p",col="6",pch=19,cex=1)
lines(xs,ywei,lty=1,col=2)#red solid
lines(xs,ynor,lty=2,col=4)#blue dashes
#van datos exactos 20 ratones machos:
points(xsort,rep(0,n),type="p",col="4",pch=19,cex=1)

#----------------------------------------- --
# DATOS CRUAPAN:
# Para Cruapán máximos las mejores fueron la Gamma y la Lognormal
X11() 
plot(xs,ygam,ylab="f(x)",yaxs="i",
     xlab="x", xlim=c(xs[1],xs[200]),
     ylim=c(0,densitop),
     main="f(x) cens:Gumbel negro, Weibull rojo, Norm azul",type="l")
points(X,rep(0,n),type="p",col="6",pch=19,cex=1)
lines(xs,yLN,lty=2,col=2)#red solid
ZZ=table(xsort) #datos exactos registrados
zuni=unique(xsort) #valores únicos
#van datos únicos con sus repeticiones de altura controlada:
points(zuni,ZZ*(0.01/25),type="p",col="4",pch=19,cex=1)


#_________________________________________ ####
#  Output file with main results:                 ####
sink(file="ConySinCensura.txt")

"EMVS Sin y Con CENSURA:"
"Exponencial media exacta y con censura:"
c(mumle,muExpC)

"Gamma forma y media exacta y con censura:"
c(alphamle,mumle)
c(alfaGamC,mugamC)
"Gamma escala exacta y con censura:"
c(mumle/alphamle,mugamC/alfaGamC)

"Normal emvs mu,sigma exactos y con cesura:"
c(mumle,sigmle)
c(munorC,signorC)

"Logmormal emvs de mu,sigma exactos y con censura:"
c(muLN,sigLN)
c(muLNC,sigLNC)

"Gaussiana Inversa media y forma, exactos y con censura:"
c(mumle,lambmle)
c(muGIC,lambGIC)

"Weibull shapebeta, scalesigma mles exact and cens:"
c(betaWei,exp(amle))
c(betaWeiC,sigWeiC)

"Gumbel minima, agum location, bgum scale mles exact and cens:"
c(amleGum,deltamle)
c(amleGumC,bGumC)

# ---------------------- ---
# AKAIKE SIN CENSURA:
"AIC sin censura PARA COMPARAR:"
cbind(sort(AICSIN))
"AIC CON CENSURA:"
cbind(sort(AICens))
#--------------------- --
"WAVD CON CENSURA:"
cbind(sort(WAVDCens))
"ANDERSON DARLING CON CENSURA"
cbind(sort(ANDCens))

"E(X) para Exponencial, Gamma,Normal y GaussInv es:"
c(muExpC,mugamC,munorC,muGIC)
"E(x) para la Lognormal es"
Elognor=exp(muLNC +(sigLNC^2)/2)
Elognor
"E(x) para la Weibull es"
EWei=sigWeiC*gamma(1+(1/betaWeiC))
EWei
"E(X) para la Gumbel es:"
EGum = amleGumC + bGumC*(0.5772156)
EGum
sink(file=NULL)

# para ver relación entre R(teta0) y el pvalor asociado.
#----------------------------- --
Rtets=seq(0,1,length.out=200)
Z=-2*log(Rtets)
pvals=1-pchisq(Z,1)
zx=-2*log(0.1465)
pz=1-pchisq(zx,1)
cbind(Rtets,pvals)

X11() 
plot(Rtets,pvals,type="l",col=1,lwd=2,
     main="R(teta0) vs pvalores asociados",xlab="R(teta0)",ylab="pvalue",
     xaxs="i",yaxs="i",xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),lty=2,col=2)
lines(c(0.1465,0.1465),c(0,1),lty=2,col=4)
axis(side=1,at=seq(0,1,by=0.1))
axis(side=2,at=seq(0,1,by=0.1))
