
#   ____________________________________________________________________________
#       INFERENCIAS PARA LA NORMAL                                          ####

#      Elaborado por Eloísa Díaz Francés, ABRIL de 2023
#  ------------------------------------- --
#    Este programa realiza los siguientes pasos:
#  1) Simula muestra normal de tamaño N O CARGA DATOS NORMALES
#  2) Calcula los emv de mu,sigma analíticamente
#  3) Calcula numéricamente los emv con OPTIM
#  4) Calcula y grafica Rp(mu) analítica 
#     con intervaloS propuestos de verosimilitud y confianza de nivel c
#  5) Ahora calcula Rp(mu) numéricamente usando ciclo For y 
#     OPTIMIZE porque sólo hay un parámetro de estorbo PARA COMPARAR.
#  6) Grafica regiones de verosimilitud y confianza para MU,SIGMA
#  7) Calcula y grafica Rp(sigma) con intervalos propuestos.
#  8) Valida el modelo normal por cuatro maneras.
#  9) Se muestra cómo partir pantallas en secciones para mostrar gráficas.
#  10) Se despliegan resultados en archivo.

##  ............................................................................
##         Limpieza del espacio de trabajo:                                ####

#   Conviene limpiar al inicio, todo el espacio de trabajo:
#   Equivale a limpiar el "environment" con la escobita.

rm(list=ls(all=TRUE))

##  ............................................................................
##      Funciones usadas:                                                   ####

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
###   LogVer(mu,sigma)                                             ####

#   Se calcula la logverosimilitud normal en [mu,sigma].
#   Los argumentos son el vector de dos parámetros, las estadísticas 
#   suficientes t1,t2 y el tamaño de muestra n.
lvnor=function(vec2,t1,t2,n){
   mu=vec2[1]
  sigma=vec2[2]
  if (sigma >0){
       s2=sigma^2
    lv=-n*log(sigma)-(t2/(2*s2))+(mu*t1/s2)-((n*mu^2)/(2*s2))
  }else{
    lv=-99999999999999999}
  return(lv)
  } 
### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
###  Cálculo de lp(mu) numérico:                                           ####
#   lpmu es idéntica que l(mu,sigma) cambiando argumentos al reducir en uno el 
#   vector parámetros. Así sólo es función de sigma, parámetro de estorbo, de
#   las estadísticas suficientes, de n y del valor fijo de MU en mufix.

lpmu=function(vec1,t1,t2,n,mufix){
  mu=mufix
  sigma=vec1
    if (sigma >0){
       s2=sigma^2
       lv=-n*log(sigma) - (t2-(2*t1*mu) + (n*(mu^2)))/(2*s2)
      }else{
     lv=-99999999999999999
    }
  return(lv)
  } 

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
###   rp(mu) analítica                               ####

rpmu=function(mu,t1,t2,n){
  HH=t2-(t1^2)/n
  muhat=t1/n
  Hmu=t2-(2*t1*mu)+(n*(mu^2))
  rp=-(n/2)*(log(Hmu)-log(HH))
  return(rp)
  }  

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
###   rp(sigma)-lnc para calcular interv perfil sigma numéricos ####

rpsigC=function(sig,SIGMV,NN,CC){
  lnC=log(CC)
  if (sig>0){
  rpsig=NN*log(SIGMV/sig)-(NN/2)*((SIGMV/sig)**2)+ (NN/2)
  rpc=rpsig-lnC
  } else {
    rpc=999
  }
return(rpc)
}


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
###       Niveles Verosimilitud: c(n) para n>=5                      ####
#  se calcula el nivel de verosimilitud como función de n
#  para intervalos de verosimilitud perfil de mu y de sigma normales
#  con confianzas del 90, 95 y 99%
  clevels=function(n){
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
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  . . ..
  ###       SIN NOMBRES                                       #### 
 # Función que ayuda a generar objeto idéntico sin nombres para facilitar 
 #  despliegue de información final:
  
sinnombre=function(object1){
  object2=object1
  names(object2)=NULL
  return(object2)
  }

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  . . ..
###   Intervalos Perfil Sigma numéricos                       #### 
#Con confianza del 90, 95 o 99%
#--------------------------- --

IntSig=function(SIGMV,NN,confi){
  #van niveles verosim para sigma:
  cpp1s=clevels(NN)$csig
  if (confi==0.90){
    cp=cpp1s[1]
  } else if (confi==0.95) {
    cp=cpp1s[2]
  } else if (confi==0.99){
    cp=cpp1s[3]
  } else {
    cat( "Mala confianza. Se calculará int sigma del 95%: ")
    print(noquote("")) 
    cp=cpp1s[2]
  }
  s11=uniroot(rpsigC,c(0.00001,SIGMV),SIGMV,NN,cp,
              lower=0.00001,upper=SIGMV,tol=0.000001)$root
  s12=uniroot(rpsigC,c(SIGMV,SIGMV*10),SIGMV,NN,cp,
              lower=SIGMV,upper=SIGMV*10,tol=0.000001)$root
   intersig=c(s11,s12,confi)
  names(intersig)=c("Sig1","Sig2","Confianza")
  return(intersig)
}


#  _________________________________________________________-___________________
#               PROGRAMA PRINCIPAL                                          ####

##  ............................................................................
##              Paso 1 DATOS                                                ####
##  ............................................................................
##              1a SIMULADOS   s                                           ####

#  SIMULADOS: 
#  Se simulan N normales con media mu y desviación estándar sigma
# N=50
# mu=10
# sigma=4
# #Los datos normales simulados se guardan en:
# DAT=rnorm(N,mean=mu,sd=sigma)

# ## se simulan datos Gamma
# gamiform=2
# gamedia=mu #la media Gamma es producto forma*escala
# gamscale=gamedia/gamiform
# DAT=rgamma(N,shape=gamiform,scale=gamscale)



##  ............................................................................
##              1b Peso Ajos                                           ####

#  Peso de una cabeza de ajo medida por alumnos DEMAT 2022 y 2023
# 
DAT=c(29,58,31,56,32,59,51,71,32,78,30,59,67, 69,40,79,72,51,12,44,
      83,32,36,37,52,34,51,72,17,7,31,60,52,48,56,102,42,62,39, 47,
      41,54,50,49)
##  ............................................................................
##              1c LOG peso AJOS                                        ####
##   (Comentar si no se usa:) 
# DAT=log(DAT)

##  ............................................................................
##              1d Núm. dientes AJo                                        ####
## ESTOS DATOS NO SON NORMALES PERO SIRVEN PARA VALORAR METODOS VALIDACION
# 
DAT=c(10,21,21,26,14,10,11,14,14,10,22,17,23,14,14,12,8,22,28,11,
      24,11,25,15,13,16,7,20,5,15,13,7,15,24,16,16,26,13,26,14,14,12)
#para comparar:
DAT=rnorm(25,mean=-0.8,sd=2)

##  ............................................................................
##              1e Bacterias por litro de agua en lago                      ####
##   Bhatacharyya y Johnson (1977, p. 261)
DAT=c(173,190,215,198,184,207,210,193,196,180)

##  ............................................................................
##              1f Gruesos de 10 láminas de plástico                              ####
##       Bhatacharyya y Johnson (1977, p. 281)

DAT=c(226,228,226,225,232,228,227,229,225,230)

##  ............................................................................
##       Paso 2 Emvs analíticos                                        ####

N=length(DAT)
datsort=sort(DAT)
#van estadísticas suficientes
t1=sum(DAT)
t2=sum(DAT^2)
#suma de cuadrados
H=t2-(t1^2)/N
#Se calculan los emv de mu y sigma^2 analíticamente:
mumv=t1/N
sig2mv=H/N
sigmamv=sqrt(sig2mv)

##  ............................................................................
##         Paso 3 Emvs numéricos                                          ####

# Van los emv con OPTIM con el método default de optimización, Nelder Mead. 
#  La salida de esta función es una lista:
#   el primer componente es el vector de emvs, donde alcanza el máximo la LV
#   el segdo comp. es el valor de la función en el emv
#   el tercer comp. es el número de llamadas a la logver.

vecini=c(mumv,sigmamv) #valores iniciales de mu y sigma
maxversali=optim(vecini,lvnor,gr="Null",t1,t2,N,method="Nelder-Mead",
          control=list(fnscale=-2) )

#van los emv (deben ser iguales a los analíticos):
emvs=maxversali$par
mumvnum=emvs[1]
sigmvnum=emvs[2]
#va el valor de la logver en este máximo:
LVmax=maxversali$value
#van número de iteraciones del proc numérico:
iteras=maxversali$counts

#va despliegue resultados en pantalla:
m="Van primero los emv analíticos de mu y sigma:"
cat(mumv,sigmamv) ##concatenate and print: cat()

m="Ahora van los emv numéricos con Optim:"
cat(mumvnum,sigmvnum)
m="Ahora va la logver en el emv y las iteraciones (llamados a logver):"
cat(LVmax,iteras[1])
##  ............................................................................
##        Paso 4 Rp(mu) e intervalos                                        ####

##  Se calcula y grafica Rp(mu) (analíticamente calculada) con intervalo
#   de verosimilitud de nivel c. Va la rejilla donde se grafica:
MM=150
mu1=mumv-1.5*sigmamv
mu2=mumv+1.5*sigmamv
mus=seq(mu1,mu2,length=MM)

#se calcula la perfil relativa en el vector de mus

rpmuss=sapply(X=mus,FUN=rpmu,t1,t2,N)
#se grafica la perfil relativa de mu en pantalla nueva:
#En RStudio sirve para esto el comando X11()

#se calculan los niveles de verosimilitud de intervalos con 90, 95 y 
# 99% de confianza
cesmu=clevels(N)$cmu
#cesmu sin nombres:
CesMuS=sinnombre(cesmu)
aux=sigmamv*sqrt((CesMuS**(-2/N))-1)
Intmus1=mumv-aux
Intmus2=mumv+aux
Intermu=cbind(Intmus1,Intmus2)
# para desplegar resultados al final:
confianzas=c("90%","95%","99%")
Intermus=cbind(Intermu,CesMuS,confianzas)

X11()
plot(mus,exp(rpmuss), type="l",col=1,lwd=1,
     main="Rp(mu) analítica",
     xlab=expression(mu),ylab=expression(R[p](mu)),xaxs="i",yaxs="i",
     ylim=c(0,1))
cvec=cesmu*cbind(rep(1,3),rep(1,3))
lines(Intermu[1,],cvec[1,],lty=1,col=2)
lines(Intermu[2,],cvec[2,],lty=1,col=2)
lines(Intermu[3,],cvec[3,],lty=1,col=2)

##  ............................................................................
##       Paso 5 Rp(mu) numérica                                             ####

# ------ Se calcula y grafica la perfil de mu pero ahora numéricamente
##      Para compararla y ratificar que son indistinguibles
#Rejilla regular para mus donde se graficará la perfil:

lvpmu=rep(1,MM) #logverosim perfil para la mu fija en mus
sigrestmv=rep(0,MM) #los emv restringidos de sigma para esas mu fijas

for(j in (1:MM)){ 
      interv=c(0.001,sigmamv*15)  #int de búsqueda para sigma
              MUFIJO=mus[j]
    auxsali=optimize(lpmu,interval=interv,t1,t2,N,MUFIJO,maximum=TRUE,
                     lower=min(interv),upper=max(interv),tol=0.0000000000001)             
    lvpmu[j]=auxsali$objective
    sigrestmv[j]=auxsali$maximum}

#se calcula verosimilitud relativa:
Rpmunum=exp(lvpmu-lpmu(sigmamv,t1,t2,N,mumv))

#se grafica la nueva perfil de mu numérica sobre la exacta para compararlas: 

X11()
plot(mus,Rpmunum, type="l",col=1,lwd=1,
     main="Rp(mu) numérica negro y analítica rojo (guiones)",
     xlab=expression(mu),ylab=expression(R[p](mu)),xaxs="i",yaxs="i",
     yaxp=c(0,1,10),ylim=c(0,1),xaxp=c(-4,6,20))
points(mus,exp(rpmuss),lty=2, type="l",col=2) 

##  ............................................................................
##      Paso 6  Regiones Verosimilitud (mu,sigma)                                                   ####
#------ Contornos de la verosimilitud relativa 
#se crea matriz aprovechando rejilla de mus y ahora
# se hace para sigmas, donde se guardar? la logverosim de mu,sigma

NN=50
mu1=mumv-1.5*sigmamv
mu2=mumv+1.5*sigmamv
mucont=seq(mu1,mu2,length=NN) 
sig1=sigmamv/2
sig2=sigmamv*2
sigis=seq(sig1,sig2,length=NN)
lvmatriz=matrix(rep(0,NN*NN),nrow=NN,ncol=NN)
#se calcula la l(mu,sigma) en esta rejilla matricial:
for(i in 1:NN){
  for(j in 1:NN){
      veci2=c(mucont[i],sigis[j])
    lvmatriz[i,j]=lvnor(veci2,t1,t2,N)
  }
}
 #Ahora va la logverosimilitud relativa en misma rejilla:
LVmax=lvnor(c(mumv,sigmamv),t1,t2,N)
Rvmatriz=exp(lvmatriz-LVmax)


#va gráfica de regiones de verosimilitud y confianza:
X11()
cesmusig=clevels(N)$cmusig
contour(mucont,sigis,Rvmatriz,
        #levels=c(0.05,0.15,0.25,0.5),
        #labels=c(0.05,0.15,0.25,0.5),
        levels=cesmusig,labels=names(cesmusig),
        labcex=1,xlim=c(mu1,mu2),ylim=c(sig1,sig2),
        xlab=expression(mu),
        ylab=expression(sigma),
        cex.lab=1,
        cex.axis=1 )
#se marca el emv sobre contornos
points(mumv,sigmamv,col=2,pch=8) 
#se marca la trayectoria de la perfil de mu sobre los contornos:
points(mus,sigrestmv,cex=0.3,col=4,type="l")
points(c(mumv,mumv),c(sig1,sig2),cex=0.3,col=5,type="l")

# OJO: COMO PUEDEN MARCAR LAS TRAYECTORIAS DE LAS DOS PERFILES
# SOBRE ESTOS CONTORNOS?  PUEDEN USAR LINES


##  ............................................................................
##       Paso 7 Rp(sigma) e intervalos                               ####
##   Intervalos para Sigma de verosimilitud perfil y confianza
cessig=clevels(N)$csig
aux=sqrt(-4*log(cessig)/N)/3
sig1=sigmamv*((1+aux)**(-3/2))
sig2=sigmamv*((1-aux)**(-3/2))
Intersig=cbind(sig1,sig2)
Intersigs=cbind(Intersig,cessig,confianzas)
S1=sigmamv/2
S2=sigmamv*2
SIGS=seq(S1,S2,length.out=200)
rpsig=N*log(sigmamv/SIGS)-(N/2)*((sigmamv/SIGS)**2)+ (N/2)
Cesigs=sinnombre(cessig)
Intersigs=cbind(Intersig,Cesigs,confianzas)

## Se calculan aquí los intervalos de verosim perfil de mismos
## niveles cessig para sigma pero de manera numérica para compararlos:
InumSig90=IntSig(sigmamv,N,0.90)
InumSig95=IntSig(sigmamv,N,0.95)
InumSig99=IntSig(sigmamv,N,0.99)
IntSigsNum=rbind(InumSig90,InumSig95,InumSig99)


X11()
plot(SIGS,exp(rpsig), type="l",col=1,lwd=1,
     main="Rp(sigma)",
     xlab=expression(sigma),ylab=expression(R[p](sigma)),xaxs="i",yaxs="i",
     ylim=c(0,1))
cvec=cessig*cbind(rep(1,3),rep(1,3))
lines(Intersig[1,],cvec[1,],lty=1,col=2)
lines(Intersig[2,],cvec[2,],lty=1,col=2)
lines(Intersig[3,],cvec[3,],lty=1,col=2)


##  ............................................................................
##    Paso 8   Validación Modelo Normal                             ####

##     8.1) Histograma vs densidad normal estimada
#      8.2) Densidad estimada y datos 
#      8.3) Distribución empirica vs estimada
#      8.4) Mi gráfica de probabilidad


#---- 8.1 Histograma y densidad estimada:--------------
X11()
hresul=hist(DAT,freq=TRUE,xlab="Datos",ylab="Frecuencias",
            main="Histograma con densidad normal estimada")
#ptos para graficar densidad normal encimada y encimar al histograma:
puntos=seq(min(hresul$breaks),max(hresul$breaks),length=25)
#se evalúan y reescalan densidades multiplicando por el número de 
# datos n y ancho barra hh: 
hh=diff(hresul$mids[1:2]) #ancho barra del histograma
nordensheight=dnorm(puntos, mean = mumv, sd =sigmamv)*N*hh
lines(puntos,nordensheight,type="l",col=2)


#---- 8.2  Gráfica Densidad Estimada y datos ------------------
X11()
x1=mumv-3*sigmamv
x2=mumv+3*sigmamv
densitop=1/(sqrt(2*pi)*sigmamv)
xs=seq(x1,x2, length.out=200)
#Lo sig. es equivalente pero más largo:
#ynor=sapply(X=xs,FUN=dnorm,mumv,sigmamv)
ynor=dnorm(xs,mumv,sigmamv)
plot(xs,ynor,ylab="f(x;mu,sigma)",yaxs="i",
     xlab="x", xlim=c(x1,x2),ylim=c(0,densitop),
     main="Densidad normal estimada y datos",type="l")
points(DAT,rep(0,N),type="p",col="6",pch=19,cex=1)


#-----  8.3 Distribución empirica vs estimada: --------------
X11()
mini=floor(datsort[1])
maxi=ceiling(datsort[N])
# type="s" grafica escaleritas:
plot(datsort,(1:N)/(N+1),type="s",ylim=c(0,1),
     main="Distribución empírica vs estimada", xaxp=c(mini,maxi,10),
     ylab="Fn(x) & Fo", xlab="Datos")
xdist=seq(mini,maxi,length=25)
nordist=pnorm(xdist, mean = mumv, sd =sigmamv)
lines(xdist,nordist,type="l",col=2)

#------  8.4  Mi gráfica de probabilidad     ------------------------
X11()
#vector de cuantiles uniformes (0,1) teóricos
enes=(1:N)
alfitas=(1:N)/(N+1)
#se transforman los datos con la distribución normal estimada.
# si el modelo normal es razonable deben ser uniformes (0,1)
uis=pnorm(datsort,mumv,sigmamv)
## Joint 95% confidence band for N prediction intervals of uniform quantiles:
tau=1-(.05/N)
tau1=(1-tau)/2
tau2=(1+tau)/2
aenes=N+1-enes
Ban1=qbeta(tau1,enes,aenes);
Ban2=qbeta(tau2,enes,aenes);
plot(alfitas,uis,pch=19,cex=.5,ylab="F(xi;mumv,sigmamv)",
     xlab="Cuantiles uniformes (0,1)", xlim=c(0,1),ylim=c(0,1),
     main="Gráfica de probabilidad normal",pty="m")
 lines(alfitas,uis,lty=1,col=1)
 lines(c(0,1),c(0,1),lty=1,col=2)
 lines(alfitas,Ban1,lty=2,col=2)
 lines(alfitas,Ban2,lty=2,col=2)
 
 
 
 ##  ............................................................................
 ##     Paso 9 Akaike Normal,                                               #### 
 #      Criterio Información de Akaike       
 # se evalúa la densidad normal estimada en los datos para
 # calcular la densidad normal conjunta de la muestra,DNC
 
 norvec=dnorm(DAT,mean=mumv,sd=sigmamv)
 DNC=prod(norvec)
 AICnorm=-2*log(DNC) + (2*2)
 ## Se calcula promedio distancias verticales (AVD en inglés) entre
 ## Puntos y la recta identidad: sum(abs(ui-E(ui)))/N
 AVD=sum(abs(uis-alfitas))/N
 
 ##  ............................................................................
 ##     Paso 10 Partición pantallas y despliegue gráficas               ####
## Se muestra cómo partir en secciones una pantalla con R:
#se usa split.screen para partir pantalla y screen para seleccionar
#en la que se quiera trabajar. 
#Primero se abre pantalla limpia:
X11()
#Se parte en tres porciones la pantalla
#primero partir toda la pantalla en dos renglones y dos columna
split.screen(c(2,2)) 
#Si se quiere subidividir más una sección se diría:
##split.screen(c(1,2),screen=1) #R da los nombres de estas nuevas
# R va numerando pantallas de manera consecutiva y las nombra así.
#-------------------------------------------- --
#se selecciona ahora la pantalla de arriba a la izquierda:
screen(1)
hresul=hist(DAT,freq=TRUE,xlab="Datos",ylab="Frecuencias",
            main="Histograma datos normales")
#ptos para graficar densidad normal encimada y encimar al histograma:
puntos=seq(min(hresul$breaks),max(hresul$breaks),length=25)
#se evalúan y reescalan densidades multiplicando por el número de 
# datos n y ancho barra hh: 
hh=diff(hresul$mids[1:2]) #ancho barra del histograma
nordensheight=dnorm(puntos, mean = mumv, sd =sigmamv)*N*hh
lines(puntos,nordensheight,type="l",col=2)
#------------------------------------------- ---
#ahora va la pantalla de arriba a la derecha:
screen(2)
mini=floor(datsort[1])
maxi=ceiling(datsort[N])
# type="s" grafica escaleritas:
plot(datsort,(1:N)/(N+1),type="s",ylim=c(0,1),
     main="Fn(x)vs dist. estimada", xaxp=c(mini,maxi,10),
     ylab="Fn(x) & Fo", xlab="Datos")
xdist=seq(mini,maxi,length=25)
nordist=pnorm(xdist, mean = mumv, sd =sigmamv)
lines(xdist,nordist,type="l",col=2)
#------------------------------------------- ---
#ahora va la pantalla de abajo a la izquierda:
screen(3)
x1=mumv-3*sigmamv
x2=mumv+3*sigmamv
densitop=1/(sqrt(2*pi)*sigmamv)
xs=seq(x1,x2, length.out=200)
#Lo sig. es equivalente pero más largo:
#ynor=sapply(X=xs,FUN=dnorm,mumv,sigmamv)
ynor=dnorm(xs,mumv,sigmamv)
plot(xs,ynor,ylab="f(x;mu,sigma)",yaxs="i",
     xlab="x", xlim=c(x1,x2),ylim=c(0,densitop),
     main="Densidad normal estimada y datos",type="l")
points(DAT,rep(0,N),type="p",col="6",pch=19,cex=1)
#------------------------------------------- ---
#ahora va la pantalla de abajo a la derecha:
screen(4)
#vector de cuantiles uniformes (0,1) teóricos
enes=(1:N)
alfitas=(1:N)/(N+1)
#se transforman los datos con la distribución normal estimada.
# si el modelo normal es razonable deben ser uniformes (0,1)
uis=pnorm(datsort,mumv,sigmamv)
## Joint 95% confidence band for N prediction intervals of uniform quantiles:
tau=1-(.05/N)
tau1=(1-tau)/2
tau2=(1+tau)/2
aenes=N+1-enes
Ban1=qbeta(tau1,enes,aenes);
Ban2=qbeta(tau2,enes,aenes);
plot(alfitas,uis,pch=19,cex=.5,ylab="F(xi;mumv,sigmamv)",
     xlab="Cuantiles uniformes (0,1)", xlim=c(0,1),ylim=c(0,1),
     main="Gráfica Probab normal",pty="m")
lines(c(0,1),c(0,1),lty=2,col=2)
lines(alfitas,Ban1,lty=2,col=2)
lines(alfitas,Ban2,lty=2,col=2)
#-------------------------------------- ---
## ___________________________ ####
## Comparación con Lognormal para los mismos datos     ####

#Se calculan estimadores de muLN y sigLN:
Y=log(datsort)
t1LN=sum(Y)
t2LN=sum(Y^2)
HLN=t2LN-(t1LN^2)/N
muLN=t1LN/N
sigLN=sqrt(HLN/N)
uisLN=pnorm(Y,muLN,sigLN)
#Promedio distancias verticales modelo Lognormal:
AVDLN=sum(abs(uisLN-alfitas))/N
"Van Prom Dist Verticales Normal y Lognormal:"
"Es mejor el modelo con el valor más chico."
c(AVD,AVDLN)
"Van AICs Normal y  Lognormal"
#se calcula aquí el AIC Lognormal con los datos originales:
LNvec=dlnorm(DAT,mean=muLN,sd=sigLN)
DNCLN=prod(LNvec) #esta es la densidad conjunta estimada Lognormal
AICLogN=-2*log(DNCLN) + (2*2)
c(AICnorm,AICLogN)

# Gráfica de probabilidad Lognormal:
X11()
plot(alfitas,uisLN,pch=19,cex=.5,ylab="F(xi;mumv,sigmamv)",
     xlab="Cuantiles uniformes (0,1)", xlim=c(0,1),ylim=c(0,1),
     main="Gráfica de probabilidad Log Normal",pty="m")
lines(alfitas,uisLN,lty=1,col=1)
lines(c(0,1),c(0,1),lty=1,col=2)
lines(alfitas,Ban1,lty=2,col=2)
lines(alfitas,Ban2,lty=2,col=2)



## ___________________________ ####

##  ............................................................................
##       Generando archivo resultados:                 ####


sink(file="NormalOutput.txt")
"LOs EMV de mu y sigma son:"
mumv
sigmamv

"Van intervalos verosimilitud perfil para MU (primeras dos columnas)"
"luego niveles verosimilitud (col 3) y confianza (col 4):"
Intermus

"--------------------"
"Van intervalos propuestos para SIGMA (primeras dos columnas)
luego niveles verosimilitud (col 3) y confianza (col 4):"
Intersigs


"--------------------"
"Van intervalos perfil exactos para SIGMA (primeras dos columnas)
luego niveles verosimilitud (col 3) y confianza (col 4):"
IntSigsNum


"--------------------"
"Se estimó también el modelo Lognormal para estos datos."
"Los AIC normal y lognormal fueron: "
c(AICnorm,AICLogN)

"Van Prom Dist Verticales Normal y Lognormal:"
"Es mejor el modelo con el valor más chico."
c(AVD,AVDLN)

sink(file=NULL)

#   ____________________________________________________________________________
#                     FIN DEL PROGRAMA                                  ####
#   ________________________________________________________________________ __






