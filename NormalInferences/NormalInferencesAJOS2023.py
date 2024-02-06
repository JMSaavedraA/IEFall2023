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
#  11) Se comparan las relaciones entre niveles de verosimilitud y n
#      para los parámetros normales, los Gamma y los Weibull.

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import norm, beta, lognorm
import AuxNormalInferences as ANI
import pickle
import matplotlib.pyplot as plt

## Cargar datos de pickle

# Ejemplo que queremos cargar
file_path = 'P1g.pickle'

# Cargar los datos
with open(file_path, 'rb') as file:
    DAT = pickle.load(file)

## Histograma de los datos

plt.hist(DAT, bins='auto', alpha=0.7, rwidth=0.85)
plt.xlabel('Datos')
plt.ylabel('Frecuencias')
plt.title('Histograma datos')
plt.grid(True)
plt.show()


##  ............................................................................
##       Paso 2 Emvs analíticos                                        ####

N = len(DAT)
datsort = np.sort(DAT)
#van estadísticas suficientes
t1 = np.sum(DAT)
t2 = np.sum(DAT ** 2)
#suma de cuadrados
H = t2 - (t1 ** 2) / N
#Se calculan los emv de mu y sigma^2 analíticamente:
mumv = t1 / N
sig2mv = H / N
sigmamv = np.sqrt(sig2mv)


##  ............................................................................
##         Paso 3 Emvs numéricos                                          ####

#   Van los emv con minimize con el método de gradiente conjugado, que no necesita definir el gradiente ni la hessiana de la logverosimilitud. 
#   La salida de esta función es una lista:
#   El primer componente es el vector de emvs, donde alcanza el mínimo de la LV
#   El segundo componente es el valor de la función en el emv
#   El tercer componente es el número de iteraciones que toma para optimizar.



#   Valores iniciales de mu y sigma
vecini = [mumv, sigmamv]

#   Optimización con Gradiente Conjugado
maxversali = minimize(ANI.lvnor, vecini, args=(t1, t2, N), method='CG')

#   EMVs numéricos (deben ser iguales a los analíticos):
emvs = maxversali.x
mumvnum = emvs[0]
sigmvnum = emvs[1]
#   Valor de la logver en este máximo:
LVmax = -maxversali.fun
#   Número de iteraciones del proc numérico:
iteras = maxversali.nfev

#   Despliegue de resultados en pantalla:
print("Van primero los emv analíticos de mu y sigma:")
print(mumv, sigmamv)
print("Ahora van los emv numéricos con Optim:")
print(mumvnum, sigmvnum)
print("Ahora va la logver en el emv y las iteraciones (llamados a logver):")
print(LVmax, iteras)

##        Paso 4 Rp(mu) e intervalos                                        ####

##  Se calcula y grafica Rp(mu) (analíticamente calculada) con intervalo
#   de verosimilitud de nivel c. Va la rejilla donde se grafica:
MM = 150
mu1 = mumv - 1.5 * sigmamv
mu2 = mumv + 1.5 * sigmamv
mus = np.linspace(mu1, mu2, num=MM)

#   Se calcula la perfil relativa en el vector de mus

rpmuss = np.array([ANI.rpmu(x,t1, t2, N) for x in mus])

#   Se calculan los niveles de verosimilitud de intervalos con 90, 95 y 
#   99% de confianza

cesmu = ANI.clevels(N)['cmu']

#   cesmu sin nombres:

CesMuS = np.array(list(cesmu.values()))  # Assuming sinnombre(cesmu) just returns cesmu

aux = sigmamv * np.sqrt((CesMuS ** (-2 / N)) - 1)
Intmus1 = mumv - aux
Intmus2 = mumv + aux
Intermu = np.column_stack((Intmus1, Intmus2))

#   Para desplegar resultados al final:
confianzas = ["90%", "95%", "99%"]
Intermus = np.column_stack((Intermu, CesMuS, confianzas))

#   Graficando Rp(mu)
plt.figure()
plt.plot(mus, np.exp(rpmuss), color='black', linewidth=1)
plt.title("Rp(mu) analítica")
plt.xlabel(r'$\mu$')
plt.ylabel(r'$R_p(\mu)$')
plt.ylim(0, 1)

cvec = np.array([x * np.ones((2)) for x in CesMuS])

for i in range(3):
    plt.plot(Intermu[i, :], cvec[i, :], linestyle='--', color='red')

plt.show()

##       Paso 5 Rp(mu) numérica                                             ####

# ------ Se calcula y grafica la perfil de mu pero ahora numéricamente
##  Para compararla y ratificar que son indistinguibles
#   Rejilla regular para mus donde se graficará la perfil:

lvpmu = np.ones(MM)  # logverosim perfil para la mu fija en mus
sigrestmv = np.zeros(MM)  # los emv restringidos de sigma para esas mu fijas

for j in range(MM):
    interv = [0.001, sigmamv * 15]  # int de búsqueda para sigma
    MUFIJO = mus[j]
    objective = lambda sigma: -ANI.lpmu(sigma, t1, t2, N, MUFIJO)
    result = minimize_scalar(objective, bounds=np.sort(interv), method='bounded')
    lvpmu[j] = -result.fun
    sigrestmv[j] = result.x

#   Se calcula verosimilitud relativa
Rpmunum = np.exp(lvpmu - ANI.lpmu(sigmamv, t1, t2, N, mumv))

# se grafica la nueva perfil de mu numérica sobre la exacta para compararlas
plt.figure()
plt.plot(mus, Rpmunum, color='black', linewidth=1, label='Rp(mu) numérica')
plt.plot(mus, np.exp(rpmuss), linestyle='--', color='red', label='Rp(mu) analítica')
plt.title("Rp(mu) numérica negro y analítica rojo (guiones)")
plt.xlabel(r'$\mu$')
plt.ylabel(r'$R_p(\mu)$')
plt.ylim(0, 1)
plt.legend()
plt.show()


##      Paso 6  Regiones Verosimilitud (mu,sigma)                                                   ####
#   Contornos de la verosimilitud relativa 
#   Se crea matriz aprovechando rejilla de mus y ahora
#   se hace para sigmas, donde se guardara la logverosimilitud de mu,sigma


NN = 50
mu1 = mumv - 1.5 * sigmamv
mu2 = mumv + 1.5 * sigmamv
mucont = np.linspace(mu1, mu2, num=NN)
sig1 = sigmamv / 2
sig2 = sigmamv * 2
sigis = np.linspace(sig1, sig2, num=NN)
lvmatriz = np.zeros((NN, NN))

# Calculate l(mu, sigma) on this grid
for i in range(NN):
    for j in range(NN):
        veci2 = np.array([mucont[i], sigis[j]])
        lvmatriz[i, j] = ANI.lvnor(veci2, t1, t2, N)

# Calculate relative log-likelihood on the grid
LVmax = ANI.lvnor(np.array([mumv, sigmamv]), t1, t2, N)
Rvmatriz = np.exp(lvmatriz - LVmax)

# Plot likelihood and confidence regions
plt.figure()
cesmusig = ANI.clevels(N)['cmusig']
inverted_cesmusig = {value: key for key, value in cesmusig.items()}
contour = plt.contour(mucont, sigis, Rvmatriz, levels= np.sort(np.array(list(cesmusig.values()))))
plt.clabel(contour, contour.levels, inline=True, fmt=inverted_cesmusig, fontsize=10)
plt.xlim(mu1, mu2)
plt.ylim(sig1, sig2)
plt.xlabel(r'$\mu$')
plt.ylabel(r'$\sigma$')
plt.scatter(mumv, sigmamv, color='red', marker='o', label='EMV')
plt.plot(mus, sigrestmv, color='blue', marker='o', label='Perfil de mu')
plt.legend()
plt.show()


##       Paso 7 Rp(sigma) e intervalos                               ####
##   Intervalos para Sigma de verosimilitud perfil y confianza

# Assuming the required functions and variables are defined

# Calculations similar to R code
cessig = np.array(list(ANI.clevels(N)['csig'].values()))
aux = np.sqrt(-4 * np.log(cessig) / N) / 3
sig1 = sigmamv * ((1 + aux) ** (-3 / 2))
sig2 = sigmamv * ((1 - aux) ** (-3 / 2))
Intersig = np.column_stack((sig1, sig2))
confianzas = ["90%", "95%", "99%"]
Intersigs = np.column_stack((sig1,sig2, cessig, confianzas))
S1 = sigmamv / 2
S2 = sigmamv * 2
SIGS = np.linspace(S1, S2, num=200)
rpsig = N * np.log(sigmamv / SIGS) - (N / 2) * ((sigmamv / SIGS) ** 2) + (N / 2)

# Calculation of interval profile likelihood numerically
InumSig90 = ANI.IntSig(sigmamv, N, 0.90)
InumSig95 = ANI.IntSig(sigmamv, N, 0.95)
InumSig99 = ANI.IntSig(sigmamv, N, 0.99)
IntSigsNum = np.vstack((InumSig90, InumSig95, InumSig99))

# Plotting
plt.figure()
plt.plot(SIGS, np.exp(rpsig), color='black', linewidth=1)
plt.title("Rp(sigma)")
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$R_p(\sigma)$')
plt.ylim(0, 1)

cvec = np.array([x * np.ones((2)) for x in cessig])
for i in range(3):
    plt.plot(Intersig[i, :], cvec[i, :], linestyle='--', color='red')

plt.show()


##    Paso 8   Validación Modelo Normal                             ####

##     8.1) Histograma vs densidad normal estimada
#      8.2) Densidad estimada y datos 
#      8.3) Distribución empirica vs estimada
#      8.4) Mi gráfica de probabilidad


#---- 8.1 Histograma y densidad estimada:--------------

plt.figure()
hresul = plt.hist(DAT, bins='auto', density=False, color='skyblue', alpha=0.7)
puntos = np.linspace(min(hresul[1]), max(hresul[1]), num=25)
hh = np.diff(hresul[1][0:2])
nordensheight = norm.pdf(puntos, loc=mumv, scale=sigmamv) * N * hh
plt.plot(puntos, nordensheight, color='red', linestyle='-', linewidth=2)
plt.xlabel('Datos')
plt.ylabel('Frecuencias')
plt.title('Histograma con densidad normal estimada')
plt.show()

#---- 8.2  Gráfica Densidad Estimada y datos ------------------
plt.figure()
x1 = mumv - 3 * sigmamv
x2 = mumv + 3 * sigmamv
densitop = 1/(np.sqrt(2*np.pi)*sigmamv)
xs = np.linspace(x1, x2, num=200)
ynor = norm.pdf(xs, loc=mumv, scale=sigmamv)
plt.plot(xs, ynor, color='blue', linewidth=2, label='Densidad Normal Estimada')
plt.scatter(DAT, np.zeros_like(DAT), color='green', marker='o', label='Datos')
plt.xlabel('x')
plt.ylabel('f(x; mu, sigma)')
plt.ylim(0, densitop)
plt.title('Densidad normal estimada y datos')
plt.legend()
plt.show()

#-----  8.3 Distribución empirica vs estimada: --------------

plt.figure()
mini = datsort[0]
maxi = datsort[-1]
xdist = np.linspace(mini, maxi, num=25)
nordist = norm.cdf(xdist, loc=mumv, scale=sigmamv)
plt.step(datsort, np.arange(1, len(datsort) + 1) / len(datsort), label='Fn(x)',where='post')
plt.plot(xdist, nordist, color='red', label='Distribución Estimada')
plt.xlabel('Datos')
plt.ylabel('Fn(x) & Fo')
plt.title('Distribución empírica vs. estimada')
plt.legend()
plt.show()

#------  8.4  Mi gráfica de probabilidad     ------------------------

plt.figure()
#vector de cuantiles uniformes (0,1) teóricos
enes = np.arange(1, len(DAT) + 1)
alfitas = enes / (len(DAT) + 1)
#   se transforman los datos con la distribución normal estimada.
#   si el modelo normal es razonable deben ser uniformes (0,1)
uis = norm.cdf(datsort, loc=mumv, scale=sigmamv)
## Joint 95% confidence band for N prediction intervals of uniform quantiles:
tau = 1 - (0.05 / len(DAT))
tau1 = (1 - tau) / 2
tau2 = (1 + tau) / 2
aenes = len(DAT) + 1 - enes

# Calculate beta distribution quantiles
Ban1 = beta.ppf(tau1, enes, aenes)
Ban2 = beta.ppf(tau2, enes, aenes)

plt.scatter(alfitas, uis, color='blue', s=5, label='Datos transformados')
plt.plot([0, 1], [0, 1], color='red', linestyle='--', label='Identidad')
plt.plot(alfitas, Ban1, color='green', linestyle='--', label='Banda de confianza')
plt.plot(alfitas, Ban2, color='green', linestyle='--')
plt.xlabel('Cuantiles uniformes (0,1)')
plt.ylabel('F(xi; mumv, sigmamv)')
plt.xlim(0,1)
plt.ylim(0,1)
plt.title('Gráfica de probabilidad normal con Beta Quantiles')
plt.legend()
plt.show()

 ##     Paso 9 Akaike Normal,                                               #### 
 #      Criterio Información de Akaike       
 # se evalúa la densidad normal estimada en los datos para
 # calcular la densidad normal conjunta de la muestra,DNC
 
norvec = norm.pdf(DAT, loc=mumv, scale=sigmamv)
DNC = np.prod(norvec)

# Calculate AIC for normal distribution
AICnorm = -2 * np.log(DNC) + (2 * 2)  # Assuming 2 parameters (mean and standard deviation)

## Se calcula promedio distancias verticales (AVD en inglés) entre
 ## Puntos y la recta identidad: sum(abs(ui-E(ui)))/N
AVD = np.sum(np.abs(uis - alfitas)) / N


##  Paso 10. Curvas niveles verosim: c(n) vs n                                                ####
#   CONCLUSIONES. LA CURVA DE SIGMA NORMAL PARECE
#   FUNCIONAR MUY BIEN PARA TODOS LOS PARAMETROS PARA
#   LOS NIVELES 90, 95 Y 99%
#   PARA OTROS NIVELES DE CONFIANZA ES FACIL DE USAR LA RELACION DE
#   LA MEDIA NORMAL DE C Y EL CUANTIL T DE STUDENT

enes = np.arange(10, 101, 1)

#   Nivel 0.90

csignor90 = 0.2585 - 1 / (0.775 + (1.641 * enes))
cwei90 = 0.2585 - 1 / (0.6818 + (1.5635 * enes))
cgam90 = 0.2585 - 1 / (0.5171 + (1.6667 * enes))
cmunor90 = 0.2585 - 1 / (-0.4743 + (1.9029 * enes))

plt.figure()
plt.plot(enes, csignor90, label='Normal', color='black')
plt.axhline(y=0.2585, linestyle='--', color='blue')
plt.plot(enes, cwei90, label='Weibull', linestyle='--', color='blue')
plt.plot(enes, cgam90, label='Gamma', linestyle='--', color='red')
plt.plot(enes, cmunor90, label='Nmugreen', color='green')
plt.title('90% LikeLevels vs n')
plt.xlabel('n')
plt.ylabel('c(n)')
plt.legend()
plt.show()

#   Nivel 0.95

csignor95 = 0.1465 - 1 / (2.757 + (2.082 * enes))
cwei95 = 0.1465 - 1 / (2.3025 + (1.9278 * enes))
cgam95 = 0.1465 - 1 / (2.029 + (2.084 * enes))
cmunor95 = 0.1465 - 1 / (0.7538 + (2.3544 * enes))

plt.figure()
plt.plot(enes, csignor95, label='Norsig', color='black')
plt.axhline(y=0.1465, linestyle='--', color='blue')
plt.plot(enes, cwei95, label='Weibull', linestyle='--', color='blue')
plt.plot(enes, cgam95, label='Gamma', linestyle='--', color='red')
plt.plot(enes, cmunor95, label='Nmugreen', color='green')
plt.title('95% LikeLevels vs n')
plt.xlabel('n')
plt.ylabel('c(n)')
plt.legend()
plt.show()

#   Nivel 0.99

csignor99 = 0.0362 - 1 / (18.294 + (5.229 * enes))
cwei99 = 0.0362 - 1 / (19.8329 + (3.8032 * enes))
cgam99 = 0.0362 - 1 / (14.609 + (4.886 * enes))
cmunor99 = 0.0362 - 1 / (8.753 + (5.551 * enes))

plt.figure()
plt.plot(enes, csignor99, label='Norsig', color='black')
plt.axhline(y=0.0362, linestyle='--', color='blue')
plt.plot(enes, cwei99, label='Weibull', linestyle='--', color='blue')
plt.plot(enes, cgam99, label='Gamma', linestyle='--', color='red')
plt.plot(enes, cmunor99, label='Nmugreen', color='green')
plt.title('99% LikeLevels vs n')
plt.xlabel('n')
plt.ylabel('c(n)')
plt.legend()
plt.show()

##      Paso 12 Curvas c(n) para regiones (mu,sigma)              ####
#Se comparan curvas del 95% confianza para Weibull y la normal

# Calculations and plotting for 90% confidence level
cregnor90 = 0.1 - 1 / (1.11 + (4.655 * enes))
cregwei90 = 0.1 - 1 / (0.9774 + (4.4508 * enes))
cregam90 = 0.1 - 1 / (1.06 + (4.617 * enes))

plt.figure()
plt.plot(enes, cregnor90, label='Normal', color='black')
plt.axhline(y=0.1, linestyle='--', color='blue')
plt.plot(enes, cregwei90, label='Weibull', linestyle='--', color='blue')
plt.plot(enes, cregam90, label='Gamma', linestyle='--', color='red')
plt.title('Likelihood Regions 90%')
plt.xlabel('n')
plt.ylabel('c(n)')
plt.legend()
plt.ylim(0.075, 0.105)
plt.show()

# Calculations and plotting for 95% confidence level
cregnor95 = 0.05 - 1 / (1.61 + (7.27 * enes))
cregwei95 = 0.05 - 1 / (3.671 + (6.9103 * enes))
cregam95 = 0.05 - 1 / (4.932 + (7.030 * enes))

plt.figure()
plt.plot(enes, cregnor95, label='Normal', color='black')
plt.axhline(y=0.05, linestyle='--', color='blue')
plt.plot(enes, cregwei95, label='Weibull', linestyle='--', color='blue')
plt.plot(enes, cregam95, label='Gamma', linestyle='--', color='red')
plt.title('Likelihood Regions 95%')
plt.xlabel('n')
plt.ylabel('c(n)')
plt.legend()
plt.ylim(0.035, 0.052)
plt.show()

# Calculations and plotting for 99% confidence level
cregnor99 = 0.01 - 1 / (19.51 + (23.59 * enes))
cregwei99 = 0.01 - 1 / (17.4989 + (25.1944 * enes))
cregam99 = 0.01 - 1 / (51.82 + (21.66 * enes))

plt.figure()
plt.plot(enes, cregnor99, label='Normal', color='black')
plt.axhline(y=0.01, linestyle='--', color='blue')
plt.plot(enes, cregwei99, label='Weibull', linestyle='--', color='blue')
plt.plot(enes, cregam99, label='Gamma', linestyle='--', color='red')
plt.title('Likelihood Regions 99%')
plt.xlabel('n')
plt.ylabel('c(n)')
plt.legend()
plt.ylim(0.005, 0.011)
plt.show()


#   Se calculan estimadores de muLN y sigLN para la distribución lognormal:
Y = np.log(datsort)
t1LN = np.sum(Y)
t2LN = np.sum(Y**2)
HLN = t2LN - (t1LN**2) / N
muLN = np.mean(Y)
sigLN = np.std(Y)
uisLN = norm.cdf(Y, loc=muLN, scale=sigLN)

#   Promedio distancias verticales modelo Lognormal:
AVDLN = np.sum(np.abs(uisLN - alfitas)) / len(Y)


print("Es mejor el modelo con el valor más chico.")
print([AVD,AVDLN])

LNvec = lognorm.pdf(DAT, s=sigLN, scale=np.exp(muLN))
DNCLN = np.prod(LNvec)
# AIC for Lognormal and Normal models
AICLogN = -2 * np.log(DNCLN) + (2 * 2)

print(f"Promedio de Distancias Verticales. Normal y Lognormal: {AVD} vs {AVDLN}")
print(f"AIC. Normal vs Lognormal: {AICnorm} vs {AICLogN}")

plt.figure()
plt.plot(alfitas, uisLN, marker='o', linestyle='', label='Log Normal')
plt.plot(alfitas, alfitas, linestyle='--', label='Identity')
plt.plot(alfitas, Ban1 * np.ones(N), linestyle='--', label='Ban1')
plt.plot(alfitas, Ban2 * np.ones(N), linestyle='--', label='Ban2')
plt.xlabel('Cuantiles uniformes (0,1)')
plt.ylabel('F(xi;mu,sigma)')
plt.title('Gráfica de probabilidad Log Normal')
plt.legend()
plt.grid(True)
plt.show()


