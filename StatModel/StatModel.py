
# Inferencias para 7 distribuciones estadísticas y selección
# del mejor modelo para datos exactos (con posibles repeticiones).
# Aparte se realizan inferencias para mismas distribuciones pero con 
# datos censurados por intervalo en StatModelCEN.R

# Revisado por Eloísa Díaz Francés en octubre de 2023

#               PROGRAMA PRINCIPAL                                    ####
#  ............................. ........

import AuxStatModel as ASM
import pickle
import math
import numpy as np
from scipy.optimize import root_scalar
from scipy.stats import beta, gamma, norm, lognorm, invgauss
import matplotlib.pyplot as plt



# Ejemplo que queremos cargar
file_path = 'Ej1.pickle'

# Cargar los datos
with open(file_path, 'rb') as f:
    X = pickle.load(f)

xsort = np.sort(X) # datos ordenados
n = len(X)

# Épsilon mitad del tamaño del intervalo de censura. En caso de desconocerlo, estimar:
# difimin = 1.0
# Diferencia entre valores sin repeticiones
Xuni = np.unique(xsort)
difis = Xuni[1:] - Xuni[:-1]
# Diferencia minima como tamaño del intervalo de censura
difimin = np.min(difis)

epsi = difimin / 2

# Intervalos en XMAT
XMAT = np.column_stack((xsort - epsi, xsort + epsi))

# ____________________________________ ####
#    Estadísticas suficientes:                                   ####
t1 = np.sum(X)  # Para distribuciones Exp, Normal, Lognormal, Gamma e Inv Gaussiana
t2 = np.sum(X ** 2) # Para la Normal y Lognormal
t3 = np.sum(np.log(X))  # Para la Gamma
t4 = np.sum(1 / X)  # Para la Gaussiana Inversa
##  ............................................................................
# EMVs Gamma                                                           ####
# Estimador momentos forma Gamma:
alphamom = (t1**2) / ((n * t2) - (t1**2))
# EMV de media Gamma, Exponencial y Normal:
mumle = t1 / n
# EMV de parámetros de forma Gamma, alphamle,  y escala, betagamle:
root_func = lambda alph: ASM.Lgampsi(alph, t1, t3, n)
result = root_scalar(root_func, bracket=[alphamom / 10, alphamom * 10], xtol=0.000001)
alphamle = result.root
betagamle = mumle / alphamle

##  ............................................................................
# EMVs  Weibull                                             ####
# Log data:
Z = np.log(X)
tz = np.sum(Z)
zsort = np.sort(Z)
zbar = np.mean(Z)
#  Estimador de momentos Gumbel para forma beta de la Weibull:
betaWeimom = np.pi * np.sqrt(len(Z)) / np.sqrt(6 * np.sum((Z - zbar)**2))
# EMV de forma beta aquí:
result = root_scalar(ASM.Weibetmle, bracket=[betaWeimom / 10, betaWeimom * 10], args=(Z, tz, len(Z)), xtol=0.000001)
betaWei = result.root
# Se calcula el emv de a, localización de los logdatos Gumbel.
# Notar que la escala WEibull sigma=exp(a).
K1 = ASM.K1bet(betaWei, Z)
K2 = ASM.K2bet(betaWei, Z)
K3 = ASM.K3bet(betaWei, Z)
amle = np.log(K1 / len(Z)) / betaWei
##  ............................................................................
#  Emvs Gumbel para datos X                                             ####
# Ahora Z=exp(X) son Weibull de parámetros amlegum,betaGum.

tzG = np.sum(xsort)
# Estimador de momentos Gumbel para beta Weibull de Z=exp(X):
zbarG = np.mean(xsort)
betaGmom = np.pi * np.sqrt(len(xsort)) / np.sqrt(6 * np.sum((xsort - zbarG)**2))
# Emv de beta:
result = root_scalar(ASM.Weibetmle, bracket=[betaGmom / 10, betaGmom * 10], args=(xsort, tzG, len(xsort)), xtol=0.000001)
betaGum = result.root
# Emv de Gumbel de localización:
K1G = ASM.K1bet(betaGum, xsort)
K2G = ASM.K2bet(betaGum, xsort)
K3G = ASM.K3bet(betaGum, xsort)
amleGum = np.log(K1G / n) / betaGum
tetamle = amleGum
deltamle = 1 / betaGum

##  ............................................................................
# EMVs Gaussiana Inversa:                                       ####
# La media mumle es el emv de mu. Este es el emv del parámetro de forma lambda:
lambmle = 1 / ((t4 / n) - (n / t1))

# EMVs Normales            ####
H = t2 - ((t1 ** 2) / n)
sigmle = np.sqrt(H / n)

# ____________________________________ ####
#  Validación de modelos                                        ####
##  ............................................................................
#  AICs:----------                                                     ####
##  ............................................................................
#   AIC Gamma                                                              ####
LnfgamY = -n * math.lgamma(alphamle) - (n * alphamle * np.log(betagamle)) + (alphamle - 1) * t3 - (t1 / betagamle)
AICgam = -2 * LnfgamY + 4

# AIC Exponencial
AICexp = 2 * n * np.log(t1 / n) + 2 * n + 2

# AIC Normal
AICnor = n * np.log(2 * np.pi * H / n) + n + 4

# AIC LogNormal
Z = np.log(xsort)
t1LN = np.sum(Z)
t2LN = np.sum(Z**2)
HLN = t2LN - (t1LN**2) / n

# Lognormal mles:
muLN = t1LN / n
sigLN = np.sqrt(HLN / n)
AICLognor = 2 * t1LN + n * np.log(2 * np.pi * HLN / n) + n + 4

#   AIC Weibull para datos X o Gumbel para logdatos Y=lnX                             ####
# La logdensidad estiamda Weibull en (amle,betaWei) es:
# LnfWeiX=lvabetaWei(c(amle,betaWei),Z,tz,n) - tz
# AICWei2= -2*LnfWeiX + 4
# se calcula directo aquí el AIC Weibull:
AICWei = 2 * t3 * (1 - betaWei) - 2 * n * np.log(betaWei) + 2 * n * amle * betaWei + 2 * np.exp(-amle * betaWei) * K1 + 4

#   AIC Gumbel para datos X:     ####
# amleGum y betaGum calculados para datos X que son Gumbel
AICGum = 2 * n * np.log(deltamle) - 2 * (t1 - (n * tetamle)) / deltamle + 2 * np.exp(-tetamle / deltamle) * np.sum(np.exp(xsort / deltamle)) + 4

#   AIC Gaussiana Inversa para datos X    ####
AICGauInv = n * np.log(2 * np.pi) - n * np.log(lambmle) + 3 * t3 - (2 * n * lambmle) / mumle + lambmle * ((t1 / (mumle**2)) + t4) + 4

print("Criterio de Info.de Akaike para datos exactos y distribuciones estimadas:")
AICs = [AICexp, AICgam, AICnor, AICLognor, AICWei, AICGum, AICGauInv]
distribution_names = ["Exponencial", "Gamma", "Normal", "Lognormal", "Weibull", "Gumbel", "Gaussiana Inversa"]

aic_dict = dict(zip(distribution_names, AICs))
sorted_aic_dict = dict(sorted(aic_dict.items(), key=lambda item: item[1], reverse=True))

max_length = 2 + max(len(distribution) for distribution in sorted_aic_dict.keys())

# Imprimir valores de AIC:
for distribution, aic_value in sorted_aic_dict.items():
    print(f"{distribution.ljust(max_length)}: {aic_value}")

##  ............................................................................
# WAVD Weighted AVerage of Vertical Distances in Probability Plots  ####
# (Estadística propuesta por Eloísa) 
# Se transforman datos ordenados con las 7 distribuciones estimadas: 
uexp = 1 - np.exp(-xsort / mumle)
ugam = gamma.cdf(xsort, a=alphamle, scale=betagamle)
unor = norm.cdf(xsort, loc=mumle, scale=sigmle)
ulognor = lognorm.cdf(xsort, s=sigLN, scale=np.exp(muLN))
#Weibull distrib para datos X es igual que Gumbel para logdatos Z
ugum = 1 - np.exp(-np.exp(betaWei * (zsort - amle)))
#GUmbel distrib directa para datos ordenados X en xsort:
ugumZG = 1 - np.exp(-np.exp(betaGum * (xsort - amleGum)))
#Gaussiana Inversa:
aux1 = np.sqrt(lambmle / xsort) * ((xsort / mumle) - 1)
aux2 = -np.sqrt(lambmle / xsort) * ((xsort / mumle) + 1)
#esta es la distribución GI con base en la distrib. normal:
uGauInv = norm.cdf(aux1) + np.exp(2 * lambmle / mumle) * norm.cdf(aux2)



# WAVDs for each distribution
WAVDExp = ASM.WAVD(uexp)
WAVDGam = ASM.WAVD(ugam)
WAVDNor = ASM.WAVD(unor)
WAVDLN = ASM.WAVD(ulognor)
WAVDWei = ASM.WAVD(ugum)
WAVDGum = ASM.WAVD(ugumZG)
WAVDGI = ASM.WAVD(uGauInv)

# Creating a dictionary with keys as distribution names and values as corresponding WAVD values
WAVDs = [WAVDExp, WAVDGam, WAVDNor, WAVDLN, WAVDWei, WAVDGum, WAVDGI]
wavds = dict(zip(distribution_names, WAVDs))

sorted_wavds = dict(sorted(wavds.items(), key=lambda item: item[1]))
print("WAVDs: Promedios ponderados de distancias verticales:")
# Imprimir valores de WAVDs:
for distribution, wavd_value in sorted_wavds.items():
    print(f"{distribution.ljust(max_length)}: {wavd_value}")

#---------------------------- --
#-------------------- --
#   Anderson Darling para 7 distribucioness             ####
#   Se calcula la estadística AD para las 7 distribs:
ADExp = ASM.AndersonDarling(uexp)
ADGam = ASM.AndersonDarling(ugam)
ADNor = ASM.AndersonDarling(unor)
ADLN = ASM.AndersonDarling(ulognor)
ADWei = ASM.AndersonDarling(ugum)
ADGum = ASM.AndersonDarling(ugumZG)
ADGI = ASM.AndersonDarling(uGauInv)

ADs = [ADExp, ADGam, ADNor, ADLN, ADWei, ADGum, ADGI]
ads = dict(zip(distribution_names, ADs))

sorted_ads = dict(sorted(ads.items(), key=lambda item: item[1]))
print("Anderson Darling para datos exactos con distribuciones")
# Imprimir valores de WAVDs:
for distribution, ad_value in sorted_ads.items():
    print(f"{distribution.ljust(max_length)}: {ad_value}")
#--------------------- --
# Gráficas Probabilidad para 7 distribuciones:    ####
# Ebetas son los cuantiles teóricos de una uniforme (0,1):
# que corresponde a cuantiles empíricos de probabilidades i/(n+1):
enes = np.arange(1, n+1, 1)
aenes = n + 1 - enes
Ebetas = enes / (n + 1)
varian = (enes * (n + 1 - enes)) / (((n + 1)**2) * (n + 2))

# Banda del 95% confianza para N intervalos de predicción simultánea para uks:
alfa = 0.05
tau = 1 - (alfa / n)
tau1 = (1 - tau) / 2
tau2 = (1 + tau) / 2
Ban1 = beta.ppf(tau1, enes, aenes)
Ban2 = beta.ppf(tau2, enes, aenes)

# Probability plots


# Exponential Probability Plot
plt.figure()
plt.plot(Ebetas, uexp, 'o', markersize=5, label='Empirical')
plt.plot([0, 1], [0, 1], 'r--', label='45-degree line')
plt.plot(Ebetas, Ban1, 'r--', label='95% Confidence Band')
plt.plot(Ebetas, Ban2, 'r--')
plt.title('Exponential Probability Plot')
plt.xlabel('Uniform (0,1) Quantiles')
plt.ylabel('F(xi; mumv)')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.legend()
plt.show()

# Gamma Probability Plot
plt.figure()
plt.plot(Ebetas, ugam, 'o', markersize=5, label='Empirical')
plt.plot([0, 1], [0, 1], 'r--', label='45-degree line')
plt.plot(Ebetas, Ban1, 'r--', label='95% Confidence Band')
plt.plot(Ebetas, Ban2, 'r--')
plt.title('Gamma Probability Plot')
plt.xlabel('Uniform (0,1) Quantiles')
plt.ylabel('F(xi; mumv, alfamv)')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.legend()
plt.show()

# Normal Probability Plot
plt.figure()
plt.plot(Ebetas, unor, 'o', markersize=5, label='Empirical')
plt.plot([0, 1], [0, 1], 'r--', label='45-degree line')
plt.plot(Ebetas, Ban1, 'r--', label='95% Confidence Band')
plt.plot(Ebetas, Ban2, 'r--')
plt.title('Normal Probability Plot')
plt.xlabel('Uniform (0,1) Quantiles')
plt.ylabel('F(xi; mumv, sigmamv)')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.legend()
plt.show()

# Lognormal Probability Plot
plt.figure()
plt.plot(Ebetas, ulognor, 'o', markersize=5, label='Empirical')
plt.plot([0, 1], [0, 1], 'r--', label='45-degree line')
plt.plot(Ebetas, Ban1, 'r--', label='95% Confidence Band')
plt.plot(Ebetas, Ban2, 'r--')
plt.title('Lognormal Probability Plot')
plt.xlabel('Uniform (0,1) Quantiles')
plt.ylabel('F(xi; muLNmv, sigmaLNmv)')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.legend()
plt.show()

# Weibull Probability Plot
plt.figure()
plt.plot(Ebetas, ugum, 'o', markersize=5, label='Empirical')
plt.plot([0, 1], [0, 1], 'r--', label='45-degree line')
plt.plot(Ebetas, Ban1, 'r--', label='95% Confidence Band')
plt.plot(Ebetas, Ban2, 'r--')
plt.title('Weibull Probability Plot')
plt.xlabel('Uniform (0,1) Quantiles')
plt.ylabel('F(xi; betaWei, aGum)')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.legend()
plt.show()

# Gumbel Probability Plot
plt.figure()
plt.plot(Ebetas, ugumZG, 'o', markersize=5, label='Empirical')
plt.plot([0, 1], [0, 1], 'r--', label='45-degree line')
plt.plot(Ebetas, Ban1, 'r--', label='95% Confidence Band')
plt.plot(Ebetas, Ban2, 'r--')
plt.title('Gumbel Probability Plot')
plt.xlabel('Uniform (0,1) Quantiles')
plt.ylabel('F(xi; betaGum, aGum)')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.legend()
plt.show()

# Gumbel Probability Plot
plt.figure()
plt.plot(Ebetas, uGauInv, 'o', markersize=5, label='Empirical')
plt.plot([0, 1], [0, 1], 'r--', label='45-degree line')
plt.plot(Ebetas, Ban1, 'r--', label='95% Confidence Band')
plt.plot(Ebetas, Ban2, 'r--')
plt.title('Inverse Gaussian Probability Plot')
plt.xlabel('Uniform (0,1) Quantiles')
plt.ylabel('F(xi; betaGum, aGum)')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.legend()
plt.show()


# Quantile functions for the distributions
FinvExp = -mumle * np.log(1 - Ebetas)
FinvGam = gamma.ppf(Ebetas, a=alphamle, scale=betagamle)
FinvNor = norm.ppf(Ebetas, loc=mumle, scale=sigmle)
FinvLN = lognorm.ppf(Ebetas, s=sigLN, scale=np.exp(muLN))
#La dist. inversa Weibull para los datos es x=exp(a+ln(-ln(1-u))/beta)
FinvWei = np.exp(amle + np.log(-np.log(1 - Ebetas)) / betaWei)
#La dist. inv. Gumbel de mínimos para los datos es x=a+ln(-ln(1-u))/beta
FinvGum = amleGum + np.log(-np.log(1 - Ebetas)) / betaGum
FinvIG = norm.ppf(Ebetas, loc=mumle, scale=np.sqrt(1 / lambmle))
miniX = xsort[0]
maxiX = xsort[-1]
# Quantile functions for confidence bands
Fban1 = -mumle * np.log(1 - Ban1)
Fban2 = -mumle * np.log(1 - Ban2)
Fban1Gam = gamma.ppf(Ban1, a=alphamle, scale=betagamle)
Fban2Gam = gamma.ppf(Ban2, a=alphamle, scale=betagamle)
Fban1Nor = norm.ppf(Ban1, loc=mumle, scale=sigmle)
Fban2Nor = norm.ppf(Ban2, loc=mumle, scale=sigmle)
Fban1LN = lognorm.ppf(Ban1, s=sigLN, scale=np.exp(muLN))
Fban2LN = lognorm.ppf(Ban2, s=sigLN, scale=np.exp(muLN))
Fban1Wei = np.exp(amle + np.log(-np.log(1 - Ban1)) / betaWei)
Fban2Wei = np.exp(amle + np.log(-np.log(1 - Ban2)) / betaWei)
Fban1Gum = amleGum + np.log(-np.log(1 - Ban1)) / betaGum
Fban2Gum = amleGum + np.log(-np.log(1 - Ban2)) / betaGum
Fban1IG = norm.ppf(Ban1, loc=mumle, scale=np.sqrt(1 / lambmle))
Fban2IG = norm.ppf(Ban2, loc=mumle, scale=np.sqrt(1 / lambmle))

#   QQ plots

#   QQ Exponencial sin censura:
plt.figure()
plt.scatter(FinvExp, xsort, s=5, label='Empirical')
plt.plot([miniX, maxiX], [miniX, maxiX], 'r--', label='45-degree line')
plt.title('Exponential QQ Plot')
plt.xlabel('Estimated Quantiles')
plt.ylabel('Ordered Data')
plt.legend()
plt.show()

#   QQ Gamma sin censura:
plt.figure()
plt.scatter(FinvGam, xsort, s=5, label='Empirical')
plt.plot([miniX, maxiX], [miniX, maxiX], 'r--', label='45-degree line')
plt.title('Gamma QQ Plot')
plt.xlabel('Estimated Quantiles')
plt.ylabel('Ordered Data')
plt.legend()
plt.show()

#   QQ Normal sin censura
plt.figure()
plt.scatter(FinvNor, xsort, s=5, label='Empirical')
plt.plot([miniX, maxiX], [miniX, maxiX], 'r--', label='45-degree line')
plt.title('Normal QQ Plot')
plt.xlabel('Estimated Quantiles')
plt.ylabel('Ordered Data')
plt.legend()
plt.show()

#   QQ LogNormal sin censura
plt.figure()
plt.scatter(FinvLN, xsort, s=5, label='Empirical')
plt.plot([miniX, maxiX], [miniX, maxiX], 'r--', label='45-degree line')
plt.title('Lognormal QQ Plot')
plt.xlabel('Estimated Quantiles')
plt.ylabel('Ordered Data')
plt.legend()
plt.show()

#   QQ Weibull sin censura
plt.figure()
plt.scatter(FinvWei, xsort, s=5, label='Empirical')
plt.plot([miniX, maxiX], [miniX, maxiX], 'r--', label='45-degree line')
plt.title('Weibull QQ Plot')
plt.xlabel('Estimated Quantiles')
plt.ylabel('Ordered Data')
plt.legend()
plt.show()

#   QQ Gumbell sin censura
plt.figure()
plt.scatter(FinvGum, xsort, s=5, label='Empirical')
plt.plot([miniX, maxiX], [miniX, maxiX], 'r--', label='45-degree line')
plt.title('Gumbel QQ Plot')
plt.xlabel('Estimated Quantiles')
plt.ylabel('Ordered Data')
plt.legend()
plt.show()

#   QQ Gaussiana Inversa sin censura
plt.figure()
plt.scatter(FinvIG, xsort, s=5, label='Empirical')
plt.plot([miniX, maxiX], [miniX, maxiX], 'r--', label='45-degree line')
plt.title('Inverse Gaussian QQ Plot')
plt.xlabel('Estimated Quantiles')
plt.ylabel('Ordered Data')
plt.legend()
plt.show()


# Evaluate densities on a grid
xs = np.linspace(0.01, xsort[-1] + 45, 200)

# Lognormal with censorship
ylnor = lognorm.pdf(xs, s=sigLN, scale=np.exp(muLN))

# Gamma
ygam = gamma.pdf(xs, a=alphamle, scale=betagamle)

# Weibull
aux = betaWei * (np.log(xs) - amle)
weidensi = (betaWei / xs) * np.exp(aux - np.exp(aux))

# Gumbel
auxG = betaGum * (xs - amleGum)
gumdensi = betaGum * np.exp(auxG - np.exp(auxG))

# Normal
ynor = norm.pdf(xs, mumle, sigmle)

# Exponential
yexp = np.exp(-xs / mumle) / mumle

# Inverse Gaussian
yGinv = invgauss.pdf(xs, mu=mumle, scale=np.sqrt(1 / lambmle))

# Maximum density value for plotting
densitop = max([max(ygam), max(ylnor), max(weidensi), max(gumdensi), max(ynor), max(yexp), max(yGinv)])


#Para datos máximos de Cruapán se grafican sólo la Gamma y Lognormal 
# pues fueron las dos mejores densidades para datos exactos:
plt.figure()
plt.plot(xs, ygam, label='Gamma (black)')
plt.plot(xs, ylnor, 'r--', label='LogNormal (red)')
plt.scatter(xsort, np.zeros_like(xsort), color='green', marker='o', label='Unique Data Points')
plt.title('Estimated Densities: Cruapán Data')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.xlim(xs[0], xs[-1])
plt.legend()
plt.show()

#Para los datos ratones Grice-Bain las mejores dist fueron Gumbel,Nor y Weibull
plt.figure()
plt.plot(xs, gumdensi, label='Gumbel (black)')
plt.plot(xs, ynor, 'r--', label='Normal (red)')
plt.plot(xs, weidensi, 'b-', label='Weibull (blue)')
plt.scatter(xsort, np.zeros_like(xsort), color='green', marker='o', label='Data Points')
plt.title('Estimated Densities: Grice-Bain Data')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.xlim(xs[0], xs[-1])
plt.legend()
plt.show()


#Para Ratones:
#Se calcula F(0):
FGumCero = 1 - np.exp(-np.exp(-betaGum * amleGum))
FNorCero = norm.cdf(0, loc=mumle, scale=sigmle)
FWeiCero = 0  # Siempre 0 para la Weibull

print("Para la Gumbel,  F(0) es: ", FGumCero)
print("Para la Normal,  F(0) es: ", FNorCero)
print("Para la Weibull, F(0) es: ", FWeiCero)

# EMVs with exact data
results = {
    "Normal mu y sigma": [mumle, sigmle],
    "LogNormal mu y sigma": [muLN, sigLN],
    "Exponential media": mumle,
    "Gamma forma, media y escala": [alphamle, mumle, betagamle],
    "Weibull forma beta y escala sigma": [betaWei, np.exp(amle)],
    "Gumbel minima, localización aGum, y escala bGum": [amleGum, deltamle],
    "Gaussiana Inversa media y forma": [mumle, lambmle]
}

max_length_results = 2 + max(len(key) for key in results.keys())

# Write results to a text file
with open("InferenciasExactas.txt", "w",encoding="utf-8") as file:
    file.write("EMVs con datos exactos:\n")
    for key, value in results.items():
        file.write(f"{key.ljust(max_length_results)}: {value}\n")

    file.write("\nCriterio de Información de Akaike:\n")
    for dist, aic in sorted_aic_dict.items():
        file.write(f"{dist.ljust(max_length)}: {aic}\n")

    file.write("\nEstadísticas de Anderson Darling:\n")
    for dist, ad in sorted_ads.items():
        file.write(f"{dist.ljust(max_length)}: {ad}\n")

    file.write("\nPromedios ponderados de distancias verticales:\n")
    for dist, wavd in sorted_wavds.items():
        file.write(f"{dist.ljust(max_length)}: {wavd}\n")
