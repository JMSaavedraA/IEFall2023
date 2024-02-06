
import numpy as np
from scipy.optimize import minimize_scalar, minimize, brentq

##  ............................................................................
##      Funciones usadas:                                                   ####

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
###   LogVer(mu,sigma)                                             ####

#   Se calcula la logverosimilitud normal en [mu,sigma].
#   Los argumentos son el vector de dos parámetros, las estadísticas 
#   suficientes t1,t2 y el tamaño de muestra n.

def lvnor(vec2, t1, t2, n):
    mu = vec2[0]
    sigma = vec2[1]
    if sigma > 0:
        s2 = sigma**2
        lv = -n * np.log(sigma) - (t2 / (2 * s2)) + (mu * t1 / s2) - ((n * mu**2) / (2 * s2))
    else:
        lv = -99999999999999999
    return lv

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
###  Cálculo de lp(mu) numérico:                                           ####
#   lpmu es idéntica que l(mu,sigma) cambiando argumentos al reducir en uno el 
#   vector parámetros. Así sólo es función de sigma, parámetro de estorbo, de
#   las estadísticas suficientes, de n y del valor fijo de MU en mufix.

def lpmu(vec1, t1, t2, n, mufix):
    mu = mufix
    sigma = vec1
    if sigma > 0:
        s2 = sigma**2
        lv = -n * np.log(sigma) - (t2 - (2 * t1 * mu) + (n * (mu**2))) / (2 * s2)
    else:
        lv = -99999999999999999
    return lv

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
###   rp(mu) analítica                               ####

def rpmu(mu, t1, t2, n):
    HH = t2 - (t1**2) / n
    muhat = t1 / n
    Hmu = t2 - (2 * t1 * mu) + (n * (mu**2))
    rp = -(n / 2) * (np.log(Hmu) - np.log(HH))
    return rp

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
###   rp(sigma)-lnc para calcular interv perfil sigma numéricos ####

def rpsigC(sig, SIGMV, NN, CC):
    lnC = np.log(CC)
    if sig > 0:
        rpsig = NN * np.log(SIGMV / sig) - (NN / 2) * ((SIGMV / sig)**2) + (NN / 2)
        rpc = rpsig - lnC
    else:
        rpc = 999
    return rpc

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
###       Niveles Verosimilitud: c(n) para n>=5                      ####
#  se calcula el nivel de verosimilitud como función de n
#  para intervalos de verosimilitud perfil de mu y de sigma normales
#  con confianzas del 90, 95 y 99%

def clevels(n):
    cc90 = 0.2585 - (1 / (-0.4743 + (1.9029 * n)))
    cc95 = 0.1465 - (1 / (0.8469 + (2.3544 * n)))
    cc99 = 0.0362 - (1 / (8.753 + (5.551 * n)))
    ccesmu = {"90%": cc90, "95%": cc95, "99%": cc99}
    
    cs90 = 0.2585 - (1 / (0.775 + (1.641 * n)))
    cs95 = 0.1465 - (1 / (2.757 + (2.082 * n)))
    cs99 = 0.0362 - (1 / (18.294 + (5.229 * n)))
    ccessig = {"90%": cs90, "95%": cs95, "99%": cs99}

    ctous90 = 0.10 - (1 / (1.11 + (4.655 * n)))
    ctous95 = 0.05 - (1 / (1.61 + (7.27 * n)))
    ctous99 = 0.01 - (1 / (19.51 + (23.59 * n)))
    ctous = {"90%": ctous90, "95%": ctous95, "99%": ctous99}
    
    cesn = {"cmu": ccesmu, "csig": ccessig, "cmusig": ctous}
    return cesn

def find_root(func, lower, upper):
    return brentq(func, lower, upper)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  . . ..
###   Intervalos Perfil Sigma numéricos                       #### 
#Con confianza del 90, 95 o 99%
#--------------------------- --

def IntSig(SIGMV, NN, confi):
    cpp1s = [0.9, 0.95, 0.99]
    
    if confi == 0.90:
        cp = cpp1s[0]
    elif confi == 0.95:
        cp = cpp1s[1]
    elif confi == 0.99:
        cp = cpp1s[2]
    else:
        print("Mala confianza. Se calculará intervalo sigma del 95%:")
        cp = cpp1s[1]

    s11 = find_root(lambda sigma: rpsigC(sigma, SIGMV, NN, cp), 0.00001, SIGMV)
    s12 = find_root(lambda sigma: rpsigC(sigma, SIGMV, NN, cp), SIGMV, SIGMV * 10)

    intersig = {"Sig1": s11, "Sig2": s12, "Confianza": confi}
    return intersig

